library(shiny)
source("utils.R") # source in the utility functions

ui <- fluidPage(
  theme = shinythemes::shinytheme("sandstone"),
  titlePanel(
    "Gene identifiers mapping within- or cross-species",
    windowTitle = "GeneMap"
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$h4("Import a data table"),
      helpText(
        glue::glue(
          "Choose a file in which the first column has gene IDs to map.",
          " Will only match the first {C_MAX_NUMBER_OF_GENES} genes if the file has more than that."
        )
      ),
      shinyfio::dataTableImportWidget(
        id = "inputGeneTable",
        label = "You can use the demo file on our server or import your own from your computer"
      )
    ),
    mainPanel(
      width = 8,
      conditionalPanel(
        condition = "output['isDataAvailable']",
        fluidRow(
          column(
            width = 6,
            wellPanel(
              tags$h4("Define the input genes"),
              selectInput(
                "inputSpecies",
                label = "Species",
                choices = speciesMappings,
                selected = "rn"
              ),
              selectInput(
                "inputIdType",
                label = "Gene ID type",
                choices = geneIdTypeMappings,
                selected = "ensemblgid"
              )
            )
          ),
          column(
            width = 6,
            wellPanel(
              tags$h4("Define the output genes"),
              selectInput(
                "outputSpecies",
                label = "Species",
                choices = speciesMappings,
                selected = "hs"
              ),
              selectInput(
                "outputIdType",
                label = "Gene ID type",
                choices = geneIdTypeMappings,
                selected = "symbol"
              )
            )
          )
        ),
        actionButton(
          "mapGenes",
          label = "Map Genes"
        ),
        textOutput("errorMsg"),
        tags$hr()
      ),
      conditionalPanel(
        condition = "output['isResultAvailable']",
        wellPanel(
          tags$h4("Mapping result"),
          checkboxInput(
            "onlyShowMappedGenes",
            label = "Only display genes that have been successfully mapped",
            value = TRUE
          ),
          DT::dataTableOutput("mappedGeneTable"),
          textOutput("mappingHint"),
          checkboxInput(
            "onlyPreserveMapped",
            label = "Only keep genes that have been successfully mapped in the downloaded files",
            value = TRUE
          ),
          downloadButton("downloadMappedTable", "Download mapped table"),
          downloadButton("downloadIdReplacedTable", "Download mapped table with mapped IDs only")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  importedGeneTable <- shinyfio::importDataTable(
    id = "inputGeneTable",
    maxNumberOfLines = C_MAX_NUMBER_OF_GENES,
    serverRootDirectories = c("sample_data" = file.path(".", "sample_data"))
  )

  inputGeneTable <- reactive({
    tableData <- importedGeneTable()
    req(tableData)
    tableData$data
  })

  output$isDataAvailable <- reactive({
    req(inputGeneTable())
    TRUE
  })
  outputOptions(output, "isDataAvailable", suspendWhenHidden = FALSE)

  errorMsg <- reactiveVal(NULL)
  output$errorMsg <- renderText({
    errorMsg()
  })

  observe({
    errorMsg(NULL)
    if (input$inputSpecies == input$outputSpecies && input$inputIdType == input$outputIdType) {
      errorMsg("No need to map genes of same species and ID types")
    }
  })

  geneMappingRes <- eventReactive(input$mapGenes, {
    fromSpecies <- input$inputSpecies
    fromType <- input$inputIdType
    toSpecies <- input$outputSpecies
    toType <- input$outputIdType

    req(is.null(errorMsg()))
    withProgress(
      message = "Mapping genes. Please wait...", {
        genesToMap <- inputGeneTable()[[1]]
        res <- genemap::map_genes(fromSpecies, genesToMap, fromType, toSpecies, toType)
        req(res)
        res
      }
    )
  })

  availableGeneMappings <- reactive({
    mappedGenes <- geneMappingRes()
    req(mappedGenes)
    if (input$onlyShowMappedGenes) {
      mappedGenes <- mappedGenes[mappedGenes %>% purrr::map_lgl(~ length(.$mapped_genes) > 0)]
    }
    mappedGenes
  })

  output$mappingHint <- renderText({
    mapped <- availableGeneMappings()
    req(mapped)

    fromSpeciesName <- isolate(getSpeciesName(input$inputSpecies))
    fromTypeName <- isolate(getGeneIdTypeName(input$inputIdType))
    toSpeciesName <- isolate(getSpeciesName(input$outputSpecies))
    toTypeName <- isolate(getGeneIdTypeName(input$outputIdType))

    numOfMappedGenes <- mapped %>% purrr::reduce(function(accu, cur) {
      accu <- accu + ifelse(length(cur$mapped_genes) > 0, 1, 0)
    }, .init = 0)

    numOfInputGenes <- isolate(nrow(inputGeneTable()))

    numOfGenesHasMapping <- length(mapped)

    validate(
      need(
        FALSE,
        glue::glue("{numOfGenesHasMapping} out of {numOfInputGenes} {fromSpeciesName} {fromTypeName} genes mapped to {numOfMappedGenes} {toSpeciesName} {toTypeName} genes")
      )
    )
  })

  geneMappingTable <- reactive({ # this table is only for data table to show (is dirty)
    mappedGenes <- availableGeneMappings()
    req(mappedGenes)

    fromSpecies <- isolate(input$inputSpecies)
    fromType <- isolate(input$inputIdType)

    toSpecies <- isolate(input$outputSpecies)
    toType <- isolate(input$outputIdType)

    addedGeneSymbols <- mappedGenes %>%
      purrr::map(~ append_gene_symbols(toSpecies, toType, .$mapped_genes, .$mapped_gene_symbols))

    list(
      append_gene_symbols(
        fromSpecies,
        fromType,
        names(mappedGenes),
        mappedGenes %>% purrr::map_chr(~ ifelse(is.null(.$gene_symbol), NA, .$gene_symbol)) # if gene has no mapping
      ),
      addedGeneSymbols %>% purrr::map_chr(~ .[1])
    ) %>%
      rlang::set_names(
        c(
          glue::glue("original_{fromSpecies}_gene_{fromType}"),
          glue::glue("mapped_{toSpecies}_gene_{toType}")
        )
      ) %>%
      tibble::as_tibble()
  })

  output$mappedGeneTable <- DT::renderDataTable({
    DT::datatable(
      geneMappingTable(),
      options = list(scrollX = TRUE),
      selection = "single",
      escape = FALSE
    )
  })

  mappedGeneIdTable <- reactive({
    mappedGenes <- geneMappingRes()
    toSpecies <- isolate(input$outputSpecies)
    toType <- isolate(input$outputIdType)
    res <- list(
      mappedGenes %>%
        purrr::map_chr(~ ifelse(length(.$mapped_genes) > 0, .$mapped_genes[[1]], NA))
    ) %>% rlang::set_names(glue::glue("mapped_{toSpecies}_gene_{toType}"))

    if (toType != "symbol") { # if not converting to symbol, then add symbol column
      res[[glue::glue("mapped_{toSpecies}_gene_symbol")]] <- mappedGenes %>%
        purrr::map_chr(~ ifelse(length(.$mapped_gene_symbols) > 0, .$mapped_gene_symbols[[1]], NA))
    } else {
      res <- rlang::set_names(res, glue::glue("mapped_{toSpecies}_gene_symbol"))
    }
    tibble::as_tibble(res)
  })

  mappedResult <- reactive({
    res <- dplyr::bind_cols(
      mappedGeneIdTable(),
      isolate(inputGeneTable())
    )
    if (input$onlyPreserveMapped) {
      res <- res %>% dplyr::filter(!is.na(.[[1]]))
    }
    res
  })

  originalIdReplaced <- reactive({
    res <- dplyr::bind_cols(
      mappedGeneIdTable() %>% dplyr::select(1),
      inputGeneTable() %>% dplyr::select(-1)
    )
    if (input$onlyPreserveMapped) {
      res <- res %>% dplyr::filter(!is.na(.[[1]]))
    }
    res
  })

  output$isResultAvailable <- reactive({
    req(mappedResult(), originalIdReplaced())
    TRUE
  })
  outputOptions(output, "isResultAvailable", suspendWhenHidden = FALSE)

  output$downloadMappedTable <- downloadHandler(
    filename = "mapped_result_with_original_ids.txt",
    content = function(file) {
      readr::write_tsv(mappedResult(), file)
    }
  )

  output$downloadIdReplacedTable <- downloadHandler(
    filename = "mapped_result_original_ids_replaced.txt",
    content = function(file) {
      readr::write_tsv(originalIdReplaced(), file)
    }
  )
}

shinyApp(ui, server)
