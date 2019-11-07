library(shiny)

`%then%` <- shiny:::`%OR%`
`%>%` <- magrittr::`%>%`

speciesMappings <- c("Human" = "hs", "Mouse" = "mm", "Rat" = "rn", "Zebrafish" = "dr", "Fruitfly" = "dm")
geneIdTypes <- c("Ensembl Gene ID" = "ensemblgid", "NCBI Gene ID" = "ncbigid", "Gene Symbol" = "symbol")

append_gene_symbols <- function(species, gene_id_type = geneIdTypes, gene_ids, gene_symbols = NULL) {
  if (gene_id_type == "symbol" || is.null(gene_symbols)) {
    gene_symbols <- gene_ids # ignore gene_symbols when gene_ids is already symbols
  }
  with_harmonizome_link <- ifelse(species == "hs", TRUE, FALSE) # harmonizome is only for human genes
  gene_id_type <- match.arg(gene_id_type)
  purrr::map2_chr(
    gene_ids,
    gene_symbols,
    function(id, symbol) {
      if (gene_id_type == "symbol") { # id is gene symbol already, only need to add harmonizome link when it's human genes
        ifelse(with_harmonizome_link && !is.na(id), add_harmonizome_link(id), id)
      } else { # ensembl or ncbi gid, add link
        id_with_link <- ifelse(!is.na(id), do.call(glue::glue("add_{gene_id_type}_link"), list(species, id)), id)
        ifelse(
          is.na(symbol), # if there is no symbol, return id
          id_with_link,
          ifelse(
            with_harmonizome_link, # if is human gene symbol, add harmonizome link
            glue::glue("{id_with_link} ({add_harmonizome_link(symbol)})"),
            glue::glue("{id_with_link} ({symbol})")
          )
        )
      }
    }
  )
}

speciesReverseMapping <- names(speciesMappings) %>% rlang::set_names(speciesMappings)

getSpeciesName <- function(speciesShortName = speciesMappings) {
  speciesShortName <- match.arg(speciesShortName)
  speciesReverseMapping[[speciesShortName]]
}

geneIdTypeReverseMappings <- names(geneIdTypes) %>% rlang::set_names(geneIdTypes)

getGeneIdTypeName <- function(geneIdTypeShortName = geneIdTypes) {
  geneIdTypeShortName <- match.arg(geneIdTypeShortName)
  geneIdTypeReverseMappings[[geneIdTypeShortName]]
}

wrap_with_link <- function(base_url, params = NULL, text = "", color = NULL) {
  href_url <- ifelse(
    is.null(params),
    base_url,
    glue::glue(
      base_url,
      "?",
      names(params) %>% purrr::reduce(function(accu, cur) {
        cur_param <- glue::glue(cur, "=", params[[cur]])
        ifelse(is.null(accu), cur_param, glue::glue(accu, "&", cur_param))
      }, .init = NULL)
    )
  )
  tags$a(
    target = "_blank",
    style = ifelse(is.null(color), "", glue::glue("color: {color};")),
    href = href_url,
    text
  ) %>% as.character()
}

add_harmonizome_link <- function(gene) {
  harmonizome_url <- glue::glue("http://amp.pharm.mssm.edu/Harmonizome/gene/", gene)
  wrap_with_link(harmonizome_url, text = gene, color = "red")
}

add_ensemblgid_link <- function(species, ensembl_id) {
  path_mappings <- list(
    "hs" = "Homo_sapiens",
    "mm" = "Mus_musculus",
    "rn" = "Rattus_norvegicus",
    "dr" = "Danio_rerio",
    "dm" = "Drosophila_melanogaster"
  )
  ensembl_url <- glue::glue("https://useast.ensembl.org/{path_mappings[species]}/Gene/Summary")
  wrap_with_link(ensembl_url, list("g" = ensembl_id), ensembl_id)
}

ui <- fluidPage(
  titlePanel("Gene identifiers mapping within- or cross-species"),
  sidebarLayout(
    sidebarPanel(
      width = 6,
      tags$h5("Step 1: import a data table file in which the first column is gene IDs to map"),
      shinywidgets::fileImportWidget(
        id = "inputGeneTable",
        shinywidgets::C_DATA_TYPE_TABLE,
        enableDataTypeSelection = FALSE
      ),
      tags$h5("Step 2: define the input genes"),
      selectInput(
        "inputSpecies",
        label = "Species",
        choices = speciesMappings,
        selected = "hs"
      ),
      selectInput(
        "inputIdType",
        label = "Gene ID type",
        choices = geneIdTypes,
        selected = "ensemblgid"
      ),
      tags$h5("Define the output genes"),
      selectInput(
        "outputSpecies",
        label = "Output species",
        choices = speciesMappings,
        selected = "hs"
      ),
      selectInput(
        "outputIdType",
        label = "Output gene id type",
        choices = geneIdTypes,
        selected = "symbol"
      ),
      actionButton(
        "mapGenes",
        label = "Map Genes"
      ),
      textOutput("errorMsg")
    ),
    mainPanel(
      width = 6,
      DT::dataTableOutput("mappedGeneTable"),
      checkboxInput(
        "onlyShowMappedGenes",
        label = "Only show genes that have been successfully mapped",
        value = TRUE
      ),
      textOutput("mappingHint"),
      tags$hr(),
      checkboxInput(
        "onlyPreserveMapped",
        label = "Preserve only genes that have mappings in the result",
        value = TRUE
      )
    )
  )
)

server <- function(input, output, session) {
  importedGeneTable <- shinywidgets::importFile("inputGeneTable", shinywidgets::C_FILE_LOCATION_LOCAL)

  inputGeneTable <- reactive({
    tableData <- importedGeneTable()
    req(tableData)
    tableData$data
  })

  output$errorMsg <- renderText({
    shinyjs::disable("mapGenes")

    fromSpecies <- input$inputSpecies
    fromType <- input$inputIdType
    toSpecies <- input$outputSpecies
    toType <- input$outputIdType

    validate(
      need(importedGeneTable(), "Please select input data") %then%
        need(fromSpecies != toSpecies || fromType != toType, "No need to map genes of same species and ID types")
    )
    shinyjs::enable("mapGenes")
  })

  geneMappingRes <- eventReactive(input$mapGenes, {
    withProgress({
      fromSpecies <- input$inputSpecies
      fromType <- input$inputIdType
      toSpecies <- input$outputSpecies
      toType <- input$outputIdType
      genesToMap <- inputGeneTable()[[1]]
      res <- genemap::map_genes(fromSpecies, genesToMap, fromType, toSpecies, toType)
      req(res)
      res
    }, message = "Mapping genes. Please wait...")
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
    res <- bind_cols(
      mappedGeneIdTable(),
      isolate(inputGeneTable())
    )
    if (input$onlyPreserveMapped) {
      res <- res %>% filter(!is.na(.[[1]]))
    }
    res
  })

  originalIdReplaced <- reactive({
    res <- bind_cols(
      mappedGeneIdTable() %>% select(1),
      inputGeneTable() %>% select(-1)
    )
    if (input$onlyPreserveMapped) {
      res <- res %>% filter(!is.na(.[[1]]))
    }
    res
  })

  dataObjectsCanBeExported <- reactive({
    list(
      newExportableDataObject(
        mappedResult(),
        type = C_DATA_TYPE_TABLE,
        name = "mapped_result_with_original_ids"
      ),
      newExportableDataObject(
        originalIdReplaced(),
        type = C_DATA_TYPE_TABLE,
        name = "mapped_result_original_ids_replaced"
      )
    )
  })
}

shinyApp(ui, server)
