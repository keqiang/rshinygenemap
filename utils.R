`%>%` <- magrittr::`%>%`

C_MAX_NUMBER_OF_GENES <- 200

speciesMappings <- c("Human" = "hs", "Mouse" = "mm", "Rat" = "rn", "Zebrafish" = "dr", "Fruitfly" = "dm")
geneIdTypeMappings <- c("Ensembl Gene ID" = "ensemblgid", "NCBI Gene ID" = "ncbigid", "Gene Symbol" = "symbol")

append_gene_symbols <- function(species, gene_id_type = geneIdTypeMappings, gene_ids, gene_symbols = NULL) {
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

geneIdTypeReverseMappings <- names(geneIdTypeMappings) %>% rlang::set_names(geneIdTypeMappings)

getGeneIdTypeName <- function(geneIdTypeShortName = geneIdTypeMappings) {
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

add_ncbigid_link <- function(species, ncbi_id) {
  ensembl_url <- "https://www.ncbi.nlm.nih.gov/gene"
  wrap_with_link(ensembl_url, list("term" = ncbi_id), ncbi_id)
}
