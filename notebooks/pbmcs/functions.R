library(tidyverse)
library(purrr)
prep_pbmc_data <- function(){
  pbmc_data <- new.env()
  load('data/pbmc-purified/deseq2-pbmc-purified.RData', envir = pbmc_data)

  #' add columns for gene name and cell type to the deseq results
  add_celltype <- function(celltype){
    pbmc_data$deseq[[celltype]] %>%
      as.data.frame() %>%
      rownames_to_column(var = 'gene') %>%
      mutate(celltype = celltype)
  }

  # make big de table
  de <- map_dfr(names(pbmc_data$deseq), add_celltype)

  # gene label mapping
  hs <- org.Hs.eg.db::org.Hs.eg.db
  gene_symbols <- unique(de$gene)
  symbol2entrez <- AnnotationDbi::select(
    hs, keys=gene_symbols,
    columns=c('ENTREZID', 'ENSEMBL'),
    keytype = 'ENSEMBL')

  #' helper function to assign names to a list in a tidy pipeline
  add_names = function(l, n){
    names(l) <- n
    return(l)
  }

  # organize data into a list of tibbles, one for each cell type
  data <- de %>%
    dplyr::rename('ENSEMBL' = gene) %>%
    left_join(symbol2entrez, by='ENSEMBL') %>%
    relocate(ENTREZID, .after=ENSEMBL) %>%
    mutate(  # set default columns
      beta = log2FoldChange,
      se = lfcSE,
      threshold.on = padj
    ) %>%
    group_by(celltype) %>%
    group_map(~ .x, .keep = T) %>%
    add_names(map_chr(., ~pluck(.x, 'celltype')[1]))
  return(data)
}

#' Construct gene list as all genes with pvalue less than quantile `q`
make_gene_list <- function(data, celltype, q){
  de <- data[[celltype]] %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::filter(!duplicated(ENTREZID))

  # background is all the genes with a pvalue reported in this cell type
  # (so it is the genes that *might have been* included in the gene list)
  gene_list <- de %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::filter(pvalue < quantile(pvalue, q, na.rm=T)) %>%
    {.$ENTREZID}
  gene_background <- de$ENTREZID
  return(list(gene_list=gene_list, gene_background=gene_background))
}

# Function to create the directory from a file path if it does not exist
create_directory_from_path <- function(file_path) {
  # Extract the directory path from the given file path
  dir_path <- dirname(file_path)

  # Check if the directory exists; if not, create it
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

save_model <- function(fit, path){
  fit$data$X <- NULL
  create_directory_from_path(path)
  saveRDS(fit, path)
}


########
# INTERACTIVE TABLE
########

library(reactable)

make_get_details <- function(tbl){
  if('details' %in% names(tbl)){
    get_details <- function(index){
      deets <- tbl[index,]
      htmltools::div(style = "padding: 1rem",
                     reactable(
                       dplyr::select(deets$details[[1]], -any_of('details')),
                       details = make_get_details(deets$details[[1]]),
                       outlined = TRUE,
                       defaultColDef = colDef(format=colFormat(digits=3))
                     )
      )
    }
  } else{
    get_details <- NULL
  }
  return(get_details)
}

#' Make a drillable html table
#'
#' @param cs_tbl_nested is a nested table, where there is a column named details. The produced html table will display all columns except details, and show details associated with a particular row when the user selects the row. The detail column may also be nested.
#' @returns an html table with the option to drill down each row, potentially multiple layers of nesting.
make_cs_reactable <- function(cs_tbl_nested){
  reactable(
    dplyr::select(cs_tbl_nested, -any_of('details')),
    details = make_get_details(cs_tbl_nested),
    defaultColDef = colDef(format=colFormat(digits=3))
  )
}

