#' Expression Heatmap
#'
#' @param data_main
#' @param col_genes
#' @param col_samples
#' @param col_expression
#' @param genes_to_include
#'
#' @return
#' @export
#'
#' @examples
gg_expression_heatmap <- function(data_main, col_genes, col_samples, col_expression, genes_to_include = NULL){
  assertthat::assert_that(is.data.frame(data_main))
  assertthat::assert_that(assertthat::is.string(col_genes))
  assertthat::assert_that(assertthat::is.string(col_samples))
  assertthat::assert_that(assertthat::is.string(col_expression))
  assertthat::assert_that(is.null(genes_to_include) | is.character(genes_to_include))



}


gg_expression_oncoplot <- function(data_main, col_genes, col_samples){

}



# brca_rnaseq_path = system.file(package = "oncoplotgg","testdata/BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.data.txt")
# brca_rnaseq_df <- data.table::fread(brca_rnaseq, header = TRUE, sep = '\t')
