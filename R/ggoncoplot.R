
#' GG oncoplot
#'
#' @param col_genes name of \strong{data_main} column containing gene names/symbols (string)
#' @param col_samples name of \strong{data_main} column containing sample identifiers (string)
#' @param col_mutation_type name of \strong{data_main} column describing mutation types (string)
#' @param col_tooltip name of \strong{data_main} column containing whatever information you want to display in (string)
#' @param topn how many of the top genes to visualise. Ignored if \code{genes_to_include} is supplied (number)
#' @param show_sample_ids show sample_ids_on_x_axis (flag)
#' @param data_main data for main oncoplot (data.frame)
#'
#' @return ggplot or ggiraph object if \code{interactive=TRUE}
#' @export
#'
#' @examples
#' GBM
#' gbm_csv <- system.file(package='oncoplotgg', "testdata/GBM_tcgamutations_mc3_maf.csv")
#' gbm_df <- read.csv(file = gbm_csv, header=TRUE)
#' ggoncoplot(gbm_df, 'Hugo_Symbol', 'Tumor_Sample_Barcode', col_mutation_type = 'Variant_Classification')
#'
#' # BRCA
#' brca_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_maf.csv")
#' brca_df <- read.csv(file = brca_csv, header=TRUE)
#' ggoncoplot(brca_df, 'Hugo_Symbol', 'Tumor_Sample_Barcode', col_mutation_type = 'Variant_Classification')
ggoncoplot <- function(data_main, col_genes, col_samples, col_mutation_type = NULL, col_tooltip = col_samples, topn=10, show_sample_ids=FALSE, interactive=TRUE, interactive_svg_width = 10, interactive_svg_height = 6, genes_to_include = NULL, xlab_title="Sample", ylab_title="Gene", sample_annotation_df = NULL){
  assertthat::assert_that(is.data.frame(data_main))
  assertthat::assert_that(assertthat::is.string(col_genes))
  assertthat::assert_that(assertthat::is.string(col_samples))
  assertthat::assert_that(is.null(col_mutation_type) | assertthat::is.string(col_mutation_type))
  assertthat::assert_that(is.null(genes_to_include) | is.character(genes_to_include))
  assertthat::assert_that(assertthat::is.string(col_tooltip))
  assertthat::assert_that(assertthat::is.number(topn))



  # Check specified columns are in data_main
  data_main_colnames <- names(data_main)

  check_valid_dataframe_column(
    data = data_main,
    colnames = c(
      col_samples,
      col_genes,
      col_tooltip
    )
  )

  # Check optional columns are in data_main
  if (!is.null(col_mutation_type))
    check_valid_dataframe_column(data = data_main, colnames = col_mutation_type)

  # Ensure Sample Column is A factor
  data_main[[col_samples]] <- as.factor(data_main[[col_samples]])

  # Look exclusively at a custom set of genes
  if(!is.null(genes_to_include)){
    genes_not_found <- genes_to_include[!genes_to_include %in% data_main[[col_genes]]]

    if(length(genes_not_found) > 0){
      cli::cli_alert_warning("Failed to find the following [{length(genes_not_found)}] genes in your dataset")

      cli::cli_alert("{genes_not_found}")
      cli::cli_alert_warning("Either no samples have mutations in the above genes, or you've got the wrong gene names")
    }

    if(length(genes_not_found) == length(genes_to_include))
      cli::cli_abort("Couldn't find any of the genes you supplied in your dataset. Either no samples have mutations in these genes, or you've got the wrong gene names")

    genes_in_selected_set <- data_main[[col_genes]] %in% genes_to_include

    data_main <- data_main |>
     dplyr::filter(.data[[col_genes]] %in% genes_to_include)
  }

  # Filter out genes_to_ignore
  ## I'll have to implement later :)


  # Identify top genes by frequency
  data_gene_counts <- data_main |>
    dplyr::distinct(.data[[col_samples]], .data[[col_genes]], .keep_all = TRUE) |> # This line stops multiple mutations of the same gene in the same sample counting multiple times towards the mutation frequency.
    dplyr::count(.data[[col_genes]])

  # If use hasn't specifically specified genes to include, filter for the 'topn' most mutated genes
  if(!is.null(genes_to_include))
    data_main_top_genes_df  <- data_gene_counts
  else{
    data_main_top_genes_df  <- data_gene_counts |>
      dplyr::slice_max(.data$n, n = topn)
  }

  data_main_top_genes <- data_main_top_genes_df |>
    dplyr::pull(.data[[col_genes]])

  data_main_top_genes_rank <- rank(data_main_top_genes_df[['n']], ties.method = "first")

  # Filter dataset to only include the topn genes
  data_main_top_df <- data_main |>
    dplyr::filter(.data[[col_genes]] %in% data_main_top_genes)

  # Order Genes Variable
  #browser()
  data_main_top_df[[col_genes]] <- forcats::fct_relevel(data_main_top_df[[col_genes]], data_main_top_genes[order(data_main_top_genes_rank)])

  # Sort Samples by mutated gene
  data_main_top_df <- data_main_top_df |>
    dplyr::group_by(.data[[col_samples]])|>
    dplyr::mutate(
      SampleRankScore = score_based_on_gene_rank(mutated_genes = .data[[col_genes]], genes_informing_score = data_main_top_genes, gene_rank = data_main_top_genes_rank) # add secondary ranking based on secondary
    ) |>
    dplyr::ungroup()
  data_main_top_df[[col_samples]] <- forcats::fct_rev(forcats::fct_reorder(data_main_top_df[[col_samples]], data_main_top_df$SampleRankScore))

  # Colour Samples by mutation type
  if(!is.null(col_mutation_type)){
    data_main_top_df <- data_main_top_df |>
      dplyr::group_by(.data[[col_samples]], .data[[col_genes]]) |>
      dplyr::mutate(MutationType = ifelse(
          dplyr::n_distinct(.data[[col_mutation_type]]) > 1,
          yes = "Multiple",
          no = unique(.data[[col_mutation_type]])
        ) |>
        forcats::fct_infreq(f = _)
      )
    col_mutation_type = "MutationType"
  }

  # Create ggplot
  gg <- ggplot2::ggplot(data = data_main_top_df, mapping = ggplot2::aes_string(
    y = col_genes,
    x = col_samples,
    fill = col_mutation_type
  ))

  # Add interactive/non-interactive geom layer
  if(interactive){
    gg <- gg + ggiraph::geom_tile_interactive(
      ggplot2::aes(
        tooltip = .data[[col_tooltip]],
        data_id = .data[[col_samples]])
      )
  }
  else {
    gg <- gg + ggplot2::geom_tile()
  }

  # Apply default theme
  gg <- gg + theme_oncoplot_default()

  # Label axis
  gg <- gg + ggplot2::xlab(xlab_title) + ggplot2::ylab(ylab_title)

  # Add fill colour
  gg <- gg + ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(n = 12, "Paired"))

 # Show/hide sample ids on x axis
  if(!show_sample_ids){
    gg = gg + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # Turn gg into an interactive ggiraph object if interactive = TRUE
  if(interactive){
    gg <- ggiraph::girafe(
      width_svg=interactive_svg_width, height_svg = interactive_svg_height,
      ggobj = gg,
      options = list(
        ggiraph::opts_tooltip(
          opacity = .8,
          css = "background-color:gray;color:white;padding:2px;border-radius:2px;"
        )
      ))
  }
  return(gg)
}

make_interactive <- function(gg){
  ggiraph::girafe(
    width_svg=10, height_svg = 6,
    ggobj = gg,
    options = list(
      ggiraph::opts_tooltip(
        opacity = .8,
        css = "background-color:gray;color:white;padding:2px;border-radius:2px;"
      )
    ))
}
#' Generate score based on genes
#'
#' Score used to sort samples based on which genes are mutated. Make sure to run on one sample at once (use grouping)
#'
#' @param mutated_genes vector of genes that are mutated for a single sample (character)
#' @param genes_informing_score which genes determine the sort order? (character)
#' @param gene_rank what is the order of importance of genes used to determine sort order. Higher number = higher in sort order (character)
#' @param debug_mode debug mode (flag)
#'
#' @return a score (higher = should be higher in the sorting order) (number)
#'
#' @examples
#' \dontrun{
#' # First set of genes has a high rank since both BRCA2 and EGFR are mutated
#' score_based_on_gene_rank(c('TERT', 'EGFR', 'PTEN', 'BRCA2'), c('EGFR', 'BRCA2'), gene_rank = 1:2)
#'
#' # If EGFR is mutated without BRCA2, we get a lower score
#' score_based_on_gene_rank(c('TERT', 'EGFR', 'PTEN', 'IDH1'), c('EGFR', 'BRCA2'), gene_rank = 1:2)
#'
#' # If BRCA2 is mutated without EGFR,
#' # we get a score lower than BRCA2+EGFR but higher than EGFR alone due to higher gene_rank of BRCA2
#' score_based_on_gene_rank(c('TERT', 'IDH1', 'PTEN', 'BRCA2'), c('EGFR', 'BRCA2'), gene_rank = 1:2)
#' }
score_based_on_gene_rank = function(mutated_genes, genes_informing_score, gene_rank, debug_mode = FALSE){
  assertthat::assert_that(is.character(mutated_genes) | is.factor(mutated_genes))
  assertthat::assert_that(is.character(genes_informing_score))
  assertthat::assert_that(is.numeric(gene_rank))
  assertthat::assert_that(length(genes_informing_score)==length(gene_rank))


  gene_order = rank(gene_rank, ties.method = "first")


  res = genes_informing_score %in% mutated_genes
  names(res) <- genes_informing_score
  base10_values = ifelse(res, yes = 2^(gene_order-1), no=0)

  total_score <- sum(base10_values)

  if(debug_mode){
    cli::cli_alert("base_values = [{names(res)}].")
    cli::cli_alert("base_values = [{base10_values}].")
    cli::cli_alert("total score = [{total_score}].")
  }
  return(total_score)
}


#' Oncoplot Theme: default
#'
#' @param ... passed to [ggplot2::theme()] theme
#' @importFrom ggplot2 %+replace%
theme_oncoplot_default <- function(...){
  ggplot2::theme_grey(...) %+replace%
    ggplot2::theme(
      panel.border = ggplot2::element_rect(size = 1, fill = NA),
      panel.grid = ggplot2::element_line(colour = "white"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}

check_valid_dataframe_column <- function(data, colnames, error_call = rlang::caller_env()){
  data_colnames = colnames(data)

  for (colname in colnames){
    if(!colname %in% data_colnames){
      cli::cli_abort(
        c(
          '!' = 'Could find column: [{colname}] in input data',
          '!' = 'Please select one of [{data_colnames}]'
        )
      )
    }
  }

}
# devtools::load_all();brca_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_maf.csv")
# brca_df <- read.csv(file = brca_csv, header=TRUE)
# ggoncoplot(brca_df, 'Hugo_Symbol', 'Tumor_Sample_Barcode', col_mutation_type = 'Variant_Classification')
