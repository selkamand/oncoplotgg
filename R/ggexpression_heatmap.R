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


#' Title
#'
#' @param data
#' @param col_sample
#'
#' @return
#' @export
#'
#' @examples
#' brca_metadata_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_clinicaldata.csv")
#' brca_metadata_df <- read.csv(file = brca_metadata_csv, header=TRUE)
#' gg_sample_annotation(brca_metadata_df, 'Tumor_Sample_Barcode')
gg_sample_annotation <- function(data, col_sample, col_sort, interactive = TRUE){

  #assertions
  assertthat::assert_that(assertthat::is.string(col_sample))
  assertthat::assert_that(is.data.frame(data))
  check_valid_dataframe_column(data, col_sample)
  check_valid_dataframe_column(data, col_sort)

  # Make sure there are enough columns in dataframe supplied
  ncols <- ncol(data)
  if (ncols <= 1)
    cli::cli_abort('Sample annotation data must have at least 2 columns. A sample identifier column, and the rest which correspond with data to visualise')

  # Figure out column classes
  data_columns <- colnames(data)[colnames(data) != col_sample]
  data_column_classes <- vapply(data_columns, function(colname){ class(data[[colname]]) }, FUN.VALUE = "")

  supported_data_classes <- c(
    numerical = 'integer',
    numerical = 'numeric',
    categorical='character',
    categorical='factor'
  )

  if(any(!data_column_classes %in% supported_data_classes)){
   cli::cli_abort(
     c(
       '!'="gg_sample_annotation does not yet support data of class: [{data_column_classes[!data_column_classes %in% supported_data_classes]}]",
       '>'="Supported classes include: [{supported_data_classes}]"
     ))
  }

  data_column_class_broad = names(supported_data_classes)[match(data_column_classes, supported_data_classes)]

  cli::cli_h1("Sample Metadata")
  cli::cli_alert_info('Found {length(data_columns)} data columns:')

  for (i in seq_along(data_columns)) {
    cli::cli_alert('{data_columns[i]} ({.cls {data_column_class_broad[i]}})')
    #cli::cli_alert('{data_columns[i]}: Specific Type = {data_column_classes[i]}')
  }

  # Now before we plot, lets sort samples based on `col_sort`
  # If col_sort is numeric, lets sort directly on the number
  col_sort_class_broad <- data_column_class_broad[match(col_sort, data_columns)]

  if(col_sort_class_broad == "categorical"){
    #data[[col_sort]] <- forcats::fct_infreq(data[[col_sort]])
    #data[[col_sample]] <- forcats::fct_reorder(data[[col_sample]], .x = data[[col_sort]])
    cli::cli_abort("Sorting based on categorical variable is not yet supported")
  }
  else if(col_sort_class_broad == "numerical"){
    cli::cli_alert_info("Sorting samples by {col_sort}")
    data[[col_sample]] <- forcats::fct_reorder(data[[col_sample]], .x = data[[col_sort]], .desc = TRUE)
  }
  else
    cli::cli_abort('col_sort_class_broad should always be categorical or numerical - please open a github issue')


  # Visualise each variable as a separate plot. Tile plots for categorical vars. Column graphs for numerics
  gglist <- list()
  for (i in seq_along(data_columns)) {
    current_col <- data_columns[i]
    current_col_data_class <- data_column_class_broad[i]

    # Categorical Data
    if(current_col_data_class == "categorical"){
      nlevels = dplyr::n_distinct(data[[current_col]])

      showlegend=TRUE

      if(nlevels > 8)
        showlegend=FALSE


      gglist[[current_col]] <- (
        ggplot2::ggplot(
          data,
          mapping = ggplot2::aes(
            x = .data[[col_sample]],
            y = {{current_col}},
            fill = .data[[current_col]],
            tooltip = paste0("Sample: ", .data[[col_sample]], "<br/>", current_col,": ", .data[[current_col]]),
            data_id = .data[[col_sample]]
            )
        ) +
          ggiraph::geom_tile_interactive() +
          sample_annotations_theme(legend.show = showlegend) +
          ggplot2::ylab("")

      )
    }
    # Numerical Data
      else if(current_col_data_class == "numerical"){
        gglist[[current_col]] <- (
          ggplot2::ggplot(
            data,
            mapping = ggplot2::aes(
              x = .data[[col_sample]],
              y = .data[[current_col]],
              tooltip = paste0("Sample: ", .data[[col_sample]], "<br/>", current_col,": ", .data[[current_col]]),
              data_id = .data[[col_sample]]
            )
          ) +
            ggiraph::geom_col_interactive() +
            sample_annotations_theme() +
            ggplot2::ylab("")

        )
      }
      else{
       cli::cli_abort('Failed to recognise current_col_data_class [{current_col_data_class}]. Expected numerical/categorical')
      }

    }


  gridplot = cowplot::plot_grid(plotlist = gglist, align = 'v', ncol = 1, axis = "lr")

  if(interactive){
    interactive_gridplot = ggiraph::girafe(
      ggobj = gridplot,
      width_svg = 8,
      height_svg = 4,
      ggiraph::opts_hover_inv(css = "opacity:0.1;")
      )
    return(interactive_gridplot)
    }
  return(gridplot)

}

sample_annotations_theme <- function(legend.show = TRUE, ...){
  legend_position = ifelse(legend.show, "right", "none")
  t <- ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = legend_position,
    ...
    )

  return(t)
}
# devtools::load_all();brca_metadata_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_clinicaldata.csv")
# brca_metadata_df <- read.csv(file = brca_metadata_csv, header=TRUE)
# gg_sample_annotation(brca_metadata_df, 'Tumor_Sample_Barcode')

# brca_rnaseq_path = system.file(package = "oncoplotgg","testdata/BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.data.txt")
# brca_rnaseq_df <- data.table::fread(brca_rnaseq, header = TRUE, sep = '\t')
