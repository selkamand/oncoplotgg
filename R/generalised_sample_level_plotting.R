
#' gg_sample_annotation
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
#'
#' # Try experimenting with sorting - you can sort on multiple columns
#' brca_metadata_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_clinicaldata.csv")
#' brca_metadata_df <- read.csv(file = brca_metadata_csv, header=TRUE)
#' gg_sample_annotation(brca_metadata_df, 'Tumor_Sample_Barcode', col_sort = c("vital_status", "age_at_initial_pathologic_diagnosis"))
gg_sample_annotation <- function(data, col_sample, col_sort, interactive = TRUE, svg_width=10, svg_heigt=4, numerical_plot_height_ratio = 2){

  #assertions
  assertthat::assert_that(assertthat::is.string(col_sample))
  assertthat::assert_that(is.data.frame(data))
  assertthat::assert_that(assertthat::is.number(numerical_plot_height_ratio))
  assertthat::assert_that(is.character(col_sort))

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
    categorical='factor',
    categorical='logical'
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
  #browser()
  sample_order_ranks =  heirarchical_ranking(data[col_sort])
  data[[col_sample]] <- forcats::fct_reorder(data[[col_sample]], .x = sample_order_ranks, .desc = TRUE)

  # Visualise each variable as a separate plot. Tile plots for categorical vars. Column graphs for numerics
  gglist <- list()
  for (i in seq_along(data_columns)) {
    current_col <- data_columns[i]
    current_col_data_class <- data_column_class_broad[i]

    # Categorical Data
    if(current_col_data_class == "categorical"){
      nlevels = dplyr::n_distinct(data[[current_col]])
      scale_fill <- ggplot2::scale_fill_discrete()
      showlegend=TRUE

      # Hide legend if there are too many levels
      if(nlevels > 6)
        showlegend=FALSE

      if(is.logical(data[[current_col]])){
        showlegend=FALSE
        scale_fill <- ggplot2::scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightgrey"))
      }


      gglist[[current_col]] <- (
        ggplot2::ggplot(
          data,
          mapping = ggplot2::aes(
            x = .data[[col_sample]],
            y = {{current_col}},
            fill = .data[[current_col]],
            tooltip = paste0(
              {{col_sample}},": ", .data[[col_sample]], "<br/>",
              {{current_col}},": <b>", .data[[current_col]]),
            data_id = .data[[col_sample]]
          )
        ) +
          ggiraph::geom_tile_interactive() +
          sample_annotations_theme(legend.show = showlegend) +
          #ggplot2::xlab(NULL) +
          ggplot2::ylab(current_col) +
          ggplot2::guides(fill = ggplot2::guide_legend(keywidth = 0.2, keyheight = .5, label.theme = ggplot2::element_text(size = 8))) +
          scale_fill
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
            tooltip = paste0({{col_sample}},": ", .data[[col_sample]], "<br/>", {{current_col}},": <b>", .data[[current_col]]),
            data_id = .data[[col_sample]]
          )
        ) +
          ggiraph::geom_col_interactive() +
          sample_annotations_theme() +
          #ggplot2::xlab(NULL) +
          ggplot2::ylab(current_col)

      )
    }
    else{
      cli::cli_abort('Failed to recognise current_col_data_class [{current_col_data_class}]. Expected numerical/categorical')
    }

  }

  # Add xlab title to bottom plot
  #gglist[[length(gglist)]] <- gglist[[length(gglist)]] + ggplot2::xlab(col_sample)

  # Make numerical Variables twice the height of categoricals
  rel_heights = ifelse(data_column_class_broad == "numerical", yes = numerical_plot_height_ratio, no = 1)

  # Combine plots
  gridplot = cowplot::plot_grid(plotlist = gglist, align = 'v', ncol = 1, axis = "lr", rel_heights = rel_heights)

  if(interactive){
    interactive_gridplot = ggiraph::girafe(
      ggobj = gridplot,
      width_svg = svg_width,
      height_svg = svg_heigt,
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
    axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5),
    axis.text.y = ggplot2::element_blank(),
    ...
  )

  return(t)
}


#' Ranking based on multiple data columns
#'
#' @param data_ls
#'
#' @return
#' @export
#'
#' @examples
heirarchical_ranking <- function(data_ls){
  if (is.vector(data_ls))
    data_ls <- list(data_ls)

  if(is.data.frame(data_ls)){
    data_ls <- as.list(data_ls)
  }

  #browser()
  # Check list length
  ncols = length(data_ls)

  if(ncols == 0)
    cli::cli_abort("List must have at least one column. Found {ncols}")

  data_vec_lengths = lengths(data_ls)

  if(dplyr::n_distinct(data_vec_lengths) > 1)
    cli::cli_abort('Data vectors that you sort on must all be the same length. Observed lengths: [{data_vec_lengths}]')


  # Refactor any categorical variables in order of sequence
  data_ls <- lapply(data_ls, FUN = function(vec){
    if(!is.numeric(vec) & !is.logical(vec)){
      return(forcats::fct_rev(forcats::fct_infreq(vec)))
    }
    else
      return(vec)
    })

  # Create Ranking
  rankings_ls = lapply(
    seq_along(data_ls),
    FUN = function(i) {
      #browser()
      scales::rescale(rank(unlist(data_ls[[i]]), ties.method = "min"), to = c(0, 1))/(10^(i-1))
      })


  #browser()
  rankings_final = numeric(data_vec_lengths[1])

  for (i in seq_along(rankings_ls)) {
    rankings_final <- rankings_final + rankings_ls[[i]]
  }

  rank(rankings_final, ties.method = "first")
}

# devtools::load_all();brca_metadata_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_clinicaldata.csv")
# brca_metadata_df <- read.csv(file = brca_metadata_csv, header=TRUE)
# gg_sample_annotation(brca_metadata_df, 'Tumor_Sample_Barcode', col_sort = c("vital_status", "gender"))

# devtools::load_all();brca_metadata_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_clinicaldata.csv")
# brca_metadata_df <- read.csv(file = brca_metadata_csv, header=TRUE)
# gg_sample_annotation(brca_metadata_df, 'Tumor_Sample_Barcode', col_sort = c("vital_status", "age_at_initial_pathologic_diagnosis"))


# devtools::load_all();brca_metadata_csv <- system.file(package='oncoplotgg', "testdata/BRCA_tcgamutations_mc3_clinicaldata.csv")
# brca_metadata_df <- read.csv(file = brca_metadata_csv, header=TRUE)
# gg_sample_annotation(brca_metadata_df, 'Tumor_Sample_Barcode', col_sort = c("vital_status", "breast_carcinoma_estrogen_receptor_status"))



