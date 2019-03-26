#' @title Plot scan1var
#'
#' @param s1v the scan1var object to be plotted
#' @param cross the cross object
#'
#' @return the plot
#' @export
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes
#'
plot_scan1var <- function(s1v,
                          cross)
{

  loc <- mvqtl_lr <- mqtl_lr <- vqtl_lr <- 'fake global for CRAN'

  s1v %>%
    insert_gmap_locs(gmap = cross$gmap) ->
    plotting_data

  # todo: check that size of plotting_data is correct
  # ie only dropped one row from s1v (the null fit row)
  # and didn't drop any markers from cross$gmap

  plotting_data %>%
    ggplot2::ggplot(mapping = aes(x = loc)) +
    ggplot2::geom_line(mapping = aes(y = mvqtl_lr), color = 'black') +
    ggplot2::geom_line(mapping = aes(y = mqtl_lr), color = 'blue') +
    ggplot2::geom_line(mapping = aes(y = vqtl_lr), color = 'red') +
    ggplot2::facet_grid(rows = ~ chr, space = 'free_x', switch = 'x') +
    ggplot2::ylab('likelihood ratio') +
    ggplot2::ggtitle(label = paste('scan1var:', attr(x = s1v, which = 'pheno_name'))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.placement = 'outside',
                   strip.background = ggplot2::element_rect(fill = 'lightgray', color = NA),
                   axis.title.x = ggplot2::element_blank())

}


#' @title Plot allele effects
#'
#' @param s1v the scan1var object with the results to be plotted
#' @param marker the marker at which the allele effects will be plotted
#'
#' @return the plot
#' @export
#' @importFrom dplyr %>%
#'
plot_allele_effects <- function(s1v,
                                marker)
{

  mean_estim_cent <- var_estim_cent <- 'fake global for CRAN'
  allele <- mean_se <- var_se <- 'fake global for CRAN'

  pull_allele_effects(s1v = s1v, markers = marker) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = mean_estim_cent,
                                           xend = mean_estim_cent,
                                           y = var_estim_cent,
                                           yend = var_estim_cent,
                                           color = allele)) +
    ggplot2::geom_hline(yintercept = 0, color = 'darkgray') +
    ggplot2::geom_vline(xintercept = 0, color = 'darkgray') +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = mean_estim_cent - mean_se,
                                                 xend = mean_estim_cent + mean_se),
                          size = 2,
                          alpha = 0.5) +
    ggplot2::geom_segment(mapping = ggplot2::aes(y = var_estim_cent - var_se,
                                                 yend = var_estim_cent + var_se),
                          size = 2,
                          alpha = 0.5) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    {
      if (identical(attr(x = s1v, which = 'alleles'), LETTERS[1:8]))
        ggplot2::scale_color_manual(values = cc_colors)
    } +
    ggplot2::xlab(label = 'mean effects') +
    ggplot2::ylab(label = 'variance effects') +
    ggplot2::ggtitle(label = paste('Allele effects on',
                                   attr(x = s1v, which = 'pheno_name'),
                                   'at',
                                   marker))
}


#' @title plot allele effects over a region
#'
#' @param s1v temp
#' @param cross temp
#' @param chr temp
#'
#' @return the plot
#' @export
#'
#' @importFrom dplyr %>%
#'
plot_allele_effects_over_region <- function(s1v,
                                            cross,
                                            # start_marker = NULL,
                                            # stop_marker = NULL,
                                            chr = NULL)
{

  `.` <- marker <- loc <- allele <- 'fake global for CRAN'
  mean_estim_cent <- var_estim_cent <- 'fake global for CRAN'
  meanvar <- estim <- 'fake global for CRAN'

  input_chr <- chr

  # s1v %>%
  #   cond_filter(chr == input_chr,
  #               execute = !is.null(input_chr)) %>%
  #   cond_slice(seq.int(from = which(marker == start_marker),
  #                      to = which(marker == stop_marker)),
  #              execute = is.null(input_chr))

  s1v %>%
    dplyr::pull(marker) %>%
    `[`(-1) %>%
    pull_allele_effects(s1v = s1v, markers = .) %>%
    dplyr::mutate(chr = input_chr) %>%
    insert_gmap_locs(gmap = cross$gmap) %>%
    dplyr::select(loc, allele, mean_estim_cent, var_estim_cent) %>%
    tidyr::gather(key = meanvar, value = estim, mean_estim_cent, var_estim_cent) %>%
    dplyr::mutate(meanvar = dplyr::case_when(
      meanvar == 'mean_estim_cent' ~ 'mean',
      meanvar == 'var_estim_cent' ~ 'var')) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = loc, color = allele)) +
    ggplot2::geom_hline(yintercept = 0, color = 'darkgray') +
    ggplot2::geom_line(mapping = ggplot2::aes(y = estim, linetype = meanvar)) +
    # ggplot2::geom_line(mapping = ggplot2::aes(y = var_estim_cent)) +
    ggplot2::theme_minimal() +
    {
      if (identical(attr(x = s1v, which = 'alleles'), LETTERS[1:8]))
        ggplot2::scale_color_manual(values = cc_colors)
    } +
    ggplot2::labs(title = 'temporary title',
                  x = 'location',
                  y = 'effect',
                  color = 'allele',
                  linetype = 'effect type')
}
