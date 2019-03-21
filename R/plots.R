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

  # combine gmap from cross with s1v to give each marker a location
  tibble::tibble(
    chr = rep(x = names(cross$gmap), times = sapply(X = cross$gmap, length)),
    marker = unlist(x = sapply(X = cross$gmap, FUN = names), use.names = FALSE),
    loc = unlist(x = cross$gmap, use.names = FALSE)
  ) %>%
    dplyr::inner_join(y = s1v, by = c('chr', 'marker')) ->
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

  input_marker <- marker
  marker <- key <- value <- 'fake global for CRAN'
  meanvar_estimse <- meanvar <- estimse <- 'fake global for CRAN'
  mean_estim <- var_estim <- allele <- mean_se <- var_se <- 'fake global for CRAN'

  purrr::cross(list(c('mean', 'var'),
                    attr(x = s1v, which = 'alleles'),
                    c('estim', 'se'))) %>%
    purrr::map_chr(paste, collapse = '_') ->
    allele_effect_names

  s1v %>%
    dplyr::filter(marker == input_marker) %>%
    dplyr::select(dplyr::matches(paste(allele_effect_names, collapse = '|'))) %>%
    tidyr::gather(key = key, value = value) %>%
    tidyr::separate(col = 'key',
                    sep = '_',
                    into = c('meanvar', 'allele', 'estimse')) %>%
    tidyr::unite(col = meanvar_estimse, meanvar, estimse) %>%
    tidyr::spread(key = meanvar_estimse, value = value) ->
    plotting_data

  ggplot2::ggplot(data = plotting_data) +
    ggplot2::geom_hline(yintercept = 0, color = 'darkgray') +
    ggplot2::geom_vline(xintercept = 0, color = 'darkgray') +
    ggplot2::geom_point(mapping = ggplot2::aes(x = mean_estim, y = var_estim, color = allele),
                        size = 2) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = mean_estim - mean_se,
                                                 xend = mean_estim + mean_se,
                                                 y = var_estim,
                                                 yend = var_estim,
                                                 color = allele),
                          size = 2,
                          alpha = 0.5) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = mean_estim,
                                                 xend = mean_estim,
                                                 y = var_estim - var_se,
                                                 yend = var_estim + var_se,
                                                 color = allele),
                          size = 2,
                          alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::xlab(label = 'mean effects') +
    ggplot2::ylab(label = 'variance effects') +
    ggplot2::ggtitle(label = paste('Allele effects on',
                                   attr(x = s1v, which = 'pheno_name'),
                                   'at',
                                   marker))


}

