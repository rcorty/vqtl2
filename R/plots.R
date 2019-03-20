#' @title Plot allele effects
#'
#' @param s1v the scan1var object with the results to be plotted
#' @param which_marker the marker at which the allele effects will be plotted
#'
#' @return the plot
#' @export
#' @importFrom dplyr %>%
#'
plot_allele_effects <- function(s1v,
                                which_marker) {

  marker <- key <- value <- 'fake global for CRAN'
  meanvar_estimse <- meanvar <- estimse <- 'fake global for CRAN'
  mean_estim <- var_estim <- allele <- mean_se <- var_se <- 'fake global for CRAN'

  purrr::cross(list(c('mean', 'var'),
                    attr(x = s1v, which = 'alleles'),
                    c('estim', 'se'))) %>%
    purrr::map_chr(paste, collapse = '_') ->
    allele_effect_names

  s1v %>%
    dplyr::filter(marker == which_marker) %>%
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
                                   which_marker))


}
