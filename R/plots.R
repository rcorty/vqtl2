plot.scan1var <- function(x,
                          genetic_map,
                          what = c('association', 'effects'),
                          association = c('LR', 'LOD', 'asymptotic_p', 'empirical_p'),
                          effects = c('alleles', 'covariates', 'both')) {
  stopifnot(is_scan1var(x))
  stopifnot(is_genetic_map(geneti_map))

  what <- match.arg(arg = what)

  # should unite gmap and x here, since both the next
  # calls will need to do it

  switch(EXPR = what,
         'association' = plot_scan1var_assoc(
           yada,
           yada,
           what = match.arg(association)
         ),
         'effects' = plot_scan1var_effects(
           yada,
           yada,
           what = match.arg(effects)
         )
  )
}


#' @title Plot scan1var
#'
#' @param s1v the scan1var object to be plotted
#' @param genetic_map the map giving he location of each marker
#'
#' @return the plot
#'
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#'
plot_scan1var_assoc <- function(s1v,
                                gmap)
{

  loc <- mvqtl_lr <- mqtl_lr <- vqtl_lr <- 'fake global for CRAN'
  chr <- qtl_type <- lr <- 'fake global for CRAN'

  join_s1v_gmap(s1v = s1v,
                gmap = genetic_map) %>%
    select(chr, loc, matches('_lr')) %>%
    gather(key = qtl_type, value = lr, matches('_lr')) %>%
    mutate(qtl_type = factor(x = qtl_type,
                             levels = c('mvqtl_lr', 'mqtl_lr', 'vqtl_lr'),
                             labels = c('mvQTL', 'mQTL', 'vQTL'))) %>%
    ggplot(mapping = aes(x = loc)) +
    geom_line(mapping = aes(y = lr, color = qtl_type), size = 1) +
    facet_grid(cols = vars(chr),
               scales = 'free_x',
               space = 'free_x',
               switch = 'x',
               labeller = label_both) +
    scale_color_manual(values = c('black', 'blue', 'red'),
                       guide = guide_legend(title = 'QTL type')) +
    labs(title = paste('scan1var:', attr(x = s1v, which = 'pheno_name')),
         y = 'likelihood ratio') +
    theme_vqtl2() +
    theme(axis.title.x = element_blank())

}


#' @title Plot allele effects
#'
#' @param s1v the scan1var object with the results to be plotted
#' @param marker the marker at which the allele effects will be plotted
#'
#' @return the plot
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#'
plot_allele_effects <- function(s1v,
                                marker)
{

  mean_estim_cent <- var_estim_cent <- 'fake global for CRAN'
  allele <- mean_se <- var_se <- 'fake global for CRAN'

  stopifnot(is_scan1var(s1v))
  stopifnot(length(marker) == 1)

  compute_allele_effects(s1v = s1v, markers = marker) %>%
    ggplot(mapping = aes(x = mean_estim_cent,
                         xend = mean_estim_cent,
                         y = var_estim_cent,
                         yend = var_estim_cent,
                         color = allele)) +
    geom_hline(yintercept = 0, color = 'darkgray') +
    geom_vline(xintercept = 0, color = 'darkgray') +
    geom_segment(mapping = aes(x = mean_estim_cent - mean_se,
                               xend = mean_estim_cent + mean_se),
                 size = 2,
                 alpha = 0.8,
                 lineend = 'round') +
    geom_segment(mapping = aes(y = var_estim_cent - var_se,
                               yend = var_estim_cent + var_se),
                 size = 2,
                 alpha = 0.8,
                 lineend = 'round') +
    geom_point(size = 4, color = 'black') +
    geom_point(size = 3) +
    {
      if (identical(attr(x = s1v, which = 'alleles'), LETTERS[1:8]))
        scale_color_manual(values = cc_colors)
    } +
    labs(title = paste('Allele effects on',
                       attr(x = s1v, which = 'pheno_name'),
                       'at',
                       marker),
         x = 'mean effects',
         y = 'variance effects') +
    theme_vqtl2()
}


#' @title plot allele effects over a region
#'
#' @param s1v temp
#' @param genetic_map temp
#' @param chr temp
#'
#' @return the plot

#' @export
#' @import ggplot2
#' @importFrom dplyr %>% mutate select pull case_when
#' @importFrom tidyr gather
#'
plot_scan1var_effects <- function(s1v,
                                  genetic_map,
                                  # start_marker = NULL,
                                  # stop_marker = NULL,
                                  chr = NULL)
{

  `.` <- marker <- loc <- allele <- 'fake global for CRAN'
  mean_estim_cent <- var_estim_cent <- 'fake global for CRAN'
  meanvar <- estim <- 'fake global for CRAN'

  input_chr <- as.character(x = chr)

  # would be better to filter first, but join_s1v_gmap()
  # can't handle that at present...todo

  s1v %>%
    pull(marker) %>%
    `[`(-1) %>%
    compute_allele_effects(s1v = s1v, markers = .) %>%
    mutate(chr = input_chr) %>%
    join_s1v_gmap(gmap = genetic_map) %>%
    select(loc, allele, mean_estim_cent, var_estim_cent) %>%
    gather(key = meanvar, value = estim, mean_estim_cent, var_estim_cent) %>%
    mutate(meanvar = case_when(
      meanvar == 'mean_estim_cent' ~ 'mean',
      meanvar == 'var_estim_cent' ~ 'variance')) %>%
    ggplot(mapping = aes(x = loc, color = allele)) +
    geom_hline(yintercept = 0, color = 'darkgray') +
    geom_line(mapping = aes(y = estim), size = 1) +
    facet_grid(rows = vars(meanvar), scales = 'free_y') +
    {
      if (identical(attr(x = s1v, which = 'alleles'), LETTERS[1:8]))
        scale_color_manual(values = cc_colors)
    } +
    labs(title = paste('Chromosome', input_chr,
                       'allele effects on',
                       attr(x = s1v, which = 'pheno_name')),
         x = 'location',
         y = NULL,
         color = 'allele') +
    theme_vqtl2()
}
