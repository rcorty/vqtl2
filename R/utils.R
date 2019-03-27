tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

tryNULL <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NULL,
                            finally = NULL))
}

fit_dglm <- function(mf, vf, locus_data, family, error_silently = TRUE) {

  # this didn't work -- some problem with how dglm eval's namespaces?
  # fit_dglm_ <- purrr::compose(ifelse(test = error_silently, yes = tryNA, no = identity),
  #                             dglm::dglm)
  #
  # fit_dglm_(formula = mf,
  #           dformula = vf,
  #           data = force(locus_data),
  #           method = 'ml',
  #           family = family,
  #           ykeep = FALSE)

  if (error_silently) {
    tryNULL(
      dglm::dglm(
        formula = mf,
        dformula = vf,
        data = locus_data,
        method = 'ml',
        family = family,
        ykeep = FALSE)
    )
  } else {
    dglm::dglm(
      formula = mf,
      dformula = vf,
      data = locus_data,
      method = 'ml',
      family = family,
      ykeep = FALSE
    )
  }
}

fit_hglm <- function(mf, df, data, glm_family) {#, obs_weights) {
  stop('commented out because hglm package archived on CRAN')
  # hglm::hglm2(meanmodel = mf, disp = df, data = data, calc.like = TRUE, family = glm_family)#, weights = obs_weights)
}

log_lik <- function(f) {

  if (inherits(x = f, what = 'dglm')) {
    if (abs(f$m2loglik) > 1e8) { return(NA) }
    return(-0.5*f$m2loglik)
  }
  if (inherits(x = f, what = 'hglm')) {
    stop('no hglm for now.')
    # if (abs(f$likelihood$hlik) > 1e8) { return(NA) }
    # return(f$likelihood$hlik)
  }
  return(stats::logLik(object = f))
}

LRT <- function(alt, null) {

  if (any(is.null(alt), is.null(null)))
    return(NA)

  if (!identical(class(alt), class(null)))
    stop('Can only calculate LOD on models of the same class.')

  if (!inherits(x = alt, what = c('dglm', 'hglm')))
    stop('Can only calcualte LOD on models of class dglm or hglm.')

  LRT <- 2*(log_lik(alt) - log_lik(null))

}

LOD <- function(alt, null) {
  return(0.5*LRT(alt = alt, null = null)/log(10))
}

LOD_from_LLs <- function(null_ll, alt_ll) {
  return((alt_ll - null_ll)/log(10))
}

LRT_from_LLs <- function(null_ll, alt_ll) {
  return(2*(alt_ll - null_ll))
}

dof <- function(f) {

  if (is.null(f)) { return(NA) }

  if (inherits(x = f, what = 'dglm')) {
    length(stats::coef(f)) + length(stats::coef(f$dispersion.fit)) - 2L
  }
}

pull_allele_names <- function(apr) {
  dimnames(x = apr)[[2]]
}

pull_marker_names <- function(apr) {
  dimnames(x = apr)[[3]]
}

make_formula <- function(response_name = NULL,
                         covar_names = '1') {

  stats::as.formula(
    paste(response_name,
          '~',
          paste(
            covar_names,
            collapse = '+'
          )
    )
  )

}

prepend_classes <- function(x, new_classes) {
  `class<-`(x, c(new_classes, class(x = x)))
}

add_attribute <- function(x, which, value) {
  `attr<-`(x, which = which, value = value)
}

conditionally <- function(fun){
  function(first_arg, ..., execute){
    if(execute) return(fun(first_arg, ...))
    else return(first_arg)
  }
}

cond_filter <- conditionally(dplyr::filter)
cond_select <- conditionally(dplyr::select)
cond_mutate <- conditionally(dplyr::mutate)
cond_slice <- conditionally(dplyr::slice)

pull_effects <- function(model,
                         which_submodel = c('mean', 'var', 'both')) {

  term <- estimate <- std.error <- 'fake global for CRAN'
  measure <- val <- united <- 'fake global for CRAN'

  if (is.null(model)) { return(tibble::tibble()) }

  which_submodel <- match.arg(arg = which_submodel)
  if (which_submodel == 'var') { model <- model$dispersion.fit }
  if (which_submodel == 'both') {
    return(
      dplyr::bind_cols(pull_effects(model = model, which_submodel = 'mean'),
                       pull_effects(model = model, which_submodel = 'var'))
    )
  }

  model %>%
    broom::tidy() %>%
    dplyr::mutate(term = dplyr::case_when(term == '(Intercept)' ~ 'intercept',
                                          TRUE ~ term)) %>%
    cond_mutate(term = paste0(which_submodel, '_', term),
                execute = !is.null(which_submodel)) %>%
    dplyr::select(term, estimate, std.error) %>%
    tidyr::gather(key = measure, value = val, estimate, std.error) %>%
    dplyr::mutate(measure = dplyr::case_when(measure == 'estimate' ~ 'estim',
                                             measure == 'std.error' ~ 'se')) %>%
    tidyr::unite(col = 'united', term, measure) %>%
    tidyr::spread(key = united, value = val)
}

compute_allele_effects <- function(s1v, markers) {

  `.` <- marker <- key <- value <- 'fake global for CRAN'
  meanvar_estimse <- meanvar <- estimse <- 'fake global for CRAN'
  mean_estim <- var_estim <- 'fake global for CRAN'

  input_markers <- markers

  cross_w_meanvar_and_estimse <- function(v) {
    purrr::cross(list(c('mean', 'var'),
                      v,
                      c('estim', 'se'))) %>%
      purrr::map_chr(paste, collapse = '_')
  }

  alleles <- attr(x = s1v, which = 'alleles')
  first_allele_names <- cross_w_meanvar_and_estimse(alleles[1])
  other_allele_names <- cross_w_meanvar_and_estimse(alleles[-1])
  regex_for_oafn <- paste(other_allele_names, collapse = '|')

  s1v %>%
    dplyr::filter(marker %in% input_markers) %>%
    dplyr::select(dplyr::matches('marker'),
                  dplyr::matches('loc'),
                  dplyr::matches(regex_for_oafn)) %>%
    tidyr::gather(key = key, value = value, dplyr::matches(regex_for_oafn)) %>%
    dplyr::bind_rows(
      tibble::tibble(marker = rep(x = markers, each = 4),
                     key = rep(x = first_allele_names, times = length(markers)),
                     value = 0),
      .
    ) %>%
    tidyr::separate(col = 'key',
                    sep = '_',
                    into = c('meanvar', 'allele', 'estimse')) %>%
    tidyr::unite(col = meanvar_estimse, meanvar, estimse) %>%
    tidyr::spread(key = meanvar_estimse, value = value) %>%
    dplyr::group_by(marker) %>%
    dplyr::mutate(mean_estim_cent = center(x = mean_estim),
                  var_estim_cent = center(x = var_estim))
}

center <- function(x) {
  x - mean(x)
}

cc_colors <- c(
  grDevices::rgb(red = 240, green = 240, blue = 000, maxColorValue = 255),
  grDevices::rgb(red = 128, green = 128, blue = 128, maxColorValue = 255),
  grDevices::rgb(red = 240, green = 128, blue = 128, maxColorValue = 255),
  grDevices::rgb(red = 016, green = 016, blue = 240, maxColorValue = 255),
  grDevices::rgb(red = 000, green = 160, blue = 240, maxColorValue = 255),
  grDevices::rgb(red = 000, green = 160, blue = 000, maxColorValue = 255),
  grDevices::rgb(red = 240, green = 000, blue = 000, maxColorValue = 255),
  grDevices::rgb(red = 144, green = 000, blue = 224, maxColorValue = 255)
)

join_s1v_gmap <- function(s1v, gmap) {

  tibble::tibble(
    chr = rep(x = names(gmap), times = sapply(X = gmap, length)),
    marker = unlist(x = sapply(X = gmap, FUN = names), use.names = FALSE),
    loc = unlist(x = gmap, use.names = FALSE)
  ) %>%
    dplyr::inner_join(y = s1v, by = c('chr', 'marker'))

}



theme_vqtl2 <- function() {
  theme_minimal() +
    theme(panel.background = element_rect(fill = '#DDDDDD', color = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = '#BBBBBB'),
          strip.placement = 'outside',
          strip.background = element_rect(fill = '#DDDDDD',
                                          color = NA),
          plot.title = element_text(size = 16,
                                    family = 'Ubuntu'))
}
