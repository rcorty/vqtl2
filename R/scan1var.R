#' Mean-Variance genome scan with a single-QTL model
#'
#' @param pheno_name name of the phenotype to scan
#' @param mean_covar_names names of the mean covariates
#' @param var_covar_names names of the variance covariates
#' @param alleleprobs Genotype probabilities as calculated by
#'     [qtl::calc_genoprob()].
#' @param non_genetic_data phenotype and covararite data.frame
#' @param model Indicates whether to use a normal model (least
#'     squares) or binary model (logistic regression) for the phenotype.
#'     If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' @param weights An optional numeric vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param num_cores Number of CPU cores to use, for parallel calculations.
#'     (If `0`, use [parallel::detectCores() - 1].)
#'
#' @param ... additional optional arguments
#'
#' @return results of the scan
#' @export
#'
#' @importFrom dplyr %>%
#'
scan1var <- function(pheno_name,
                     mean_covar_names = NULL,
                     var_covar_names = NULL,
                     alleleprobs,
                     non_genetic_data,
                     model = c("normal", "binary"),
                     weights = NULL,
                     num_cores = 1,
                     ...)
{

  model <- match.arg(arg = model)
  family <- switch(EXPR = model,
                   'normal' = stats::gaussian,
                   'binary' = stats::binomial)

  null_fit <- fit_dglm(mf = make_formula(response_name = pheno_name, covar_names = mean_covar_names),
                       vf = make_formula(covar_names = var_covar_names),
                       locus_data = non_genetic_data,
                       family = family,
                       error_silently = FALSE)

  num_cores <- ifelse(test = num_cores == 0,
                      yes = parallel::detectCores() - 1,
                      no = num_cores)

  if (num_cores == 1) {
    result <- dplyr::bind_rows(
      lapply(X = alleleprobs,
             FUN = scan1var_onechr,
             pheno_name = pheno_name,
             mean_covar_names = mean_covar_names,
             var_covar_names = var_covar_names,
             non_genetic_data = non_genetic_data,
             family = family,
             null_fit = null_fit)
    )
  } else {

    # would have preferred to use mclapply but it didn't work
    cl <- parallel::makeCluster(spec = num_cores)
    doParallel::registerDoParallel(cl = cl)
    `%dopar%` <- foreach::`%dopar%`   # necessary to get the loop to parse
    result <-  foreach::foreach(i = 1:length(alleleprobs), .combine = dplyr::bind_rows) %dopar% {

      scan1var_onechr(pheno_name = pheno_name,
                      mean_covar_names = mean_covar_names,
                      var_covar_names = var_covar_names,
                      alleleprobs = alleleprobs[[i]],
                      non_genetic_data = non_genetic_data,
                      family = family,
                      null_fit = null_fit)
    }
    parallel::stopCluster(cl = cl)
  }

  null_mean_effects <- pull_effects(model = null_fit, effect_name_prefix = 'mean')
  null_var_effects <- pull_effects(model = null_fit$dispersion.fit, effect_name_prefix = 'var')

  dplyr::bind_cols(null_mean_effects, null_var_effects) %>%
    dplyr::bind_rows(result) %>%
    dplyr::select(marker,
           dplyr::ends_with(match = '_lr'),
           dplyr::ends_with(match = '_dof'),
           dplyr::everything()) %>%
    prepend_class(new_class = 'scan1') %>%
    prepend_class(new_class = 'scan1var')
}


scan1var_onechr <- function(pheno_name,
                            mean_covar_names,
                            var_covar_names,
                            alleleprobs,
                            pheno,
                            non_genetic_data,
                            family,
                            null_fit) {

  allele_names <- pull_allele_names(apr = alleleprobs)

  mean_alt_formula <- make_formula(response_name = pheno_name, covar_names = c(allele_names[-1], mean_covar_names))
  var_alt_formula <- make_formula(covar_names = c(allele_names[-1], var_covar_names))
  mean_null_formula <- make_formula(response_name = pheno_name, covar_names = mean_covar_names)
  var_null_formula <- make_formula(covar_names = var_covar_names)

  formulae <- list(mean_alt = mean_alt_formula,
                   var_alt = var_alt_formula,
                   mean_null = mean_null_formula,
                   var_null = var_null_formula)

  marker_names <- pull_marker_names(apr = alleleprobs)

  results <- list()

  for (mn in marker_names) {

    results[[mn]] <- scan1var_onelocus(
      marker_name = mn,
      formulae = formulae,
      locus_data = dplyr::bind_cols(as.data.frame(alleleprobs[,,mn]),
                                    non_genetic_data),
      family = family,
      null_fit = null_fit)
  }

  return(dplyr::bind_rows(results))

}


scan1var_onelocus <- function(marker_name,
                              formulae,
                              locus_data,
                              family,
                              null_fit) {

  mv = fit_dglm(mf = formulae$mean_alt,
                vf = formulae$var_alt,
                locus_data = locus_data,
                family = family)

  m = tryNA(
    fit_dglm(mf = formulae$mean_alt,
             vf = formulae$var_null,
             locus_data = locus_data,
             family = family)
  )

  v = tryNA(
    fit_dglm(mf = formulae$mean_null,
             vf = formulae$var_alt,
             locus_data = locus_data,
             family = family)
  )

  tibble::tibble(marker = marker_name,
                 mvqtl_lr = LRT(alt = mv, null = null_fit),
                 mvqtl_dof = dof(f = mv) - dof(null_fit),
                 mqtl_lr = LRT(alt = mv, null = v),
                 mqtl_dof = dof(f = mv) - dof(v),
                 vqtl_lr = LRT(alt = mv, null = m),
                 vqtl_dof = dof(f = mv) - dof(m)) %>%
    dplyr::bind_cols(pull_effects(model = mv, effect_name_prefix = 'mean')) %>%
    dplyr::bind_cols(pull_effects(model = mv$dispersion.fit, effect_name_prefix = 'var'))
}


#' Test whether an R object is a scan1var
#'
#' @param x the object to test
#'
#' @return TRUE if [x] is a scan1var, FALSE otherwise
#' @export
#'
is_scan1var <- function(x) {

  if (!identical(class(x), c('scan1var', 'scan1', 'tbl_df', 'tbl', 'data.frame')))
    return(FALSE)

  if (!identical(x = sapply(X = x, FUN = class)[1:7],
                 y = c(marker = 'character',
                       mvqtl_lr = 'numeric',
                       mqtl_lr = 'numeric',
                       vqtl_lr = 'numeric',
                       mvqtl_dof = 'integer',
                       mqtl_dof = 'integer',
                       vqtl_dof = 'integer')))
    return(FALSE)

  if (any(with(data = x, expr = stats::na.omit(c(mvqtl_lr, mqtl_lr, vqtl_lr))) < 0))
    return(FALSE)

  return(TRUE)
}
