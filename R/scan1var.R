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
                     mean_covar_names = '1',
                     var_covar_names = '1',
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

  i <- marker <- 'fake global for CRAN'

  # doesn't work -- I don't know why.
  result <- parallel::mclapply(
    X = alleleprobs,
    FUN = scan1var_onechr,
    pheno_name = pheno_name,
    mean_covar_names = mean_covar_names,
    var_covar_names = var_covar_names,
    non_genetic_data = non_genetic_data,
    family = family,
    null_fit = null_fit,
    mc.cores = num_cores
  )

  # make result -- first row is null_fit, rest is locus fits
  dplyr::bind_cols(
    tibble::tibble(marker = 'null_fit',
                   mvqtl_lr = NA,
                   mqtl_lr = NA,
                   vqtl_lr = NA,
                   mvqtl_dof = NA,
                   mqtl_dof = NA,
                   vqtl_dof = NA),
    pull_effects(model = null_fit, which_submodel = 'mean'),
    pull_effects(model = null_fit, which_submodel = 'var')
  ) %>%
    dplyr::bind_rows(result) %>%
    prepend_class(new_class = 'scan1') %>%
    prepend_class(new_class = 'scan1var')
}


scan1var_onechr <- function(alleleprobs,
                            pheno_name,
                            mean_covar_names,
                            var_covar_names,
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
      locus_data = dplyr::bind_cols(
        as.data.frame(alleleprobs[,,mn]),
        non_genetic_data
      ),
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

  mv <- fit_dglm(mf = formulae$mean_alt,
                 vf = formulae$var_alt,
                 locus_data = locus_data,
                 family = family)

  m <- fit_dglm(mf = formulae$mean_alt,
                vf = formulae$var_null,
                locus_data = locus_data,
                family = family)

  v <- fit_dglm(mf = formulae$mean_null,
                vf = formulae$var_alt,
                locus_data = locus_data,
                family = family)

  tibble::tibble(marker = marker_name,
                 mvqtl_lr = LRT(alt = mv, null = null_fit),
                 mvqtl_dof = dof(f = mv) - dof(null_fit),
                 mqtl_lr = LRT(alt = mv, null = v),
                 mqtl_dof = dof(f = mv) - dof(v),
                 vqtl_lr = LRT(alt = mv, null = m),
                 vqtl_dof = dof(f = mv) - dof(m)) %>%
    dplyr::bind_cols(pull_effects(model = mv, which_submodel = 'mean')) %>%
    dplyr::bind_cols(pull_effects(model = mv, which_submodel = 'var'))
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
