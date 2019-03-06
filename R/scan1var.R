#' Mean-Variance genome scan with a single-QTL model
#'
#' @param alleleprobs Genotype probabilities as calculated by
#'     [qtl::calc_genoprob()].
#' @param pheno A data.frame of phenotypes, individuals x phenotypes.
#' @param mean_addcovar An optional data.frame of additive covariates.
#' @param var_addcovar An optional data.frame of additive covariates.
#' @param weights An optional numeric vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param model Indicates whether to use a normal model (least
#'     squares) or binary model (logistic regression) for the phenotype.
#'     If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' @param ... additional optional arguments
#'
#' @return results of the scan
#' @export
#'
scan1var <- function(pheno_name,
                     mean_covar_names = NULL,
                     var_covar_names = NULL,
                     alleleprobs,
                     non_genetic_data,
                     model = c("normal", "binary"),
                     weights = NULL,
                     cores = 1,
                     ...)
{

  model <- match.arg(arg = model)
  family <- switch(EXPR = model,
                   'normal' = stats::gaussian,
                   'binary' = stats::binomial)

  mean_null_formula <- make_formula(response_name = pheno_name, covar_names = mean_covar_names)
  var_null_formula <- make_formula(covar_names = var_covar_names)

  null_fit <- scan1var_nullfit(formulae = list(mean_null = mean_null_formula,
                                               var_null = var_null_formula),
                               data = non_genetic_data,
                               family = family)

  dplyr::bind_rows(
    lapply(X = alleleprobs,
           FUN = scan1var_onechr,
           pheno_name = pheno_name,
           mean_covar_names = mean_covar_names,
           var_covar_names = var_covar_names,
           non_genetic_data = non_genetic_data,
           family = family,
           null_fit = null_fit)
  )
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
      data = dplyr::bind_cols(as.data.frame(alleleprobs[,,mn]),
                              non_genetic_data),
      family = family,
      null_fit = null_fit)
  }

  return(dplyr::bind_rows(results))

}

scan1var_nullfit <- function(formulae,
                             data,
                             family) {

  fit_dglm(mf = formulae$mean_null,
           df = formulae$var_null,
           data = data,
           family = family)

}

scan1var_onelocus <- function(marker_name,
                              formulae,
                              data,
                              family,
                              null_fit) {

  mv = tryNA(
    fit_dglm(mf = formulae$mean_alt,
             df = formulae$var_alt,
             data = data,
             family = family)
  )

  m = tryNA(
    fit_dglm(mf = formulae$mean_alt,
             df = formulae$var_null,
             data = data,
             family = family)
  )

  v = tryNA(
    fit_dglm(mf = formulae$mean_null,
             df = formulae$var_alt,
             data = data,
             family = family)
  )

  tibble::tibble(marker = marker_name,
                 mv_lr = LRT(alt = mv, null = null_fit),
                 mv_dof = dof(f = mv) - dof(null_fit),
                 m_lr = LRT(alt = mv, null = v),
                 m_dof = dof(f = mv) - dof(v),
                 v_lr = LRT(alt = mv, null = m),
                 v_dof = dof(f = mv) - dof(m))
}

