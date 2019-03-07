#
# pull_additive_component <- function(gp) {
#
# 	tibble::tibble(add = switch(
# 		attr(x = gp, which = 'crosstype'),
# 		'f2' = gp[,2,] + 2*gp[,3,],
# 		'bc' = gp[,2],
# 		'do' =
# 	)
# 	)
# }


tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

fit_dglm <- function(mf, df, data, family, wts = NULL) {
  dglm::dglm(formula = mf,
             dformula = df,
             data = data,
             method = 'ml',
             family = family,
             ykeep = FALSE)
}

fit_hglm <- function(mf, df, data, glm_family) {#, obs_weights) {
  stop('commented out because hglm package archived on CRAN')
  # hglm::hglm2(meanmodel = mf, disp = df, data = data, calc.like = TRUE, family = glm_family)#, weights = obs_weights)
}

fit_dhglm <- function(mf, df, data) {
  stop('dhglm not yet implemented.')
}

fit_model <- function(formulae,
                      data,
                      mean = c('alt', 'null'),
                      var = c('alt', 'null'),
                      model = c('dglm', 'hglm', 'dhglm'),
                      glm_family = c('gaussian', 'poisson'),
                      permute_what = c('none', 'mean', 'var', 'both'),
                      the.perm = seq(from = 1, to = nrow(data)),
                      obs_weights = rep(1, nrow(data))) {

  mean <- match.arg(arg = mean)
  var <- match.arg(arg = var)
  model <- match.arg(arg = model)
  glm_family <- match.arg(arg = glm_family)
  permute_what <- match.arg(arg = permute_what)

  mf <- switch(mean, alt = formulae[['mean.alt.formula']], null = formulae[['mean.null.formula']])
  vf <- switch(var, alt = formulae[['var.alt.formula']], null = formulae[['var.null.formula']])

  fit_model <- switch(EXPR = model,
                      dglm = fit_dglm,
                      hglm = fit_hglm,
                      dhglm = fit_dhglm)

  glm_family <- switch(EXPR = glm_family,
                       gaussian = stats::gaussian,
                       poisson = stats::poisson)

  data <- switch(EXPR = permute_what,
                 none = data,
                 mean = permute.mean.QTL.terms_(df = data, the.perm = the.perm),
                 var = permute.var.QTL.terms_(df = data, the.perm = the.perm),
                 both = permute.QTL.terms_(df = data, the.perm = the.perm))

  tryNA(fit_model(mf = mf, df =  vf, data = data, glm_family = glm_family))#, obs_weights = obs_weights)
  # tryNA(do.call(what = fit_model,
  # args = list(mf = mf, df =  vf, data = data, glm_family = glm_family, weights = obs_weights)))
}



log_lik <- function(f) {

  if (inherits(x = f, what = 'dglm')) {
    if (abs(f$m2loglik) > 1e8) { return(NA) }
    return(-0.5*f$m2loglik)
  }
  if (inherits(x = f, what = 'hglm')) {
    if (abs(f$likelihood$hlik) > 1e8) { return(NA) }
    return(f$likelihood$hlik)
  }
  return(stats::logLik(object = f))
}

LRT <- function(alt, null) {

  if (any(identical(alt, NA), identical(null, NA))) {
    return(NA)
  }

  if (!identical(class(alt), class(null))) {
    stop('Can only calculate LOD on models of the same class.')
  }

  if (!inherits(x = alt, what = c('dglm', 'hglm'))) {
    stop('Can only calcualte LOD on models of class dglm or hglm.')
  }

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
  if (inherits(x = f, what = 'dglm')) {
    length(coef(f)) + length(coef(f$dispersion.fit)) - 2L
  }
}

pull_allele_names <- function(apr) {
  dimnames(x = apr)[[2]]
}

pull_marker_names <- function(apr) {
  dimnames(x = apr)[[3]]
}

make_formula <- function(response_name = NULL,
                         covar_names) {

  if (all(is.null(response_name), is.null(covar_names))) {
    covar_names <- 1
  }

  formula(
    paste(response_name,
          '~',
          paste(
            covar_names,
            collapse = '+'
          )
    )
  )

}
