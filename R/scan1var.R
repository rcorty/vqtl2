#' Genome scan with a single-QTL model
#'
#' Genome scan with a single-QTL model by Haley-Knott regression or a
#' linear mixed model, with possible allowance for covariates.
#'
#' @md
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [qtl::calc_genoprob()].
#' @param pheno A numeric matrix of phenotypes, individuals x phenotypes.
#' @param addcovar An optional numeric matrix of additive covariates.
#' @param intcovar An numeric optional matrix of interactive covariates.
#' @param weights An optional numeric vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param model Indicates whether to use a normal model (least
#'     squares) or binary model (logistic regression) for the phenotype.
#'     If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#'
#' @export
#'
#' @return A matrix of LOD scores, positions x phenotypes.
#' Also contains one or more of the following attributes:
#' * `sample_size` - Vector of sample sizes used for each
#'    phenotype
scan1var <-
	function(genoprobs,
			 pheno,
			 addcovar = NULL,
			 Xcovar = NULL,
			 intcovar = NULL,
			 weights = NULL,
			 reml = TRUE,
			 model = c("normal", "binary"),
			 cores = 1,
			 ...)
	{
		scan1(genoprobs = genoprobs,
			  pheno = pheno,
			  addcovar = addcovar,
			  Xcovar = Xcovar,
			  intcovar = intcovar,
			  weights = weights,
			  model = model,
			  cores = cores,
			  ...)
	}