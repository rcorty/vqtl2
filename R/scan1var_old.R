#' #' Genome scan with a single-QTL model
#' #'
#' #' Genome scan with a single-QTL model by Haley-Knott regression or a
#' #' linear mixed model, with possible allowance for covariates.
#' #'
#' #' @md
#' #'
#' #' @param alleleprobs Genotype probabilities as calculated by
#' #' [qtl::calc_genoprob()].
#' #' @param pheno A data.frame of phenotypes, individuals x phenotypes.
#' #' @param mean_addcovar An optional data.frame of additive covariates.
#' #' @param var_addcovar An optional data.frame of additive covariates.
#' #' @param weights An optional numeric vector of positive weights for the
#' #' individuals. As with the other inputs, it must have `names`
#' #' for individual identifiers.
#' #' @param model Indicates whether to use a normal model (least
#' #'     squares) or binary model (logistic regression) for the phenotype.
#' #'     If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' #' @param cores Number of CPU cores to use, for parallel calculations.
#' #' (If `0`, use [parallel::detectCores()].)
#' #'
#' #' @export
#' #'
#' #' @return A matrix of LOD scores, positions x phenotypes.
#' #' Also contains one or more of the following attributes:
#' #' * `sample_size` - Vector of sample sizes used for each
#' #'    phenotype
#' scan1var <- function(alleleprobs,
#' 					 pheno,
#' 					 mean_addcovar = NULL,
#' 					 var_addcovar = NULL,
#' 					 weights = NULL,
#' 					 reml = TRUE,
#' 					 model = c("normal", "binary"),
#' 					 cores = 1,
#' 					 ...)
#' {
#'
#' 	# check arguments
#' 	if(is.null(alleleprobs)) stop("alleleprobs is NULL")
#' 	if(is.null(pheno)) stop("pheno is NULL")
#' 	if(ncol(pheno) > 5) warning('More than 5 phenotypes is going to take a while...')
#'
#' 	model <- match.arg(model)
#'
#' 	dotargs <- list(...)
#' 	if ('n_perm' %in% names(dotargs))
#' 		stop("You included n_perm as an argument; you probably want to run scan1varperm not scan1var.")
#' 	if ('kinship' %in% names(dotargs))
#' 		stop("You included kinship as an argument; scan1var can't handle that.")
#'
#' 	# read and check dotargs
#' 	tol <- grab_dotarg(dotargs = dotargs, argname = 'tol', default = 1e-12)
#' 	if(!is_pos_number(tol)) stop("tol should be a single positive number")
#'
#' 	maxit <- grab_dotarg(dotargs = dotargs, argname = 'maxit', default = 100)
#'
#' 	quiet <- grab_dotarg(dotargs = dotargs, argname = 'quiet', default = TRUE)
#'
#' 	# check that the objects have rownames
#' 	check4names(pheno, addcovar, Xcovar, intcovar)
#'
#' 	# # force things to be matrices
#' 	# if(!is.matrix(pheno)) {
#' 	# 	pheno <- as.matrix(pheno)
#' 	# 	if(!is.numeric(pheno)) stop("pheno is not numeric")
#' 	# }
#' 	# if(is.null(colnames(pheno))) # force column names
#' 	# 	colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))
#' 	# if(!is.null(addcovar)) {
#' 	# 	if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
#' 	# 	if(!is.numeric(addcovar)) stop("addcovar is not numeric")
#' 	# }
#' 	# if(!is.null(Xcovar)) {
#' 	# 	if(!is.matrix(Xcovar)) Xcovar <- as.matrix(Xcovar)
#' 	# 	if(!is.numeric(Xcovar)) stop("Xcovar is not numeric")
#' 	# }
#' 	# if(!is.null(intcovar)) {
#' 	# 	if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)
#' 	# 	if(!is.numeric(intcovar)) stop("intcovar is not numeric")
#' 	# }
#'
#'
#' 	# find individuals in common across all arguments
#' 	# and drop individuals with missing covariates or missing *all* phenotypes
#' 	# ind2keep <- get_common_ids(genoprobs, addcovar, Xcovar, intcovar,
#' 	# 						   weights, complete.cases=TRUE)
#' 	# ind2keep <- get_common_ids(ind2keep, rownames(pheno)[rowSums(is.finite(pheno)) > 0])
#' 	# if(length(ind2keep)<=2) {
#' 	# 	if(length(ind2keep)==0)
#' 	# 		stop("No individuals in common.")
#' 	# 	else
#' 	# 		stop("Only ", length(ind2keep), " individuals in common: ",
#' 	# 			 paste(ind2keep, collapse=":"))
#' 	# }
#' 	#
#' 	# # make sure addcovar is full rank when we add an intercept
#' 	# addcovar <- drop_depcols(addcovar, TRUE, tol)
#' 	#
#' 	# # make sure columns in intcovar are also in addcovar
#' 	# addcovar <- force_intcovar(addcovar, intcovar, tol)
#' 	#
#' 	# # drop things from Xcovar that are already in addcovar
#' 	# Xcovar <- drop_xcovar(addcovar, Xcovar, tol)
#' 	#
#' 	# # batch phenotypes by missing values
#' 	# phe_batches <- batch_cols(pheno[ind2keep,,drop=FALSE], max_batch)
#' 	#
#' 	# # drop cols in genotype probs that are all 0 (just looking at the X chromosome)
#' 	# genoprob_Xcol2drop <- genoprobs_col2drop(genoprobs)
#' 	# is_x_chr <- attr(genoprobs, "is_x_chr")
#' 	# if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))
#'
#' 	# set up parallel analysis
#' 	cores <- setup_cluster(cores)
#' 	if(!quiet && n_cores(cores)>1) {
#' 		message(" - Using ", n_cores(cores), " cores")
#' 		quiet <- TRUE # make the rest quiet
#' 	}
#'
#'
#' 	#
#' 	# # batches for analysis, to allow parallel analysis
#' 	# run_batches <- data.frame(chr = rep(x = seq_len(length(genoprobs)), length(phe_batches)),
#' 	# 						  phe_batch = rep(x = seq_along(phe_batches), each = length(genoprobs)))
#' 	#
#' 	# run_indexes <- seq_len(length(genoprobs)*length(phe_batches))
#'
#' 	# the function that does the work
#' 	# by_group_func <- function(i) {
#' 	#
#' 	# 	# deal with batch information, including individuals to drop due to missing phenotypes
#' 	# 	chr <- run_batches$chr[i]
#' 	# 	chrnam <- names(genoprobs)[chr]
#' 	# 	phebatch <- phe_batches[[run_batches$phe_batch[i]]]
#' 	# 	phecol <- phebatch$cols
#' 	# 	omit <- phebatch$omit
#' 	# 	these2keep <- ind2keep # individuals 2 keep for this batch
#' 	# 	if(length(omit) > 0) these2keep <- ind2keep[-omit]
#' 	# 	if(length(these2keep)<=2) return(NULL) # not enough individuals
#' 	#
#' 	# 	# subset the genotype probabilities: drop cols with all 0s, plus the first column
#' 	# 	Xcol2drop <- genoprob_Xcol2drop[[chrnam]]
#' 	# 	if(length(Xcol2drop) > 0) {
#' 	# 		pr <- genoprobs[[chr]][these2keep,-Xcol2drop,,drop=FALSE]
#' 	# 		pr <- pr[,-1,,drop=FALSE]
#' 	# 	}
#' 	# 	else
#' 	# 		pr <- genoprobs[[chr]][these2keep,-1,,drop=FALSE]
#' 	#
#' 	# 	# subset the rest
#' 	# 	ac <- addcovar; if(!is.null(ac)) { ac <- ac[these2keep,,drop=FALSE]; ac <- drop_depcols(ac, TRUE, tol) }
#' 	# 	Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
#' 	# 	ic <- intcovar; if(!is.null(ic)) { ic <- ic[these2keep,,drop=FALSE]; ic <- drop_depcols(ic, TRUE, tol) }
#' 	# 	ph <- pheno[these2keep,phecol,drop=FALSE]
#' 	# 	wts <- weights[these2keep]
#' 	#
#' 	# 	# if X chr, paste X covariates onto additive covariates
#' 	# 	# (only for the null)
#' 	# 	if(is_x_chr[chr]) ac0 <- drop_depcols(cbind(ac, Xc), add_intercept=FALSE, tol)
#' 	# 	else ac0 <- ac
#' 	#
#' 	# 	if(add_intercept) {
#' 	# 		addcovar <- cbind(rep(1,n), addcovar) # add intercept
#' 	# 	}
#' 	#
#' 	#
#' 	#
#' 	# 	scan1var_clean(genoprobs, pheno, addcovar, intcovar, weights, tol)
#' 	#
#' 	# }
#'
#' 	# # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
#' 	# npos_by_chr <- dim(genoprobs)[3,]
#' 	# totpos <- sum(npos_by_chr)
#' 	# pos_index <- split(seq_len(totpos), rep(seq_len(length(genoprobs)), npos_by_chr))
#' 	#
#' 	# # object to contain the LOD scores; also attr to contain sample size
#' 	# result <- matrix(nrow=totpos, ncol=ncol(pheno))
#' 	# n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)
#' 	# if(totpos==0) { # edge case of no genoprobs
#' 	# 	colnames(result) <- colnames(pheno)
#' 	# 	attr(result, "sample_size") <- n
#' 	# 	class(result) <- c("scan1", "matrix")
#' 	# 	return(result)
#' 	# }
#'
#'
#' 	# null fit
#' 	# null_fit <- scan1var_onelocus()
#'
#' 	if(n_cores(cores) == 1) { # no parallel processing
#'
#' 		sapply(X = alleleprobs,
#' 			   FUN = scan1var_onechr)
#' 		#
#' 		# for(i in run_indexes) {
#' 		# 	chr <- run_batches$chr[i]
#' 		# 	chrnam <- names(genoprobs)[chr]
#' 		# 	phebatch <- phe_batches[[run_batches$phe_batch[i]]]
#' 		# 	phecol <- phebatch$cols
#' 		#
#' 		# 	this_result <- by_group_func(i)
#' 		# 	if(!is.null(this_result)) {
#' 		# 		result[pos_index[[chr]], phecol] <- t(this_result$lod)
#' 		# 		if(chr==1) n[phecol] <- this_result$n
#' 		# 	}
#' 		# }
#' 	}
#' 	# else {
#' 	# 	# calculations in parallel
#' 	# 	list_result <- cluster_lapply(cores, run_indexes, by_group_func)
#' 	#
#' 	# 	# check for problems (if clusters run out of memory, they'll return NULL)
#' 	# 	result_is_null <- vapply(list_result, is.null, TRUE)
#' 	# 	if(any(result_is_null))
#' 	# 		stop("cluster problem: returned ", sum(result_is_null), " NULLs.")
#' 	#
#' 	# 	# reorganize results
#' 	# 	for(i in run_indexes) {
#' 	# 		chr <- run_batches$chr[i]
#' 	# 		chrnam <- names(genoprobs)[chr]
#' 	# 		phebatch <- phe_batches[[run_batches$phe_batch[i]]]
#' 	# 		phecol <- phebatch$cols
#' 	#
#' 	# 		if(!is.null(list_result[[i]])) {
#' 	# 			result[pos_index[[chr]], phecol] <- t(list_result[[i]]$lod)
#' 	# 			if(chr==1) n[phecol] <- list_result[[i]]$n
#' 	# 		}
#' 	# 	}
#' 	# }
#'
#' 	# pos_names <- unlist(dimnames(genoprobs)[[3]])
#' 	# names(pos_names) <- NULL # this is just annoying
#' 	# dimnames(result) <- list(pos_names, colnames(pheno))
#' 	#
#' 	# # add some attributes with details on analysis
#' 	# attr(result, "sample_size") <- n
#' 	#
#' 	# class(result) <- c("scan1", "matrix")
#' 	# result
#' }
#'
#'
#' # turns 3d array into list of 2d arrays
#'
#' # scan1var function taking nicely aligned data with no missing values
#' # Here genoprobs is a list of data.frame's
#' #
#' scan1var_onechr <- function(alleleprobs,
#' 							pheno,
#' 							mean_addcovar,
#' 							var_addcovar,
#' 							family,
#' 							weights)
#' {
#'
#' 	for(p in names(pheno)) {
#'
#' 		mean_formula <- paste(
#' 			paste(p,
#' 				  '~',
#' 				  paste(
#' 				  	pull_allele_names(alleleprobs = alleleprobs),
#' 				  	mean_addcovar,
#' 				  	sep = '+'
#' 				  )
#' 			)
#' 		)
#'
#' 		browser()
#'
#' 	}
#'
#' 	fit_locus <- dglm(formula = mean_formula,
#' 					  dformula = disp_formula,
#' 					  family = family,
#' 					  data = this_locus_data,
#' 					  weights = weights)
#'
#' 	lapply(X = locus_list, FUN = fit_locus)
#'
#'
#'
#' }
#'
#'
#'
#'
#' #
#' #
#' # # work in progress
#' # DGLM_norm <- function(y,
#' # 					  X_mean,
#' # 					  X_disp,
#' # 					  w = rep(1, length(y)),
#' # 					  maxiter=20,
#' # 					  conv=1e-6)
#' # {
#' # 	convergence <- 1
#' # 	iter <- 0
#' # 	while (convergence > conv & iter < maxiter) {
#' # 		iter <- iter + 1
#' # 		w.old <- w
#' # 		mean_submodel <- lm.fit(x = w*X_mean,
#' # 								y = w*y)
#' # 		q <- hatvalues(mean_submodel)
#' # 		y2 <- resid(mean_submodel)^2/(1-q)
#' # 		glm2 <- glm.fit(x = X_disp,
#' # 						y = y2,
#' # 						family=Gamma(link=log),
#' # 						weights=(1-q)/2)
#' # 		w <- 1/fitted(glm2)
#' # 		convergence <- (max(abs(w.old-w)) + (summary(glm1)$sigma-1) )
#' # 	}
#' # 	return(list(mean=glm1, disp=glm2, iter=iter))
#' # }