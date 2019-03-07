# utility for grabbing "..." args
#
# dotargs = list(...) from function call
# argname = character string of the argument to grab
# default = default value for argument
# values  = optional vector of character strings with possible values
grab_dotargs <-
	function(dotargs, argname, default, values=NULL)
	{
		if(argname %in% names(dotargs)) {
			arg <- dotargs[[argname]]
			if(!is.null(values) && !(arg %in% values)) {
				warning(argname, ' "', arg, '" not valid; using "',
						default, '".')
				arg <- default
			}
		}
		else arg <- default
		arg
	}

# is a number? (from assertthat)
is_number <- function(x) is.numeric(x) && length(x)==1
is_nonneg_number <- function(x) is_number(x) && x >= 0
is_pos_number <- function(x) is_number(x) && x > 0


# check that objects have rownames, if they are matrices
#   (or names, if not matrices)
check4names <-
	function(pheno=NULL, addcovar=NULL, Xcovar=NULL, intcovar=NULL, nullcovar=NULL)
	{
		args <- list(pheno=pheno,
					 addcovar=addcovar,
					 Xcovar=Xcovar,
					 intcovar=intcovar,
					 nullcovar=nullcovar)

		for(i in seq_along(args)) {
			a <- args[[i]]
			if(!is.null(a)) {
				if(is.matrix(a)) {
					if(is.null(rownames(a)))
						stop(names(args)[i], " has no rownames")
				}
				else {
					if(is.null(names(a)))
						stop(names(args)[i], " has no names")
				}
			} # end if(!is.null)
		} # end loop

		TRUE
	}

sqrt_weights <-
	function(weights, tol=1e-12)
	{
		if(is.null(weights)) return(weights)

		if(any(weights <= 0))
			stop("weights must all be positive")

		if(all(!is.na(weights) & abs(weights - 1)<tol))
			return(NULL)

		return(sqrt(weights))
	}


get_common_ids <-
	function(..., complete.cases=FALSE)
	{
		args <- list(...)
		if(length(args)==0) {
			return(character(0))
		}

		# find the IDs in common across all
		id <- NULL
		for(i in seq_along(args)) {
			if(is.null(args[[i]])) next

			if(is.matrix(args[[i]]) || is.data.frame(args[[i]]) || is.array(args[[i]])) {
				if(length(dim(args[[i]])) > 3)
					stop("Can't handle arrays with >3 dimensions")
				these <- rownames(args[[i]])
				if(complete.cases && (is.matrix(args[[i]]) || is.data.frame(args[[i]])))
					these <- these[rowSums(!is.finite(args[[i]]))==0]
			}
			else if(is.list(args[[i]]) && !is.null(rownames(args[[i]][[1]]))) {
				these <- rownames(args[[i]][[1]])
			}
			else if(is.vector(args[[i]])) {
				if(is.character(args[[i]])) {
					if(is.null(names(args[[i]]))) {
						these <- args[[i]]
					} else {
						these <- names(args[[i]])
						if(complete.cases) {
							these <- these[!is.na(args[[i]])]
						}
					}
				}
				else {
					these <- names(args[[i]])
					if(complete.cases) {
						these <- these[is.finite(args[[i]])]
					}
				}
			}
			else if(is.character(args[[i]])) { # character but not vector
				if(is.null(names(args[[i]]))) {
					these <- args[[i]]
				} else {
					these <- names(args[[i]])
					if(complete.cases) {
						these <- these[!is.na(args[[i]])]
					}
				}
			}
			else if(!is.null(names(args[[i]]))) { # not a vector but has names
				these <- names(args[[i]])
			}
			else {
				stop("Not sure what to do with object of class ", class(args[[i]]))
			}

			if(length(unique(these)) != length(these))
				stop("Duplicate names in argument ", i)

			if(is.null(id)) id <- these
			else id <- id[id %in% these]
		}

		id
	}


setup_cluster <-
	function(cores, quiet=TRUE)
	{
		if(is_cluster(cores)) return(cores)

		if(is.null(cores) || is.na(cores)) cores <- 1
		if(cores==0) cores <- parallel::detectCores() # if 0, detect cores
		if(is.na(cores)) cores <- 1

		if(cores > 1 && Sys.info()[1] == "Windows") { # windows doesn't support mclapply
			cores <- parallel::makeCluster(cores)
			# the following calls on.exit() in the function that called this one
			# see http://stackoverflow.com/a/20998531
			do.call("on.exit",
					list(quote(parallel::stopCluster(cores))),
					envir=parent.frame())
		}
		cores
	}

is_cluster <-
  function(cores)
  {
    "cluster" %in% class(cores) && "SOCKcluster" %in% class(cores)
  }
