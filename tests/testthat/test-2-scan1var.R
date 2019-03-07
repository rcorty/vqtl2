context('Testing scan1var')

library(qtl2)
library(vqtl2)


testthat::test_that(
  desc = 'F2 experiment',
  code = {

    iron <- read_cross2(file = system.file("extdata", "iron.zip",
                                           package = "qtl2"))
    iron <- subset(x = iron, chr = c(17, 18, 19))
    iron_map <- insert_pseudomarkers(map = iron$gmap, step = 1)
    iron_gp <- calc_genoprob(cross = iron, map = iron_map, error_prob = 0.002)
    iron_ap <- genoprob_to_alleleprob(probs = iron_gp)

    s1v <- scan1var(pheno_name = 'liver',
                    mean_covar_names = 'spleen',
                    alleleprobs = iron_ap,
                    non_genetic_data = as.data.frame(iron$pheno))

    expect_true(object = is_scan1var(x = s1v))

    expect_equal(object = nrow(x = s1v),
                 expected = sum(sapply(X = iron_map, FUN = length)))

    s1v <- scan1var(pheno_name = 'liver',
                    mean_covar_names = 'spleen',
                    alleleprobs = iron_ap,
                    non_genetic_data = as.data.frame(iron$pheno),
                    num_cores = 0)

    expect_true(object = is_scan1var(x = s1v))

    expect_equal(object = nrow(x = s1v),
                 expected = sum(sapply(X = iron_map, FUN = length)))
  }
)



# # do population
# gatti_file <- 'https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Gatti2014/do.zip'
#
# gatti_cross <- read_cross2(file = gatti_file)
#
# small_do_cross <- subset(x = gatti_cross, ind = 1:100, chr = 1:5)
#
# map <- insert_pseudomarkers(small_do_cross$gmap, step = 10)
#
# pr <- calc_genoprob(cross = small_do_cross, map = map, quiet = FALSE)
#
# apr <- genoprob_to_alleleprob(probs = pr, quiet = FALSE)
#
# so <- scan1(genoprobs = apr, pheno = small_do_cross$pheno)
#
# plot(x = so, map = small_do_cross$pmap)
#
#
# sov <- scan1var(pheno_name = 'WBC',
#                 mean_covar_names = 'NEUT',
#                 var_covar_names = 'NEUT',
#                 alleleprobs = apr,
#                 non_genetic_data = as.data.frame(x = small_do_cross$pheno))
#
# plot(x = sov, map = small_do_cross$pmap)
#
#
# # riself
# grav2 <- read_cross2( system.file("extdata", "grav2.zip", package="qtl2") )
