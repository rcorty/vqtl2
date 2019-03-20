context('Testing scan1var')

testthat::test_that(
  desc = 'F2 experiment',
  code = {

    s1v <- scan1var(pheno_name = 'liver',
                    mean_covar_names = 'spleen',
                    var_covar_names = 'spleen',
                    alleleprobs = tiny_F2_ap,
                    non_genetic_data = tibble::as_tibble(tiny_F2_cross$pheno))

    expect_true(object = is_scan1var(x = s1v))

    # result should have one row per locus, plus one for null fit
    expect_equal(object = nrow(x = s1v),
                 expected = sum(sapply(X = tiny_F2_ap, FUN = dim)[3,]) + 1)

    # some NA is expected in results, but should be < 10%
    expect_lt(object = mean(is.na(s1v$mvqtl_lr)), expected = 0.1)

    s1v <- scan1var(pheno_name = 'liver',
                    mean_covar_names = 'spleen',
                    alleleprobs = tiny_F2_ap,
                    non_genetic_data = tibble::as_tibble(tiny_F2_cross$pheno),
                    num_cores = 2)

    expect_true(object = is_scan1var(x = s1v))

    # result should have one row per locus, plus one for null fit
    expect_equal(object = nrow(x = s1v),
                 expected = sum(sapply(X = tiny_F2_ap, FUN = dim)[3,]) + 1)

    # some NA is expected in results, but should be < 10%
    expect_lt(object = mean(is.na(s1v$mvqtl_lr)), expected = 0.1)
  }
)


testthat::test_that(
  desc = 'DO experiment',
  code = {

    testthat::skip_on_cran()

    s1v <- scan1var(pheno_name = 'WBC',
                    mean_covar_names = 'NEUT',
                    var_covar_names = 'NEUT',
                    alleleprobs = tiny_DO_ap,
                    non_genetic_data = tibble::as_tibble(x = tiny_DO_cross$pheno))

    expect_true(object = is_scan1var(x = s1v))

    # result should have one row per locus, plus one for null fit
    expect_equal(object = nrow(x = s1v),
                 expected = sum(sapply(X = tiny_DO_ap, FUN = dim)[3,]) + 1)

    # some NA is expected in results, but should be < 10%
    expect_lt(object = mean(is.na(s1v$mvqtl_lr)), expected = 0.1)
  }
)

#
# # riself
# grav2 <- read_cross2( system.file("extdata", "grav2.zip", package="qtl2") )
