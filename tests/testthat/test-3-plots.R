context('Testing scan1var')

library(qtl2)
library(vqtl2)

testthat::test_that(
  desc = 'DO experiment',
  code = {

    testthat::skip_on_cran()

    gatti_file <- 'https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Gatti2014/do.zip'

    gatti_cross <- read_cross2(file = gatti_file)

    small_do_cross <- subset(x = gatti_cross, ind = 1:100, chr = 19)

    map <- insert_pseudomarkers(small_do_cross$gmap, step = 10)

    pr <- calc_genoprob(cross = small_do_cross, map = map, quiet = FALSE)

    apr <- genoprob_to_alleleprob(probs = pr, quiet = FALSE)

    s1v <- scan1var(pheno_name = 'WBC',
                    mean_covar_names = 'NEUT',
                    var_covar_names = 'NEUT',
                    alleleprobs = apr,
                    non_genetic_data = tibble::as_tibble(x = small_do_cross$pheno))

    s1v %>% glimpse()

    effects_plot(s1v = s1v)


  }
)
