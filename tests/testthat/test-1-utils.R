context('Testing utilities')

library(qtl2)
library(vqtl2)

testthat::test_that(
  desc = 'make_formula',
  code = {

    expect_equal(object = vqtl2:::make_formula(),
                 expected = stats::formula(~ 1))

    expect_equal(object = vqtl2:::make_formula(response_name = 'a'),
                 expected = stats::formula(a ~ 1))

    expect_equal(object = vqtl2:::make_formula(covar_names = c('b', 'c')),
                 expected = stats::formula(~ 0 + b + c))

    expect_equal(object = vqtl2:::make_formula(response_name = 'a',
                                               covar_names = c('b', 'c')),
                 expected = stats::formula(a ~ 0 + b + c))

  }
)
