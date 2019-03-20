context('Testing plotting')

test_that(
  desc = 'plot_allele_effects',
  code = {

    F2_ae_plot <- plot_allele_effects(s1v = tiny_F2_s1v,
                                      which_marker = 'D19Mit68')

    expect_is(object = F2_ae_plot, class = 'ggplot')

    DO_ae_plot <- plot_allele_effects(s1v = tiny_DO_s1v,
                                      which_marker = 'UNC180249267')

    expect_is(object = DO_ae_plot, class = 'ggplot')
  }
)
