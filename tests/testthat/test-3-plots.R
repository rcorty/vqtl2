context('Testing plotting')

test_that(
  desc = 'plot_scan1var',
  code = {

    F2_plot <- plot_scan1var(s1v = tiny_F2_s1v,
                             genetic_map = tiny_F2_cross$gmap)

    expect_is(object = F2_plot, class = 'ggplot')


    DO_plot <- plot_scan1var(s1v = tiny_DO_s1v,
                             genetic_map = tiny_DO_cross$gmap)

    expect_is(object = DO_plot, class = 'ggplot')
  }
)



test_that(
  desc = 'plot_allele_effects',
  code = {

    F2_ae_plot <- plot_allele_effects(s1v = tiny_F2_s1v,
                                      marker = 'D19Mit68')

    expect_is(object = F2_ae_plot, class = 'ggplot')

    DO_ae_plot <- plot_allele_effects(s1v = tiny_DO_s1v,
                                      marker = 'UNC180249267')

    expect_is(object = DO_ae_plot, class = 'ggplot')
  }
)


test_that(
  desc = 'plot_allele_effects_over_region',
  code = {

    F2_ae_plot <- plot_allele_effects_scan(s1v = tiny_F2_s1v,
                                           genetic_map = tiny_F2_cross$gmap,
                                           chr = '19')
    # start_marker = 'D18Mit20',
    # stop_marker = 'D18Mit186')

    expect_is(object = F2_ae_plot, class = 'ggplot')

    DO_ae_plot <- plot_allele_effects_scan(s1v = tiny_DO_s1v,
                                           genetic_map = tiny_DO_cross$gmap,
                                           chr = '18')

    expect_is(object = DO_ae_plot, class = 'ggplot')
  }
)
