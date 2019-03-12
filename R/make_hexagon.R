# install.packages('hexSticker')
# install.packages('extrafont')
remove.packages('sysfonts'); install.packages('sysfonts')

library(sysfonts)
library(extrafont)
library(hexSticker)
library(ggplot2)

#
# font_import()
fonts()
# fonttable()
# fonttable()[40:45,]

# RColorBrewer::display.brewer.all()

pal <- RColorBrewer::brewer.pal(n = 3, name = 'Set1')
pal2 <- RColorBrewer::brewer.pal(n = 8, name = 'YlGnBu')

p <- ggplot(mapping = aes(x = x),
       data = tibble::tibble(x = c(-2, 2))) +
  stat_function(fun = dnorm, n = 101, geom = 'area', fill = pal[3], alpha = 0.8, args = list(sd = 1.3)) +
  stat_function(fun = dnorm, n = 101, geom = 'area', fill = pal[2], alpha = 0.6, args = list(sd = 0.8)) +
  stat_function(fun = dnorm, n = 101, geom = 'area', fill = pal[1], alpha = 0.3, args = list(sd = 0.5)) +
  # annotate(geom = 'text', x = 0, y = 0.55, label = 'vqtl2', size = 10, family = 'Ubuntu Mono', color = pal2[8]) +
  theme_void() +
  theme_transparent()

p

sticker(subplot = p,
        package = 'vqtl2',
        p_size = 8,
        s_x = 1,
        s_y = 1.15,
        s_width = 1.8,
        s_height = 1.4,
        p_x = 1,
        p_y = 1.3,
        p_color = 'darkblue',
        h_fill = 'white',
        h_color = pal[2],
        url = 'https://github.com/rcorty/vqtl2',
        filename = 'man/figures/hex_logo.png')
