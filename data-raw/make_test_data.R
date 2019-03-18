library(qtl2)

# todo: add X chromosomes in to sample data whenever functionality to handle them gets added

# get F2 cross directly from qtl2
F2_cross <- read_cross2(file = system.file('extdata', 'iron.zip', package = 'qtl2'))
print(x = F2_cross)
print(x = object.size(x = F2_cross), units = 'Mb')

tiny_F2_cross <- subset(x = F2_cross, chr = c(17, 18, 19))
print(x = object.size(x = F2_cross), units = 'Mb')
usethis::use_data(tiny_F2_cross)

tiny_F2_map <- insert_pseudomarkers(map = tiny_F2_cross$gmap, step = 10)
print(x = object.size(x = tiny_F2_map), units = 'Mb')

tiny_F2_gp <- calc_genoprob(cross = tiny_F2_cross, map = tiny_F2_map, error_prob = 0.002)
tiny_F2_ap <- genoprob_to_alleleprob(probs = tiny_F2_gp)
print(x = object.size(x = tiny_F2_ap), units = 'Mb')
usethis::use_data(tiny_F2_ap)



# get DO "cross" from qtl2data
# accessed from commit 74baa4509a2671bc707cac560fc7cbf4ceb6c0c0
DO_cross <- read_cross2(file = paste0('https://raw.githubusercontent.com/rqtl/',
                                      'qtl2data/master/DO_Gatti2014/do.zip'))
print(x = DO_cross)
print(x = object.size(x = DO_cross), units = 'Mb')

tiny_DO_cross <- subset(x = DO_cross, chr = c('18', '19'), ind = 1:100)
print(x = tiny_DO_cross)
print(x = object.size(x = tiny_DO_cross), units = 'Mb')
usethis::use_data(tiny_DO_cross)

tiny_DO_map <- insert_pseudomarkers(map = tiny_DO_cross$gmap, step = 10)
print(x = object.size(x = tiny_DO_map), units = 'Mb')

tiny_DO_gp <- calc_genoprob(cross = tiny_DO_cross, map = tiny_DO_map, error_prob = 0.002)
tiny_DO_ap <- genoprob_to_alleleprob(probs = tiny_DO_gp)
print(x = object.size(x = tiny_DO_ap), units = 'Mb')
usethis::use_data(tiny_DO_ap)