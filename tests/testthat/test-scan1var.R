context('Testing scan1var')

gatti_file <- 'https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Gatti2014/do.zip'

gatti_cross <- read_cross2(file = gatti_file)

small_do_cross <- subset(x = gatti_cross, ind = 1:100, chr = 1:5)

map <- insert_pseudomarkers(small_do_cross$gmap, step = 10)

pr <- calc_genoprob(cross = small_do_cross, map = map, quiet = FALSE)

apr <- genoprob_to_alleleprob(probs = pr, quiet = FALSE)

so <- scan1(genoprobs = apr, pheno = small_do_cross$pheno)

plot(x = so, map = small_do_cross$pmap)


sov <- scan1var(genoprobs = apr, pheno = small_do_cross$pheno)

plot(x = sov, map = small_do_cross$pmap)
