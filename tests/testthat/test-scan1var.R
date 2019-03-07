context('Testing scan1var')

library(qtl2)
library(vqtl2)

# f2 cross
# iron <- read_cross2(file = "https://kbroman.org/qtl2/assets/sampledata/iron/iron.zip")
iron <- read_cross2(file = system.file("extdata", "iron.zip", package = "qtl2"))
map <- insert_pseudomarkers(map = iron$gmap, step = 1)
pr <- calc_genoprob(cross = iron, map = map, error_prob = 0.002)
apr <- genoprob_to_alleleprob(probs = pr)
# out <- scan1(genoprobs = apr, pheno = iron$pheno)
#
# plot(out, map, lodcolumn = 1, col="slateblue")
# plot(out, map, lodcolumn = 2, col="violetred", add=TRUE)
# legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out), bg="gray90")

system.time(
  s1v_iron <- scan1var(pheno_name = 'liver',
                       mean_covar_names = 'spleen',
                       alleleprobs = apr,
                       non_genetic_data = as.data.frame(iron$pheno),
                       num_cores = 1)
)

system.time(
  s1v_iron <- scan1var(pheno_name = 'liver',
                       mean_covar_names = 'spleen',
                       alleleprobs = apr,
                       non_genetic_data = as.data.frame(iron$pheno),
                       num_cores = 4)
)

# do population
gatti_file <- 'https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Gatti2014/do.zip'

gatti_cross <- read_cross2(file = gatti_file)

small_do_cross <- subset(x = gatti_cross, ind = 1:100, chr = 1:5)

map <- insert_pseudomarkers(small_do_cross$gmap, step = 10)

pr <- calc_genoprob(cross = small_do_cross, map = map, quiet = FALSE)

apr <- genoprob_to_alleleprob(probs = pr, quiet = FALSE)

so <- scan1(genoprobs = apr, pheno = small_do_cross$pheno)

plot(x = so, map = small_do_cross$pmap)


sov <- scan1var(pheno_name = 'WBC',
                mean_covar_names = 'NEUT',
                var_covar_names = 'NEUT',
                alleleprobs = apr,
                non_genetic_data = as.data.frame(x = small_do_cross$pheno))

plot(x = sov, map = small_do_cross$pmap)


# riself
grav2 <- read_cross2( system.file("extdata", "grav2.zip", package="qtl2") )
