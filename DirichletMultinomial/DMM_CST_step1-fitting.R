t0 <- Sys.time()

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=FALSE, help = "The microbiome input data. So far, only phyloseq format is supported"),
  make_option(c("-o", "--output"), type="character", default=FALSE, help = "the output directory where you want results to be stored"),
  make_option(c("-p", "--prevalence"), type="double", default=0.01, help = "Feature's minimum prevalence across samples"),
  make_option(c("-d", "--detection"), type="integer", default=10, help = "minimum counts for each feature, based on the input prevalence"),
  make_option(c("-c", "--communities"), type="integer", default=10, help = "The total amount of communities you want to try to model"),
  make_option(c("-r", "--min_reads_per_sample"), type="integer", default=1000, help = "The total amount of communities you want to try to model"),
  make_option(c("-t", "--taxonomic_level"), type="character", default="ASV", help = "The total amount of communities you want to try to model"),
  make_option(c("-v", "--verbose"), type="logical", default="TRUE", help = "The total amount of communities you want to try to model")
)

opt_parser = OptionParser(option_list=option_list)

opts_args = parse_args(opt_parser)


print(opts_args)

library(gautils2)
chrismb_phy <- readRDS(file = opts_args$input)

# more than 1000 reads per sample
chrismb_phy_CST <- subset_samples(chrismb_phy, sample_sums(chrismb_phy) > opts_args$min_reads_per_sample) %>% 
  core(detection = opts_args$detection, prevalence = opts_args$prevalence) %>% 
  phy_tax_glom(opts_args$taxonomic_level)



# only ASV more abundant than 0.1% - this is totally arbitrary. bear in mind that this is useful if you want to run the algorithm on high resultion taxonomies(species or otu)
# abundance_crit <- sum(taxa_sums(chrismb_phy_CST))*0.001
# chrismb_phy_CST <- prune_taxa(taxa_sums(chrismb_phy_CST)> abundance_crit, chrismb_phy_CST)


dir.create(file.path(opts_args$output, "DMM_output_data"), showWarnings = FALSE, recursive = T)
dir.create(file.path(opts_args$output, "DMM_output_figures"), showWarnings = FALSE, recursive = T)

# Dirichilet part 1: prepare ASV data matrix

#The output of the following chunk should be fed to the script for Dirichlet multinomial on the high computational power cluster

### necessary for example dataset - for certain functions below, tells R to load pre-made data instead.

#prepare ASV table with correct formatting
chrismb_asv_4_CST <- chrismb_phy_CST %>% 
  abundances() %>% 
  t()

cnts <- log10(colSums(chrismb_asv_4_CST))
# deal with -Inf counts of 0 

#plot the taxon representation and save
library(lattice)

if(opts_args$verbose){
  message("Plot 1...")
}

{pdf(file.path(opts_args$output, "DMM_output_figures", "01_taxon-counts.pdf"))
  densityplot(cnts, xlim=range(cnts), xlab="Taxon representation (log 10 count)")
  dev.off()}

rm(cnts)

###### Identifying number of community types in study population ######
### fit count data to a series of models with 1 to 10 community types

#############################################
#NB: Do this step on a high-computing power machine, which can run ideally one model per thread, at least 4G memory each

library(DirichletMultinomial)
library(parallel)

if(opts_args$verbose){
  message("\nFitting the DMM models\n")
}

fit <- mclapply(X = 1:10, FUN = dmn, count=chrismb_asv_4_CST,mc.cores = data.table::getDTthreads()*2, verbose=TRUE)
  
### The return value can be queried for measures of fit (Laplace, AIC, BIC); these are plotted for different k
### plot laplace approximation (i.e. model fit) by number of community types
# put laplace values in a vector - you can use other measures of fit such as AIC or BIC

lplc <- sapply(fit, laplace)

{pdf(file.path(opts_args$output, "DMM_output_figures", "02_laplace_curve.pdf"))
  plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
  dev.off()}

save.image(file=file.path(opts_args$output, "DMM_output_data", "DMM_CST_Robj.Rda"))

if(opts_args$verbose){
message("Modeling finished")
  message(sprintf(Sys.time()-t0))
  }

