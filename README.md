# Microbiome Community Clustering

methods to cluster microbiome data, both for counts and for relative abundance data. The methods are mainly **Dirichlet Multinomial Clustering**, **WGCNA**, a more classic hierarchical clustering that relies on beta diversity estimation, for a start.


# Download Repository 

## UNIX & UNIX-like Command line

  1. Make soure you have git installed in your system 
  2. Navigate to where you want to save this repository (eg: `cd ~/path/to/binaries/`)
  3. `git clone https://github.com/g-antonello/microbiome_clustering`
 


# Dirichlet Multinomial (on COUNTS data)

This is based on the following [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030126), in my implementation, I have a three-step approach that aims to make a command line interface to run. CST = "Community State Type", basically the microbial cluster

  1. Makes the actual fit for each pre-assumed number of communities
  2. Extracts CSTs from the "best fit" (of your choosing)
  3. Makes some figures based on the number of CSTs chosen (still WIP)
  
 An example run would be:
 **STEP 1** 
 ```
Rscript ~/path/to/binaries/microbiome_clustering/DirichletMultinomial/DMM_CST_step1-fitting.R \
      -i input_phyloseq \
      -o output_directory \
      -p 0.01 \ # prevalence threshold \
      -d 10 \ # minimum detection threshold \
      -c 10 \ # number of communities to model (from 1 to c)
      -r 1000 \ # minimum reads per sample, all samples below this threshold will be removed
      -t ASV \ # the taxonomic level you want to model at. NB: this must be the column name of the `tax_table(input_phyloseq)`
      -v TRUE \ # verbosity
 ```
NB: except the I/O parameters, all the values put here are arbitrary, default values. After this analysis has finished, you can run the second step 

**STEP 2**
```
Rscript  ~/path/to/binaries/microbiome_clustering/DirichletMultinomial/DMM_CST_step2-BestFitCSTExtraction.R \
      -i output_directory/DMM_output_data/01_DMM_CST_Robj.Rda \
      -n 5 \ 
```

This `-n` parameter can be derived or is arbitrarily chosen by the script, see `Rscript DirichletMultinomial/DMM_CST_step1-fitting.R -h` for the help page, that tells you how it is chosen. you can choose to run it with a value of your choosing, based on the observation of the figure stored in `output_directory/DMM_output_figures/02_laplace_curve.pdf`

If you are interested in an actual visualization of the composition of each CST, you can run step 3, but it's not crucial:

**STEP 3**
```
Rscript  ~/path/to/binaries/microbiome_clustering/DirichletMultinomial/DMM_CST_step2-BestFitCSTExtraction.R \
      -i output_directory/DMM_output_data/ 
```

After `step 2`, you can use the `output_directory/DMM_output_data/06_CSTs_ready_for_metadata.csv` file to add to your metadata. Since you will need to install my  `gautils2` package, you will probably do something like:

```
# devtools::install_github("g-antonello/gautils2") # uncomment in case you need to install it
library(gautils2)

phyloseq_no_CST <- readRDS("phyloseq_object.RDS")

cst_info <- read.csv("output_directory/DMM_output_data/06_CSTs_ready_for_metadata.csv", header = TRUE)

phyloseq_with_CST <- gautils2::phy_add_metadata_variables(physeq = phyloseq_no_CST,
                                                          df = cst_info,
                                                          by = "IDs") # make sure you merge using the correct IDS!                      
```

# Latent Dirichlet Allocation

One of the new kids in the block by Susan and her group, see this [link](https://academic.oup.com/biostatistics/article/20/4/599/5032578) for the statistical underpinning.
