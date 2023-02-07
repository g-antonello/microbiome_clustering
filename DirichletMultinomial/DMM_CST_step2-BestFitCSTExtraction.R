
t0 <- Sys.time()

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=FALSE, help = "The result of step 1 Dirichlet Multinomial Modeling"),
  make_option(c("-n", "--number_components_bestFit"), type="integer", default=NA, help = "The number of communities that minimize the laplace curve, see the output figures after step 1 script. By default it will choose the CST after which the slope decrement does not half anymore")
)


opt_parser = OptionParser(option_list=option_list)

opts_args_st2 = parse_args(opt_parser)

print(opts_args_st2)

load(opts_args_st2$input)


if(is.na(opts_args_st2$number_components_bestFit)){
  # if the number of clusters is not specified, use the BIC Goodness of Fit method, which is arguably better at resolving the number of clusters in a DMM fit
 opts_args_st2$number_components_bestFit <- match(min(sapply(fit, DirichletMultinomial::BIC)), sapply(fit, DirichletMultinomial::BIC))
  
}

library(lattice)
library(DirichletMultinomial)
library(gautils2)


ls()
###### Characterizing best model ######
###Once you have inspected your Laplace plot, pull out best model and save as "best"
best <- fit[[opts_args_st2$number_components_bestFit]]

######
# important: make sure to EXAMINE the plot since you need to consider whether the difference in model fit is worth
#            complicating your analysis with additional community types. For example, let's say your plot looks
#            like an L with the "turn" at the model with 3 community types. Even though the best model fit may be
#            for the model with 10 community types you probably want to stick with 3 instead. This is a similar 
#            concept as forward selection model building. Do you really want to add parameters that don't improve 
#            your model by much? Code for project data will include this consideration.
######

### check weight (pi) and homogeneity (theta) of fitted model
mixturewt(best)
# weight (pi) = proportion of population with this community type
# homogenity (theta) = larger values are more homogenous. Not sure what the values determine whether something is 
#              is homogenous vs. hetergenous. Original Holmes paper describes 30.2 and 18.7 as "highly variable"
#              and 52.0 and 53.3 as "homogenous". Seems rather arbitrary.


### check probability of samples being each community type
head(mixture(best), 5)

# save probability table
write.csv(mixture(best), file.path(opts_args$output, "DMM_output_data", "03_bestcomponent_probab_table.csv"))

library(lattice)
### plot contribution of each taxonomic group to community types
{pdf(file.path(opts_args$output, "DMM_output_data", "contrib_taxon_for_each_CST.pdf"))
  splom(log(fitted(best)))
  dev.off()}

# diagonal nature of points suggest commnity types are correlated

### poster mean difference - how different are your community types from the population average (i.e. one community type)
p0 <- fitted(fit[[1]], scale=TRUE) #scale by theta
pN <- fitted(best, scale=TRUE) #scale by theta
colnames(pN) <- paste("m", 1:opts_args_st2$number_components_bestFit, sep="") # change "6" to whatever your maximum  number of community types is.

meandiff <- colSums(abs(pN - as.vector(p0)))
sum(meandiff)
# Range: [0 - 2] where 0 is completely identical and 2 is completely dissimilar to the population average.

write.csv(meandiff, file=file.path(opts_args$output, "DMM_output_data", "02_meandiff.csv"))

### taxonomic contributions to community types
diff <- rowSums(abs(pN - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
differ_k1_kN_coeff <- cumsum(diff[o]) / sum(diff)
df <- (cbind(Mean=p0[o], pN[o,], diff=diff[o], differ_k1_kN_coeff))
colnames(df) <- c("k1_coeff", 
                  paste0("k",opts_args_st2$number_components_bestFit, "_com", 1:opts_args_st2$number_components_bestFit), 
                  paste0("diff_k1_k", opts_args_st2$number_components_bestFit), 
                  "cumulat_difference_taxon")
head(df, 5)

# mean = relative abundance of taxa in population average
# m1-m4 = relative abundance of taxa in each community type
# diff = contribution of taxa to the total mean difference (i.e. meandiff) 
# differ_k1_k6_coeff = cumulative contribution
# original paper was able to get 95% credible intervals which were "calculated as the maximum posterior estimate minus/plus 
# two standard deviations as calculated from the inverse Hessian" - not sure how to do that with this package.

# save table as csv
write.csv(df %>% rownames_to_column("microbial_feature"), file=file.path(opts_args$output, "DMM_output_data", paste0("04_taxon_contribution_", opts_args_st2$number_components_bestFit, "_model.csv")))


### plot heatmap of read counts for samples and community types

{pdf(file.path(opts_args$output, "DMM_output_figures", "03_heatmap_reads_Samples.pdf"))
  heatmapdmn(chrismb_asv_4_CST, fit[[1]], best, 30)
  dev.off()
  }

# Explanation: narrow columns are samples. Broader columns are the community types. Color represents square-root counts,
#              with darker colors representing larger counts.

# Dirichilet part 2 

###### Assigning Community Types & Consider Misclassification
### assign a community type to each sample based on maximum probability
components <- data.frame(
  CST = apply(mixture(best), 1, function(x) which(x==max(x))),
  CST_assign_prob = apply(mixture(best), 1, function(x) max(x, na.rm = TRUE))
  )



### compare the best community type with the next best community
# transform to matrix
probmatrix<-as.matrix(mixture(best)[,1:opts_args_st2$number_components_bestFit])            # change "6" to the total number of community types in your dataset

# sort the matrix so the last column has decreasing probability
probmatrix_sorted<-as.data.frame(t(apply(probmatrix, 1, sort)))
head(probmatrix_sorted)

# subtract probability of best community type by probability of next best community type
probmatrix_sorted$probdiff<-probmatrix_sorted[[paste0("V",opts_args_st2$number_components_bestFit)]] - probmatrix_sorted[[paste0("V",opts_args_st2$number_components_bestFit-1)]]    # change "6" to the last column and change "5" to second to last 


### misclassification criteria (arbitrary: >=80% probability of classification. >=10% probability than any other community type)
# check how many are <90%
table(components$CST_assign_prob<.9)
prop.table(table(components$CST_assign_prob<.9))

# check how many are <80%
table(components$CST_assign_prob<.8)
prop.table(table(components$CST_assign_prob<.8))

# check how many are >10% difference
table(probmatrix_sorted$probdiff<.1)
prop.table(table(probmatrix_sorted$probdiff<.1))

# make misclassified = unknown when:
# 1. probability of assignment is lower than 0.8
# 2. difference in prob. of 2 most probable communities is less than 0.1
excluded_IDs <- rownames(probmatrix_sorted)[probmatrix_sorted[[paste0("V",opts_args_st2$number_components_bestFit)]] < .8 | probmatrix_sorted$probdiff<.1]

components$CST[row.names(components) %in% excluded_IDs] <-"Other"
components$CST_assign_prob[row.names(components) %in% excluded_IDs]<-NA
head(components, 10)

###############
# Important: there is no established exclusion criteria. For this template, we are assuming any sample which has a probability of 
#            <80% may be a misclassification. We are also assuming potential misclassification if the probability of being 
#            one community type is <10% greater than another community type. Discuss this with your advisor.
###############

### updated summaries
# summary

# histogram
#hist(components$CST_assign_prob, breaks=100, main="Density Plot Probability of Dirichlet Component", xlab="Probability", ylab="Number of Samples")

# boxplot
#boxplot(components$CST_assign_prob)

# boxplot by community type
#plot(maxcomp~as.factor(CST),data=components, main="Probability of Dirichlet Component by Component",xlab="Component", ylab="Probability")

# distribution of composition across all samples

###### Merge community type to metadata (study data) #######

### rename community type probabilities so they make sense - 
# note: this step is up to you. You can just remove the columns if you want. If you do rename, make sure you're assigning
#       new names to the correct columns
components_nicerNames <- mixture(best) %>% 
  as.data.frame() %>% 
  set_names(paste0("CST", 1:opts_args_st2$number_components_bestFit, "prob")) %>% 
  rownames_to_column("ID")

write.csv(components_nicerNames, file = file.path(opts_args$output, "DMM_output_data", "05_CST_assignment_probab.csv"))
write.csv(components %>% rownames_to_column("ID"), file = file.path(opts_args$output, "DMM_output_data", "06_CSTs_ready_for_metadata.csv"))


t0 - Sys.time()
