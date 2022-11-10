
#########################################################################################################
#                                  MERGE RESULTS IN PHYLOSEQ OF INTEREST (NAFLD)                        #
#########################################################################################################
#add to metadata
chrismb_phy2 <- readRDS("~/CHRISMB/saliva_overview/Data analysis/Microbiome_data/02_chrismb_phyloseq_newBetsyMetadata.rds")

m <- microbiome::meta(chrismb_phy)
components$michigan_codes <- rownames(components)
components <- components %>% select(c("michigan_codes", "CST", "CST_assign_prob"))
m_merged_CST <- merge(m, components, by = "michigan_codes")
m_merged_CST_sorted <- m_merged_CST[match(m$michigan_codes,m_merged_CST$michigan_codes),]
rownames(m_merged_CST_sorted) <- m_merged_CST_sorted$michigan_codes

sample_data(chrismb_phy) <- sample_data(m_merged_CST_sorted)

saveRDS(chrismb_phy, "~/CHRISMB/saliva_overview/Data analysis/Microbiome_data/02_chrismb_phyloseq_newBetsyMetadata.rds")
#########################################################################################################

```
## Looking at distribution of taxa in each CST

```{r}

#read.csv("project_genus_long.csv")

# library(reshape2)
# library(ggplot2)
taxcontribution <- read.csv(file="../output/DMM_output_tables/taxon_contribution_k6_model.csv", row.names = 1 ) 


library(ggplot2)
library(reshape2)
library(dplyr)
taxcontribution$genus<-rownames(taxcontribution)
taxcontribution$genus<-ifelse(taxcontribution$differ_k1_k6_coeff > 0.90, "Other", taxcontribution$genus) #differ_k1_k6_coeff is cumulative contribution, those with a higher differ_k1_k6_coeff means they contribute less (can look at data if confused, see how differ_k1_k6_coeff is a cumulative sum as you go down the rows)
taxcontribution <- taxcontribution %>% group_by(genus) %>% summarise(m1 = sum(k6_com1), m2=sum(k6_com2), m3=sum(k6_com3), m4=sum(k6_com4), m5 = sum(k6_com5), m6 = sum(k6_com6))

taxcontribution1 <- melt(taxcontribution, id = "genus")

#get colors for plotting
genus_distinct <- length(unique(taxcontribution1$genus)) 
library(randomcoloR)
genus_colors <- randomColor(genus_distinct) 

CST.tax <- ggplot(taxcontribution1, aes(x=variable, y = value, fill=genus)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = genus_colors)+
  xlab("Community State Type") + ylab("Proportion of Total Community Composition") + 
  theme(legend.position = "bottom") + theme(legend.text = element_text(face = "italic"))+theme(legend.text=element_text(size=8))+
  theme(legend.key.size = unit(.3, "cm"))+ggtitle("CST relative abundance of genera")#+theme(plot.margin = unit(c(.5,3,1,.1), "cm"))
print(CST.tax)
ggsave("../output/CST_genus compositions.png", dpi = 600, height = 10, width = 12)

# interested in seeing which genus varies across CST? see output below this chunk
library(plotly)
ggplotly(CST.tax)

```

## Waffle plot
```{r}
# first make data backup, because this chunk overwrites on the same object
df_backup <- df

smart_round <- function(x, digits = 0) { # somewhere on SO
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

#make the Relative Abundances into percents (only of CST, not the last 2 cols)
df <- df[,2:7]
df <- df*100

df <- smart_round(df)

#remove columns with 0
df <- df[rowSums(df)>0,]
#separate into different CST types
CST1 <-df[,1]
names(CST1) <- rownames(df)
CST2 <-df[,2]
names(CST2) <- rownames(df)
CST3 <-df[,3]
names(CST3) <- rownames(df)
CST4 <-df[,4]
names(CST4) <- rownames(df)
CST5 <-df[,5]
names(CST5) <- rownames(df)
CST6 <-df[,6]
names(CST6) <- rownames(df)


library(RColorBrewer)
library(grDevices)

#set diverse colors for ease of distinction
waffle_colors <- randomcoloR::distinctColorPalette(length(CST1))

# Waffle plots
CST1_waffle <- waffle::waffle(CST1, colors = waffle_colors, title = "CST1", xlab ="1 tile is 1% relative abundance")
ggsave(plot = CST1_waffle, filename = "../output/CST1_wafflePlot.png", dpi =600, height = 10, width = 10)

CST2_waffle <- waffle::waffle(CST2, colors = waffle_colors, title = "CST2", xlab ="1 tile is 1% relative abundance")
ggsave(plot = CST2_waffle, filename = "../output/CST2_wafflePlot.png", dpi =600, height = 10, width = 10)

CST3_waffle <- waffle::waffle(CST3, colors = waffle_colors, title = "CST3", xlab ="1 tile is 1% relative abundance")
ggsave(plot = CST3_waffle, filename = "../output/CST3_wafflePlot.png", dpi =600, height = 10, width = 10)

CST4_waffle <- waffle::waffle(CST4, colors = waffle_colors, title = "CST4", xlab ="1 tile is 1% relative abundance")
ggsave(plot = CST4_waffle, filename = "../output/CST4_wafflePlot.png", dpi =600, height = 10, width = 10)

CST5_waffle <- waffle::waffle(CST5, colors = waffle_colors, title = "CST5", xlab ="1 tile is 1% relative abundance")
ggsave(plot = CST5_waffle, filename = "../output/CST5_wafflePlot.png", dpi =600, height = 10, width = 10)

CST6_waffle <- waffle::waffle(CST6, colors = waffle_colors, title = "CST6", xlab ="1 tile is 1% relative abundance")
ggsave(plot = CST6_waffle, filename = "../output/CST6_wafflePlot.png", dpi =600, height = 10, width = 10)


library(ggpubr)
waffles_combined <- ggarrange(CST1_waffle, CST2_waffle, CST3_waffle, CST4_waffle, CST5_waffle,CST6_waffle, ncol = 3, nrow =2, legend = "bottom", common.legend = TRUE)
print(waffles_combined)
ggsave(plot = waffles_combined, filename = "../output/wafflePlots_combined.png", dpi = 600, height = 10, width = 10)
```



###### Descriptive Analysis ######
### check frequency of community types by variables

```{r}

#Let's say I wanna see how the CST distribution looks within my variable "Fake" (which has two groups: Sun and No Sun) 

#comp is what CST they fit to
#don't forget some samples weren't assigned to a CST
table(meta_merged_CST$comp, useNA="always")

ggplot(meta_merged_CST, aes(x = x0bp01, fill = as.factor(CST))) + geom_bar() + xlab("systolic blood pressure (mmHg)")
ggsave(filename = "../output/CST_syst_bp_histogram.png", dpi = 600, width = 8, height = 6)
ggplot(meta_merged_CST, aes(x = x0bp02, fill = as.factor(CST))) + geom_bar()+ xlab("diastolic blood pressure (mmHg)")
ggsave(filename = "../output/CST_diast_bp_histogram.png", dpi = 600, width = 8, height = 6)

#get proportions
meta_merged_CST<- mutate(meta_merged_CST, prop = n / sum(n))

#ggplot(meta_merged_CST, aes(x = x0bp01, fill = as.factor(CST))) + geom_bar(aes(y = (..count..)/sum(..count..)), position = "fill" ) + ylab("Proportion")

```

# Preparation 100 most discriminative taxa

Merged with metadata to make heritability analysis on the first 100 genera used to make CSTs

```{r}
library(phyloseq)
first_100_taxa <- rownames(taxcontribution)[1:100]
chrismb_100discrimTaxa <- phyloseq::subset_taxa(chrismb_phy_CST_genus, taxa_names(chrismb_phy_CST_genus) %in% first_100_taxa)
x <- data.frame(t(otu_table(chrismb_100discrimTaxa)))
x$michigan_codes <- rownames(x)

#take metadata
shrunk_metadata <- microbiome::meta(chrismb_100discrimTaxa)[, c("gender", "age", "michigan_codes", "AID")]
names(shrunk_metadata) <- c("sex", "age", "michigan_codes","AID")

shrunk_metadata$AID <- paste("00", shrunk_metadata$AID, sep = "")
#merge the two datasets
merged_meta_100taxa <- merge(shrunk_metadata,x, by = "michigan_codes")
merged_meta_100taxa$michigan_codes <- NULL
merged_meta_100taxa <- merged_meta_100taxa[, c(3,1,2,4:ncol(merged_meta_100taxa))]

write.table(merged_meta_100taxa, file = "../output/100_Influential_taxa_Age_Sex.txt", quote = F, sep = "\t", row.names = FALSE)
```

# Alpha diversity of CSTs

And significance! This link is gold to learn how to plot significance on ggplots:
  http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  
  ```{r}
theme_set(theme_bw())
library(phyloseq)
library(ggplot2)
sample_data(chrismb_phy) <- sample_data(meta)

library(ggpubr)
alpha_div <- estimate_richness(chrismb_phy, measures = "Shannon")
alpha_div$michigan_codes <- rownames(alpha_div)
tmp <- microbiome::meta(chrismb_phy_CST)
tmp_meta <- merge(tmp, alpha_div, by = "michigan_codes")
tmp_meta$CST <- as.character(tmp_meta$CST)

my_stat_comparisons <- combn(sort(unique(tmp_meta$CST)), 2, simplify = FALSE)

ggplot(data = tmp_meta, aes(x = CST, y = Shannon, color = CST)) + geom_boxplot() + stat_compare_means(comparisons = my_stat_comparisons, method = "wilcox.test")

ggsave("../.output/Shannon_diversity_community_types.png", dpi = 600, height = 5, width = 6)


```
