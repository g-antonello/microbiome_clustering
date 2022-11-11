
t0 <- Sys.time()

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=FALSE, help = "Directory where to find step 1 and 2 output data")
)


opt_parser = OptionParser(option_list=option_list)

opts_args_st3 = parse_args(opt_parser)

print(opts_args_st3)

file_input0 <- list.files(path = opts_args_st3$input, pattern = "04_taxon", full.names = TRUE)

taxcontribution = read.csv(file_input0, header = TRUE)

n_communities <- colnames(taxcontribution) %>% str_extract("k[2-9]_") %>% table() %>% names() %>% parse_number()

data_for_stackplots <- taxcontribution %>% 
  select(starts_with(c("microbial", paste0("k", n_communities)))) %>% 
  pivot_longer(cols = starts_with(paste0("k", n_communities)), names_to = "CST") %>% 
  mutate(microbial_feature = case_when(value < 0.01 ~ "Other", 
                                        TRUE ~ microbial_feature), 
         CST = paste0("CST ", gsub(paste0("k", n_communities, "_com"), "", CST))) %>% 
  group_by(microbial_feature, CST) %>% 
  summarise(value = sum(value))

set.seed(2993)
microbial_feature_colors <- randomColor(length(unique(data_for_stackplots$microbial_feature)))

stackplots_CST_mean_composition <- ggplot(data_for_stackplots, aes(x = CST, y = value, fill = microbial_feature))+ 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = microbial_feature_colors)+ 
  labs(x = NULL, 
       y = "Taxon Mean Relative abundance",
       fill = NULL)+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=4,
                           byrow=TRUE))


ggsave(plot = stackplots_CST_mean_composition, 
       filename = file.path(opts_args_st3$input, "../DMM_output_figures/03_MeanTaxContrib_per_CST.pdf"),
       height = 7, 
       width = 9
       )



################################################################################
###############################  WAFFLE PLOTS ##################################
################################################################################

others_row <- data_for_stackplots %>% 
  mutate(value = round(value*100, 0)) %>% 
  filter(microbial_feature != "Other") %>% 
  group_by(CST) %>% 
  summarise(tot = 100 - sum(value)) %>% 
  .$tot

others_df <- data.frame(microbial_feature = rep("Other", length(unique(data_for_stackplots$CST))),
                        CST = unique(data_for_stackplots$CST),
                        value = others_row)

tmp.list <- data_for_stackplots %>% 
  mutate(value = round(value*100, 0)) %>% 
  filter(microbial_feature != "Other") %>% 
  rbind.data.frame(others_df) %>% 
  split.data.frame(x = ., f = .$CST) %>% 
  lapply(select, -CST) 

all_ASVs <- lapply(tmp.list, "[[", "microbial_feature") %>% 
  Reduce(union, .) %>% 
  sort()

tmp.list2 <- list()

for(i in 1:length(tmp.list)){
  
  tmp.list2[[names(tmp.list)[[i]]]] <- ifelse(is.na(match(all_ASVs, tmp.list[[i]]$microbial_feature)), 
                                              yes = 0, 
                                              no = tmp.list[[i]]$value[match(all_ASVs, tmp.list[[i]]$microbial_feature)]
                                              )
  names(tmp.list2[[names(tmp.list)[[i]]]]) <- all_ASVs
}


waffle_colors <- randomcoloR::distinctColorPalette(length(all_ASVs))

# Waffle plots

waffle_plots_list <- sapply(names(tmp.list2), function(n) 
  suppressMessages(
    waffle::waffle(tmp.list2[[n]], 
                 colors = waffle_colors, 
                 title = n, 
                 xlab = NULL,
                 flip = TRUE,
                 reverse = T,
                 size = 1,
                 keep = T)+
      guides(fill=guide_legend(nrow=3,
                                   byrow=TRUE))
  #  scale_fill_manual(breaks = all_ASVs, values = waffle_colors)
    ),
  USE.NAMES = TRUE,
  simplify = FALSE
  )

plots_arranged <- ggpubr::ggarrange(plotlist = waffle_plots_list, common.legend = T) %>% 
  ggpubr::annotate_figure(bottom = "1 tile is 1% relative abundance")

ggsave(plot = plots_arranged,
       filename = file.path(opts_args_st3$input, "../DMM_output_figures/04_Waffle_plots_CST_composition.pdf"),
       height = 8, 
       width = 11)


