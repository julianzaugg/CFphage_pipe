# Generate summary figures from viral summary spreadsheet

library(tidyverse)
library(readxl)
library(ggplot2)

setwd("/Users/julianzaugg/Desktop/ACE/major_projects/CF_phage/minion/analysis/")

source("code/utility.R")

colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")

# Load viral summary sheet
viral_summary.df <- readxl::read_excel("data/viral_summary/viral_summary.xlsx") %>% as.data.frame()

# Create run variable
viral_summary.df$Run <- gsub("(.*)_barcode.*", "\\1", viral_summary.df$Sample)

# Assign value of zero for Number_of_genes blank
viral_summary.df$Number_of_genes[is.na(viral_summary.df$Number_of_genes)] <- 0

# Assign values for empty Majority_taxonomy and create Cluster_label
viral_summary.df <-
  viral_summary.df %>%
  mutate(Majority_taxonomy =
           case_when(Majority_taxonomy_level == "Not enough genes with hits to reference"~"Not enough genes with hits to reference",
                     Majority_taxonomy_level == "No majority lineage"~"No majority lineage",
                     is.na(Majority_taxonomy_level) & Number_of_genes == 0~"No genes",
                     TRUE~as.character(.$Majority_taxonomy)
                     )
         ) %>%
  group_by(Cluster) %>%
  mutate(Cluster_size = n(),
         Cluster_label = case_when(Cluster == "Singleton"~paste0(Cluster, "\nn = ", Cluster_size),
                                   TRUE~paste0("Cluster ", gsub("cluster_", "", Cluster), "\nn = ", Cluster_size))) %>%
  as.data.frame()

# Order Majority taxa
taxa <-
  viral_summary.df %>%
  filter(!Majority_taxonomy %in% c("No majority lineage",
                                   "Not enough genes with hits to reference",
                                   "No genes")) %>%
  pull(Majority_taxonomy) %>%
  as.character() %>%
  unique() %>%
  sort()

viral_summary.df$Majority_taxonomy <- factor(viral_summary.df$Majority_taxonomy, levels = c(taxa, c("No majority lineage", "Not enough genes with hits to reference","No genes")))


# Calculate the size (number of members) of each cluster
Cluster_sizes.df <-
  viral_summary.df %>%
  group_by(Cluster) %>%
  summarise(Cluster_size = n()) %>%
  arrange(desc(Cluster_size))

# Calculate the mean sequence length of each cluster
Cluster_mean_sequence_lengths.df <-
  viral_summary.df %>%
  group_by(Cluster) %>%
  summarise(Cluster_mean_sequence_length = mean(Sequence_length)) %>%
  arrange(desc(Cluster_mean_sequence_length))

# Calculate the mean number of genes of each cluster
Cluster_mean_snumber_of_genes.df <-
  viral_summary.df %>%
  group_by(Cluster) %>%
  summarise(Cluster_mean_number_of_genes = mean(Number_of_genes)) %>%
  arrange(desc(Cluster_mean_number_of_genes))

# Factorise the Cluster and label variable by order of size
Cluster_sizes.df$Cluster <- factor(Cluster_sizes.df$Cluster, as.character(Cluster_sizes.df$Cluster))
viral_summary.df$Cluster  <- factor(viral_summary.df$Cluster, as.character(Cluster_sizes.df$Cluster))
viral_summary.df$Cluster_label <- factor(viral_summary.df$Cluster_label, levels = unique(as.character(viral_summary.df$Cluster_label[order(viral_summary.df$Cluster)])))

# Factorise other variables
viral_summary.df$Run <- factor(viral_summary.df$Run, levels = unique(as.character(viral_summary.df$Run[order(as.numeric(gsub(".*batch(.*)", "\\1",viral_summary.df$Run)))])))
viral_summary.df$checkv_quality <- factor(viral_summary.df$checkv_quality, c("Not-determined","Low-quality", "Medium-quality","High-quality", "Complete"))

# ------------------------------------------------------------------------------
# Assign colours to discrete variables

checkv_lq_col <- rgb(red=121,green=168,blue=122,maxColorValue = 255)
checkv_mq_col <- rgb(red=121,green=175,blue=201,maxColorValue = 255)
checkv_hq_col <- rgb(red=255,green=196,blue=0,maxColorValue = 255)
checkv_c_col_col <- rgb(red=179,green=19,blue=19,maxColorValue = 255)

viral_summary.df <- assign_colours_to_df(viral_summary.df,
                     columns = c("Viral_tool",
                                 "Derived_from",
                                 "Run",
                                 "checkv_quality"),
                     auto_assign = T,
                     my_palette = list(checkv_quality = list("Low-quality" = checkv_lq_col,
                                                             "Medium-quality" = checkv_mq_col,
                                                             "High-quality" = checkv_hq_col,
                                                             "Complete" = checkv_c_col_col)),
                     my_default_palette = colour_palette_soft_8
                     )
viral_summary.df <- assign_colours_to_df(viral_summary.df,columns = c("Majority_taxonomy","Sample"))

# Create palettes
tax_palette <- setNames(unique(viral_summary.df$Majority_taxonomy_colour),unique(viral_summary.df$Majority_taxonomy))
tax_palette["No genes"] <- "grey50"
tax_palette["Not enough genes with hits to reference"] <- "grey90"
# tax_palette["No majority lineage"] <- darken(rgb(red=179,green=19,blue=19,maxColorValue = 255))
tax_palette["Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned;Unassigned"] <- "orange"
  # lighten(rgb(red=179,green=19,blue=19,maxColorValue = 255))

viral_tool_palette <- setNames(unique(viral_summary.df$Viral_tool_colour),unique(viral_summary.df$Viral_tool))
derived_from_palette <- setNames(unique(viral_summary.df$Derived_from_colour),unique(viral_summary.df$Derived_from))
sample_palette <- setNames(unique(viral_summary.df$Sample_colour),unique(viral_summary.df$Sample))
run_palette <- setNames(unique(viral_summary.df$Run_colour),unique(viral_summary.df$Run))
checkv_quality_palette <- setNames(unique(viral_summary.df$checkv_quality_colour),unique(viral_summary.df$checkv_quality))
checkv_quality_palette["Not-determined"] <- "grey50"

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Pie charts

# Percentage taxonomy string each cluster
cluster_tax_percentage.df <-
  viral_summary.df %>%
  select(Majority_taxonomy_level, Majority_taxonomy, Cluster, Cluster_size, Cluster_label) %>%
  group_by(Cluster, Majority_taxonomy) %>%
  mutate(Taxonomy_count = n(),
         Taxonomy_percentage = round(n()/Cluster_size,4)*100,
         Taxonomy_percentage_label = case_when(Taxonomy_percentage > 10~paste0(round(Taxonomy_percentage,0), "%"),
                                               TRUE~"")
         ) %>%
  unique() %>%
  as.data.frame()

myplot <-
  ggplot(cluster_tax_percentage.df, aes(x = "",y = Taxonomy_percentage, fill = Majority_taxonomy)) +
  geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
  geom_text(size = 2,aes(label = Taxonomy_percentage_label),
            position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(title = "Composition of viral clusters: Majority taxonomy") +
  coord_polar("y") +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~Cluster_label) +
  scale_fill_manual(values = tax_palette, name = "Majority taxonomy") +
  theme(panel.background = element_blank(),
        legend.title.align = 0.5,
        strip.text = element_text(face = "bold",size = 10),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        )
ggsave(filename = "results/cluster_majority_taxonomy_composition.pdf",
       plot = myplot,
       height = 50,
       width = 70,
       units = "cm",
       device = "pdf")


# Percentage viral tool each cluster
cluster_viral_tool_percentage.df <-
  viral_summary.df %>%
  select(Viral_tool, Cluster, Cluster_size, Cluster_label) %>%
  group_by(Cluster, Viral_tool) %>%
  mutate(Viral_tool_count = n(),
         Viral_tool_percentage = round(n()/Cluster_size,4)*100,
         Viral_tool_percentage_label = case_when(Viral_tool_percentage > 10~paste0(round(Viral_tool_percentage,0), "%"),
                                                 TRUE~"")
  )%>%
  unique() %>%
  as.data.frame()

myplot <-
  ggplot(cluster_viral_tool_percentage.df, aes(x = "",y = Viral_tool_percentage, fill = Viral_tool)) +
  geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
  geom_text(size = 2,aes(label = Viral_tool_percentage_label),
            position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(title = "Composition of viral clusters: Viral tool") +
  coord_polar("y") +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~Cluster_label) +
  scale_fill_manual(values = viral_tool_palette, name = "Viral tool") +
  theme(panel.background = element_blank(),
        legend.title.align = 0.5,
        strip.text = element_text(face = "bold",size = 10),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
  )
ggsave(filename = "results/cluster_viral_tool_composition.pdf",
       plot = myplot,
       height = 50,
       width = 70,
       units = "cm",
       device = "pdf")


# Percentage derived from each cluster
cluster_derived_from_percentage.df <-
  viral_summary.df %>%
  select(Derived_from, Cluster, Cluster_size, Cluster_label) %>%
  group_by(Cluster, Derived_from) %>%
  mutate(Derived_from_count = n(),
         Derived_from_percentage = round(n()/Cluster_size,4)*100,
         Derived_from_percentage_label = case_when(Derived_from_percentage > 10~paste0(round(Derived_from_percentage,0), "%"),
                                                 TRUE~"")
  )%>%
  unique() %>%
  as.data.frame()

myplot <-
  ggplot(cluster_derived_from_percentage.df, aes(x = "",y = Derived_from_percentage, fill = Derived_from)) +
  geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
  geom_text(size = 2,aes(label = Derived_from_percentage_label),
            position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(title = "Composition of viral clusters: Derived from") +
  coord_polar("y") +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~Cluster_label) +
  scale_fill_manual(values = derived_from_palette, name = "Derived from") +
  theme(panel.background = element_blank(),
        legend.title.align = 0.5,
        strip.text = element_text(face = "bold",size = 10),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
  )
ggsave(filename = "results/cluster_derived_from_composition.pdf",
       plot = myplot,
       height = 50,
       width = 70,
       units = "cm",
       device = "pdf")


# Percentage sample each cluster
cluster_sample_percentage.df <-
  viral_summary.df %>%
  select(Sample, Cluster, Cluster_size, Cluster_label) %>%
  group_by(Cluster, Sample) %>%
  mutate(Sample_count = n(),
         Sample_percentage = round(n()/Cluster_size,4)*100,
         Sample_percentage_label = case_when(Sample_percentage > 10~paste0(round(Sample_percentage,0), "%"),
                                             TRUE~"")
  )%>%
  unique() %>%
  as.data.frame()

myplot <-
  ggplot(cluster_sample_percentage.df, aes(x = "",y = Sample_percentage, fill = Sample)) +
  geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
  geom_text(size = 2,aes(label = Sample_percentage_label),
            position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(title = "Composition of viral clusters: Sample") +
  coord_polar("y") +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~Cluster_label) +
  scale_fill_manual(values = sample_palette, name = "Sample") +
  theme(panel.background = element_blank(),
        legend.title.align = 0.5,
        strip.text = element_text(face = "bold",size = 10),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
  )
ggsave(filename = "results/cluster_sample_composition.pdf",
       plot = myplot,
       height = 50,
       width = 70,
       units = "cm",
       device = "pdf")


# Percentage run each cluster
cluster_run_percentage.df <-
  viral_summary.df %>%
  select(Run, Cluster, Cluster_size, Cluster_label) %>%
  group_by(Cluster, Run) %>%
  mutate(Run_count = n(),
         Run_percentage = round(n()/Cluster_size,4)*100,
         Run_percentage_label = case_when(Run_percentage > 10~paste0(round(Run_percentage,0), "%"),
                                          TRUE~"")
  )%>%
  unique() %>%
  as.data.frame()

myplot <-
  ggplot(cluster_run_percentage.df, aes(x = "",y = Run_percentage, fill = Run)) +
  geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
  geom_text(size = 2,aes(label = Run_percentage_label),
            position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(title = "Composition of viral clusters: Run") +
  coord_polar("y") +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~Cluster_label) +
  scale_fill_manual(values = run_palette, name = "Run") +
  theme(panel.background = element_blank(),
        legend.title.align = 0.5,
        strip.text = element_text(face = "bold",size = 10),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
  )
ggsave(filename = "results/cluster_run_composition.pdf",
       plot = myplot,
       height = 50,
       width = 70,
       units = "cm",
       device = "pdf")


# Percentage checkv_quality each cluster
cluster_checkv_quality_percentage.df <-
  viral_summary.df %>%
  select(checkv_quality, Cluster, Cluster_size, Cluster_label) %>%
  group_by(Cluster, checkv_quality) %>%
  mutate(checkv_quality_count = n(),
         checkv_quality_percentage = round(n()/Cluster_size,4)*100,
         checkv_quality_percentage_label = case_when(checkv_quality_percentage > 10~paste0(round(checkv_quality_percentage,0), "%"),
                                                     TRUE~"")
  )%>%
  unique() %>%
  as.data.frame()

myplot <-
  ggplot(cluster_checkv_quality_percentage.df, aes(x = "",y = checkv_quality_percentage, fill = checkv_quality)) +
  geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
  geom_text(size = 2,aes(label = checkv_quality_percentage_label),
            position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(title = "Composition of viral clusters: Checkv quality") +
  coord_polar("y") +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~Cluster_label) +
  scale_fill_manual(values = checkv_quality_palette, name = "Checkv quality") +
  theme(panel.background = element_blank(),
        legend.title.align = 0.5,
        strip.text = element_text(face = "bold",size = 10),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
  )
ggsave(filename = "results/cluster_checkv_quality_composition.pdf",
       plot = myplot,
       height = 50,
       width = 70,
       units = "cm",
       device = "pdf")


# temp <- viral_summary.df %>%
#   mutate(Viral_tool_checkv_quality = paste0(Viral_tool, ":",checkv_quality)) %>%
#   select(Viral_tool_checkv_quality,Viral_tool, checkv_quality, Cluster, Cluster_size, Cluster_label) %>%
#   group_by(Cluster, Viral_tool_checkv_quality) %>%
#   mutate(checkv_quality_count = n(),
#          checkv_quality_percentage = round(n()/Cluster_size,4)*100,
#          checkv_quality_percentage_label = case_when(checkv_quality_percentage > 10~paste0(round(checkv_quality_percentage,0), "%"),
#                                                      TRUE~"")
#   )%>%
#   unique() %>%
#   as.data.frame()
#
# ggplot(temp, aes(x = "",y = checkv_quality_percentage, fill = Viral_tool_checkv_quality)) +
#   geom_bar(stat = "identity", colour = "grey20", lwd = .1) +
#   geom_text(size = 2,aes(label = checkv_quality_percentage_label),
#             position = position_stack(vjust = 0.5)) +
#   xlab("") +
#   ylab("") +
#   labs(title = "Composition of viral clusters: checkv_quality") +
#   coord_polar("y") +
#   scale_x_discrete(drop=FALSE) +
#   facet_wrap(~Cluster_label) +
#   # scale_fill_manual(values = checkv_quality_palette, name = "checkv_quality") +
#   theme(panel.background = element_blank(),
#         legend.title.align = 0.5,
#         strip.text = element_text(face = "bold",size = 10),
#         axis.line = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
#   )
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Box plots

viral_summary.df$Cluster <- factor(viral_summary.df$Cluster, rev(Cluster_mean_sequence_lengths.df$Cluster))
myplot <-
  ggplot(viral_summary.df,aes(x = Cluster, y = Sequence_length)) +
  geom_boxplot(outlier.shape = NA, fill = "grey20", colour = "grey40", lwd = .2) +
  geom_jitter(size = 1,width = .1, shape = 21, aes(fill = Majority_taxonomy, colour = Majority_taxonomy), alpha= .6) +
  ylab("Sequence length") +
  # coord_flip() +
  scale_fill_manual(values = tax_palette, name = "Majority taxonomy") +
  scale_colour_manual(values = lapply(tax_palette, darken,3), name = "Majority taxonomy") +
  scale_y_continuous(breaks = seq(0, 115000, 5000), limits= c(0,115000),expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, colour = "black"),
        legend.title = element_text(face = "bold", size = 10, hjust = 0.5),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.box.background = element_blank(),
        panel.grid.major.y = element_line(colour="grey70", size=0.1),
        panel.grid.major.x = element_blank()
        )
ggsave(filename = "results/cluster_sequence_lengths.pdf",
       plot = myplot,
       height = 12,
       width = 100,
       units = "cm",
       device = "pdf")

viral_summary.df$Cluster <- factor(viral_summary.df$Cluster, rev(Cluster_mean_snumber_of_genes.df$Cluster))
myplot <-
  ggplot(viral_summary.df,aes(x = Cluster, y = Number_of_genes)) +
  geom_boxplot(outlier.shape = NA, fill = "grey20", colour = "grey40", lwd = .2) +
  geom_jitter(size = 1,width = .1, shape = 21, aes(fill = Majority_taxonomy, colour = Majority_taxonomy), alpha= .6) +
  ylab("Number of genes") +
  # coord_flip() +
  scale_fill_manual(values = tax_palette, name = "Majority taxonomy") +
  scale_colour_manual(values = lapply(tax_palette, darken,3), name = "Majority taxonomy") +
  scale_y_continuous(breaks = seq(0, 200, 20), limits= c(-1,201),expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, colour = "black"),
        legend.title = element_text(face = "bold", size = 10, hjust = 0.5),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.box.background = element_blank(),
        panel.grid.major.y = element_line(colour="grey70", size=0.1),
        panel.grid.major.x = element_blank()
  )
ggsave(filename = "results/cluster_number_of_genes.pdf",
       plot = myplot,
       height = 10,
       width = 100,
       units = "cm",
       device = "pdf")



