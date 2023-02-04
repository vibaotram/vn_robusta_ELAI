#!/usr/bin/env Rscript

### scripts to estimate ancestral proportion of all the African and VN individuals

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)
library(StAMPP)
library(huxtable)

# accessions information
seq_info_file <- "sequence_info_final.tsv"
seq_info <- read_tsv(seq_info_file)

# final vcf file containing 1.1M SNPs of all individuals
allvcf_file <- "/data3/projects/vietcaf/baotram/vietcaf_final_biallelic_random_10perc.vcf"

# extract genotypes and convert it to geno format on LEA package
all_vcfR <- read.vcfR(allvcf_file)
all_vcfR@fix[,1] <- gsub(".", "_", all_vcfR@fix[,1], fixed = T)
all_gl <- vcfR2genlight(all_vcfR, n.cores = 20)
all_geno <- as.matrix(all_gl)
all_geno[is.na(all_geno)] <- 9

# run snmf for only african
snmf_dir <- "snmf_random_10perc_SNPs_AF"
dir.create(snmf_dir)
af_geno_file <- file.path(snmf_dir, "random_10perc_SNPs_AF.geno")
write.geno(all_geno[!grepl("TR", rownames(all_geno)),], af_geno_file)

af_snmf <- snmf(af_geno_file, K = 1:10, CPU = 20, repetitions = 10, project = "new", entropy = T, seed = 987654321)

k <- 5
best_run <- which.min(cross.entropy(all_snmf, K = k))

## plot ancestry proportion
ind_order <- barchart(all_snmf, K = k, run = best_run, plot = F)
prop <- as.data.frame(Q(all_snmf, K = k, run = best_run))
prop$Label <- grep("TR", rownames(all_geno), value = T, invert = T)


## snmf for all inds

snmf_dir <- "snmf_random_10perc_SNPs"
dir.create(snmf_dir)
all_geno_file <- file.path(snmf_dir, "random_10perc_SNPs.geno")
write.geno(all_geno, all_geno_file)

all_snmf <- snmf(all_geno_file, K = 1:10, CPU = 20, repetitions = 10, project = "new", entropy = T, seed = 987654321)

plot_dir <- "../plots/genetic_structure"

tiff(file = file.path(plot_dir, "cross_entropy.tiff"), 
           width = 4, height = 2, units = "in", res = 600)
par(mar = c(5,5,2,2), cex = 0.5)
plot(all_snmf, pch = 16)
dev.off()

## best k = 5
## choose best run based on cross entropy
k <- 5
best_run <- which.min(cross.entropy(all_snmf, K = k))

## plot ancestry proportion
# ind_order <- barchart(all_snmf, K = k, run = best_run, plot = F)
prop_mt <- as.data.frame(Q(all_snmf, K = k, run = best_run))
prop_mt$Label <- seq_info$Label
colnames(prop_mt)[1:5] <- c("C", "D", "OB", "AG", "ER")
# label_ord <- prop_mt %>% 
#   rowwise() %>%
#   mutate(group = colnames(prop_mt)[which.max(c(C, D, OB, AG, ER))]) %>% 
#   # arrange(group, AG, C, ER, D, OB) %>% 
#   arrange(group, ER, AG, D, C, OB) %>% 
#   pull(Label)

prop_mt <- prop_mt %>% 
  rowwise() %>%
  mutate(group = colnames(prop_mt)[which.max(c(C, D, OB, AG, ER))]) 


prop <- reshape2::melt(prop_mt, variable.name = "ancestry", value.name = "proportion")
prop <- left_join(prop, seq_info, by = "Label")

## reorder the individuals
labelnotOB_ord <- prop %>% 
  group_by(group) %>% 
  filter(ancestry == group) %>% 
  filter(group != "OB") %>% 
  mutate(group = factor(group, levels = c("D", "C", "AG", "ER"))) %>% 
  arrange(group, desc(proportion)) %>% 
  pull(Label)

labelOB_ord <- prop %>% 
  group_by(group) %>% 
  filter(ancestry == group) %>% 
  filter(group == "OB") %>% 
  arrange(proportion) %>% 
  pull(Label)

label_ord <- c(labelnotOB_ord, labelOB_ord)
# prop$Label <- factor(prop$Label, levels = unique(prop$Label)[ind_order$order], ordered = F)
# levels(prop$ancestry) <- c("C", "D", "OB", "AG", "ER")
# prop$ancestry <- factor(prop$ancestry, levels = c("ER", "OB", "C", "AG", "D"))
# # prop$Country <- factor(prop$Country, levels = c("CI", "GH", "GN", "CF", "CM", "CR", "GA", "TG", "ID", "AO", "UG", "VN"))
# prop$Country <- factor(prop$Country, levels = c("Ivory Coast", "Ghana", "Guinea", "RCA", "Cameroon", "Congo", "RDC", "Gabon", "Togo", "Angola", "Uganda", "Indonesia", "Vietnam"))

group_col <- c("#FB8022FF", "#A2FC3CFF", "#7A0403FF", "#30123BFF", "#28BBECFF")
names(group_col) <- c("AG", "C", "D", "ER", "OB")

# tiff(file = file.path(plot_dir, "population_structure_k5.tiff"), 
#     # colormodel = "rgb", 
#     width = 12, height = 5, units = "in", res = 1200)
# ggplot(prop) +
#   facet_grid(cols = vars(Country), scales = "free", space = "free", switch = "x") +
#   geom_col(aes(x = Label, y = proportion, fill = ancestry), position = position_stack(), width = .95, size = 0) + 
#   # scale_fill_viridis_d(end = .9)+
#   scale_fill_manual(values = group_col) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   ylab("Ancestry proportion")  +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_blank(),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.key.size = unit(5, "mm"),
#         axis.ticks.y = element_line(size = .1),
#         axis.text.y = element_text(size = 8),
#         panel.spacing = unit(2, "mm"),
#         strip.text.x = element_text(angle = 90, size = 10, hjust = 1)) 
# dev.off()


tiff(file = file.path(plot_dir, "population_structure_k5_v2_rv.tiff"), 
    # colormodel = "rgb", 
    width = 8, height = 3, units = "in", res = 1200)
prop %>% 
  mutate(type = ifelse(Country == "Vietnam", "Vietnam", "Reference"),
         Label = fct_relevel(Label, label_ord), 
         ancestry = fct_relevel(ancestry, c("D", "C", "AG", "ER", "OB"))) %>% 
  ggplot() +
  facet_grid(cols = vars(type), scales = "free", space = "free", switch = "x") +
  geom_col(aes(x = Label, y = proportion, fill = ancestry), position = position_stack(), width = .95, size = 0) + 
  # scale_fill_viridis_d(end = .9)+
  scale_fill_manual(values = group_col) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Ancestry proportion")  +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(5, "mm"),
        axis.ticks.y = element_line(linewidth = .1),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(3, "mm"),
        panel.border = element_rect(linewidth = 2, color = "gray60", fill = NA),
        # panel.background = element_blank(),
        strip.text.x = element_text(size = 12)) 
dev.off()


## assign groups based on snmf results

prop <- prop %>% 
  group_by(Label) %>% 
  arrange(desc(proportion)) %>%
  slice_head(n = 1) %>% 
  rowwise() %>% 
  mutate(group = ifelse(proportion > 0.8, as.character(ancestry), "hybrid")) %>% 
  mutate(group = ifelse(grepl("VN", as.character(Label)), "Vietnam", group))

all_gl@pop <- as.factor(prop$group)

group_gl <- all_gl[all_gl@pop != "hybrid"]
group_fst <- stamppFst(group_gl)


matrix_file <- tempfile(fileext = ".lfmm")
write.lfmm(all_geno, matrix_file)
all_pca <- pca(matrix_file)

plot(all_pca$projections)


## fst
source("./group_fst.R")
group_fst <- readRDS("./group_fst_1M_SNPs.RDS")

group_fst$Fsts[-1,-6] %>% 
  as.data.frame() %>% 
  mutate(across(OB:AG, ~ ifelse(is.na(.x), "", sprintf(.x, fmt = "%.2f")))) %>% 
  rownames_to_column("Group") %>%
  hux() %>% 
  set_right_border(everywhere, 1, brdr(3, "solid", "grey")) %>% 
  set_bottom_border(1, everywhere, brdr(3, "solid", "grey")) %>% 
  set_bold(row = 1, col = everywhere) %>% 
  set_bold(row = everywhere, col = 1) %>% 
  quick_html(.)


## PCA 
all_pca <- pca(geno2lfmm(all_geno_file))

all_pca_proj <- all_pca$projections[,1:2]
all_pca_proj <- all_pca_proj %>% 
  as.data.frame() %>% 
  mutate(Label = prop_mt$Label,
         group = prop_mt$group) %>% 
  mutate(ancestry = if_else(grepl("TR", Label), "Vietnam", group))

variance <- all_pca$eigenvalues/sum(all_pca$eigenvalues)

group_col <- c("#FB8022FF", "#A2FC3CFF", "#7A0403FF", "#30123BFF", "#28BBECFF", "gray50")
names(group_col) <- c("AG", "C", "D", "ER", "OB", "Vietnam")

tiff(file = file.path(plot_dir, "pca12.tiff"), 
     width = 6, height = 5, units = "in", res = 1000)
ggplot(all_pca_proj) +
  geom_point(aes(x = V1, y = V2, color = ancestry, shape = ancestry), 
             size = 5, alpha = 0.8) +
  xlab(paste0("PC1 (", sprintf("%.1f", variance[1]*100), "%)")) + 
  ylab(paste0("PC2 (", sprintf("%.1f", variance[2]*100), "%)")) + 
  scale_color_manual(values = group_col) + 
  scale_shape_manual(values = c(rep(19, 5), 6)) +
  guides(color = guide_legend(title="group"),
         shape = guide_legend(title="group")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 5),
        panel.grid.minor = element_blank()) 
dev.off()
