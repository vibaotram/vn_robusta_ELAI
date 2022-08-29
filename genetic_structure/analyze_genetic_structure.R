#!/usr/bin/env Rscript

### scripts to estimate ancestral proportion of all the African and VN individuals

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)

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
ind_order <- barchart(all_snmf, K = k, run = best_run, plot = F)
prop <- as.data.frame(Q(all_snmf, K = k, run = best_run))
prop$Label <- seq_info$Label
prop <- reshape2::melt(prop, variable.name = "ancestry", value.name = "proportion")
prop <- left_join(prop, seq_info, by = "Label")

prop$Label <- factor(prop$Label, levels = unique(prop$Label)[ind_order$order], ordered = F)
levels(prop$ancestry) <- c("C", "D", "OB", "AG", "ER")
prop$ancestry <- factor(prop$ancestry, levels = c("ER", "OB", "C", "AG", "D"))
# prop$Country <- factor(prop$Country, levels = c("CI", "GH", "GN", "CF", "CM", "CR", "GA", "TG", "ID", "AO", "UG", "VN"))
prop$Country <- factor(prop$Country, levels = c("Ivory Coast", "Ghana", "Guinea", "RCA", "Cameroon", "Congo", "RDC", "Gabon", "Togo", "Angola", "Uganda", "Indonesia", "Vietnam"))

group_col <- c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF", "#7A0403FF")

tiff(file = file.path(plot_dir, "population_structure_k5.tiff"), 
    # colormodel = "rgb", 
    width = 12, height = 5, units = "in", res = 1200)
ggplot(prop) +
  facet_grid(cols = vars(Country), scales = "free", space = "free", switch = "x") +
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
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(5, "mm"),
        axis.ticks.y = element_line(size = .1),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(2, "mm"),
        strip.text.x = element_text(angle = 90, size = 10, hjust = 1)) 
dev.off()


tiff(file = file.path(plot_dir, "population_structure_k5_v2.tiff"), 
    # colormodel = "rgb", 
    width = 8, height = 3, units = "in", res = 1200)
prop %>% 
  mutate(type = ifelse(Country == "Vietnam", "Vietnam", "Reference")) %>% 
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
        axis.ticks.y = element_line(size = .1),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(3, "mm"),
        panel.border = element_rect(size = 2, color = "gray60", fill = NA),
        # panel.background = element_blank(),
        strip.text.x = element_text(size = 12)) 
dev.off()
