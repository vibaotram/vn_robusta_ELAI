#!/usr/bin/env Rscript

### simulate ancestral genotypes using 100k random SNPs on chromosome 1

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)

## accessions information
seq_info_file <- "../../genetic_structure/sequence_info_final.tsv"
seq_info <- read_tsv(seq_info_file)


## convert vcf to geno format
chr1_vcf_file <- "vietcaf_final_chr01.recode.vcf"
chr1_vcfR <- read.vcfR(chr1_vcf_file)
chr1_vcfR@fix[,1] <- gsub(".", "_", chr1_vcfR@fix[,1], fixed = T)
chr1_gl <- vcfR2genlight(chr1_vcfR, n.cores = 18)
chr1_geno <- as.matrix(chr1_gl)
chr1_geno[is.na(chr1_geno)] <- 9

## run snmf with 100k random SNPs
k <- 5
test_snps <- sample(locNames(chr1_gl), 100000)
small_gl <- chr1_gl[,test_snps]
small_geno <- as.matrix(small_gl)
small_geno[is.na(small_geno)] <- NA

snmf_dir <- "snmf_chr1_100kSNPs"
dir.create(snmf_dir, showWarnings = F)
test_geno_file <- file.path(snmf_dir, "chr1_100kSNPs.geno")
write.geno(small_geno, test_geno_file)
small_snmf <- snmf(test_geno_file, K = k,
                   CPU = 18, repetitions = 10,
                   project = "new", entropy = T, seed = 987654321)

best_run <- which.min(cross.entropy(small_snmf, K = k))

ind_order <- barchart(small_snmf, K = k, run = best_run, plot = F)
prop_small <- as.data.frame(Q(small_snmf, K = k, run = best_run))
prop_small$Label <- seq_info$Label
prop_small <- reshape2::melt(prop_small, variable.name = "ancestry", value.name = "proportion")
prop_small <- left_join(prop_small, seq_info, by = "Label")

prop_small$code <- factor(prop_small$code, levels = unique(prop_small$code)[ind_order$order], ordered = F)
prop_small$Label <- factor(prop_small$Label, levels = unique(prop_small$Label)[ind_order$order], ordered = F)
levels(prop_small$ancestry) <- c("C", "AG", "D", "ER", "OB")
prop_small$ancestry <- factor(prop_small$ancestry, levels = c("ER", "OB", "C", "D", "AG"))


## get genotypic frequencies
small_freq <- G(small_snmf, K = 5, run = best_run)
n_inds <- 100
n_snps <- nrow(small_freq)/3

## simulate 100 genotype each group 
sml_geno <- data.frame() # genotype matrix of all groups
for (v in 1:ncol(small_freq)) { # for each group
  cat(v)
  freq <- small_freq[,v] # get allele freq of the group
  geno <- matrix(nrow=n_inds, ncol=n_snps) # genotype matrix of the group
  # for (i in 1:n_snps) { 
  #   # for locus i, frequencies of 0, 1, 2 genotype are in row 3i-2, 3i-1, 3i of freq matrix
  #   for (j in 1:n_inds) {
  #     geno[j,i] = sample(c(0,1,2), 1, prob = c(freq[3*i-2], freq[3*i-1], freq[3*i]))
  #   }
  # }
  rownames(geno) <- paste0("simulated_group_", v, "-", 1:n_inds) # name the generated indv
  colnames(geno) <- small_gl@loc.names # SNP IDs
  sml_geno <- rbind(sml_geno, geno)
}

## snmf with only simulated genotypes
sml_geno_file <- file.path(snmf_dir, "sml_chr1_100kSNPs.geno")
write.geno(sml_geno, sml_geno_file)
sml_snmf <- snmf(sml_geno_file, K = k,
                 CPU = 18, repetitions = 10,
                 project = "new", entropy = T, seed = 987654321)

best_run <- which.min(cross.entropy(sml_snmf, K = k))

ind_order <- barchart(sml_snmf, K = k, run = best_run, plot = F)
prop_sml <- as.data.frame(Q(sml_snmf, K = k, run = best_run))
prop_sml$Label <- rownames(sml_geno)
prop_sml <- reshape2::melt(prop_sml, variable.name = "ancestry", value.name = "proportion")
prop_sml <- left_join(prop_sml, seq_info, by = "Label")

# prop_sml$code <- factor(prop_sml$code, levels = unique(prop_sml$code)[ind_order$order], ordered = F)
prop_sml$Label <- factor(prop_sml$Label, levels = unique(prop_sml$Label)[ind_order$order], ordered = F)
levels(prop_sml$ancestry) <- c("C", "AG", "ER", "OB", "D")
prop_sml$ancestry <- factor(prop_sml$ancestry, levels = c("ER", "OB", "C", "D", "AG"))
prop_sml <- prop_sml %>% 
  mutate(sml_group = ifelse(grepl("simulated_", Label), 
                            trimws(gsub("_|-\\d+", " ", Label)), 
                            "real set")) 

prop_sml %>% 
  filter(is.na(code)) %>% 
  filter(proportion > .8) %>% 
  # group_by(sml_group, ancestry) %>% 
  summarise(mean = mean(proportion),
            sd = sd(proportion)) 


## snmf with all real and simulated genotypes
all_geno <- rbind(small_geno, sml_geno)
all_geno_file <- file.path(snmf_dir, "real_sml_chr1_100kSNPs.geno")
write.geno(all_geno, all_geno_file)
all_snmf <- snmf(all_geno_file, K = k,
                 CPU = 18, repetitions = 10,
                 project = "new", entropy = T, seed = 987654321)

best_run <- which.min(cross.entropy(all_snmf, K = k))

# plot ancestry proportion
ind_order <- barchart(all_snmf, K = k, run = best_run, plot = F)
prop_all <- as.data.frame(Q(all_snmf, K = k, run = best_run))
prop_all$Label <- c(seq_info$Label, rownames(sml_geno))
prop_all <- reshape2::melt(prop_all, variable.name = "ancestry", value.name = "proportion")
prop_all <- left_join(prop_all, seq_info, by = "Label")

# prop_all$code <- factor(prop_all$code, levels = unique(prop_all$code)[ind_order$order], ordered = F)
prop_all$Label <- factor(prop_all$Label, levels = unique(prop_all$Label)[ind_order$order], ordered = F)
levels(prop_all$ancestry) <- c("OB", "AG", "D", "C", "ER")
prop_all$ancestry <- factor(prop_all$ancestry, levels = c("ER", "OB", "C", "AG", "D"))
prop_all <- prop_all %>% 
  mutate(sml_group = ifelse(grepl("simulated_", Label), 
                            trimws(gsub("_|-\\d+", " ", Label)), 
                            "real set")) 

plot_dir <- "../../plots/validate_elai"

group_col <- c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF", "#7A0403FF")

tiff(file = file.path(plot_dir, "simulated_100kSNPs_chr1_k5.tiff"),
     # colormodel = "rgb",
     width = 12, height = 4, units = "in", res = 500)
ggplot(prop_all) +
  facet_grid(cols = vars(sml_group), space = "free", scales = "free") + 
  geom_col(aes(x = Label, y = proportion, fill = ancestry), position = position_stack(), width = 1, size = 0) + 
  # scale_fill_viridis_d(end = .9)+
  scale_fill_manual(values = group_col) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Ancestry proportion")  +
  xlab("Individual") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 14),
        # axis.title.x = element_blank(),
        axis.ticks.y = element_line(size = .1),
        axis.text.y = element_text(size = 10),
        axis.line.y = element_line(color = "black", size = .2),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(size = 2, color = "gray60", fill = NA),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(6, "mm"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_blank())
dev.off()


## mean ancestry proportion = 0.9817867, sd = 0.004372262
prop_all %>% 
  filter(is.na(code)) %>% 
  filter(proportion > .8) %>% 
  # group_by(sml_group, ancestry) %>% 
  summarise(mean = mean(proportion),
            sd = sd(proportion)) 

