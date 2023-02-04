#!/usr/bin/env Rscript

### scripts to estimate ancestral proportion of all the African and VN individuals
setwd("/data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/genetic_structure")

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)
library(StAMPP)
library(gt)

# accessions information
seq_info_file <- "sequence_info_final.tsv"
seq_info <- read_tsv(seq_info_file)

# final vcf file containing 1.1M SNPs of all individuals
allvcf_file <- "/data3/projects/vietcaf/baotram/vietcaf_final_biallelic_random_10perc.vcf"

# extract genotypes and convert it to geno format on LEA package
all_vcfR <- read.vcfR(allvcf_file)
all_vcfR@fix[,1] <- gsub(".", "_", all_vcfR@fix[,1], fixed = T)
all_gl <- vcfR2genlight(all_vcfR, n.cores = 20)


## snmf for all inds

all_snmf <- load.snmfProject("./snmf_random_10perc_SNPs/random_10perc_SNPs.snmfProject")


## best k = 5
## choose best run based on cross entropy
k <- 5
best_run <- which.min(cross.entropy(all_snmf, K = k))


## assign groups based on snmf results

prop <- as.data.frame(Q(all_snmf, K = k, run = best_run))
prop$Label <- seq_info$Label
prop <- reshape2::melt(prop, variable.name = "ancestry", value.name = "proportion")
prop <- left_join(prop, seq_info, by = "Label")

levels(prop$ancestry) <- c("C", "D", "OB", "AG", "ER")


prop <- prop %>% 
  group_by(Label) %>% 
  arrange(desc(proportion)) %>%
  slice_head(n = 1) %>% 
  rowwise() %>% 
  mutate(group = ifelse(proportion > 0.8, as.character(ancestry), "hybrid")) %>% 
  mutate(group = ifelse(grepl("VN", as.character(Label)), "Vietnam", group))

all_gl@pop <- as.factor(prop$group)

group_gl <- all_gl[all_gl@pop != "hybrid", sample(all_gl@n.loc, 100000)]
group_fst_100K_SNPs <- stamppFst(group_gl, nclusters = 20)
write_rds(group_fst_100K_SNPs, "./group_fst_100K_SNPs.RDS")


group_gl <- all_gl[all_gl@pop != "hybrid"]
group_fst_1M_SNPs <- stamppFst(group_gl, nclusters = 20)
write_rds(group_fst_1M_SNPs, "./group_fst_1M_SNPs.RDS")

group_fst_1M_SNPs$Fsts[-1,-6] %>% 
  as.data.frame() %>% 
  mutate(across(OB:AG, ~ ifelse(is.na(.x), "", sprintf(.x, fmt = "%.2f")))) %>% 
  rownames_to_column() %>% 
  gt(rowname_col = "rowname", caption = "Pairwise Fst")

table(all_gl@pop)
