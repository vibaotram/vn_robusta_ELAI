#!/usr/bin/env Rscript

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)
library(viridis)
library(aplot)

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

## run snmf with all SNPs


snmf_dir <- "snmf_chr1"
chr1_geno_file <- file.path(snmf_dir, "chr1.geno")
write.geno(chr1_geno, chr1_geno_file)
chr1_snmf <- snmf(chr1_geno_file, K = k,
                  CPU = 18, repetitions = 10,
                  project = "new", entropy = T, seed = 987654321)

best_run <- which.min(cross.entropy(chr1_snmf, K = k))


ind_order <- barchart(chr1_snmf, K = k, run = best_run, plot = F)
prop_real <- as.data.frame(Q(chr1_snmf, K = k, run = best_run))
prop_real$Label <- seq_info$Label
prop_real <- reshape2::melt(prop_real, variable.name = "ancestry", value.name = "proportion")
prop_real <- left_join(prop_real, seq_info, by = "Label")

prop_real$code <- factor(prop_real$code, levels = unique(prop_real$code)[ind_order$order], ordered = F)
levels(prop_real$ancestry) <- c("C", "AG", "OB", "D", "ER")
prop_real$ancestry <- factor(prop_real$ancestry, levels = c("ER", "OB", "C", "D", "AG"))



plot_dir <- "../../plots/validate_elai"

tiff(file = file.path(plot_dir, "population_structure_chr1_k5.tiff"),
     # colormodel = "rgb",
     width = 12, height = 5, units = "in", res = 1200)
ggplot(prop_real) +
  geom_col(aes(x = code, y = proportion, fill = ancestry), position = position_stack(), width = .95, size = 0) + 
  scale_fill_viridis_d(end = .9)+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Ancestry proportion")  +
  xlab("Individual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(5, "mm"),
        axis.ticks.y = element_line(size = .1),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(1, "mm"),
        strip.text.x = element_text(angle = 90, size = 10, hjust = 1)) 
dev.off()


# get ancestral geno freq
chr1_freq <- G(chr1_snmf, K = 5, run = best_run)
colnames(chr1_freq) <- levels(prop_real$ancestry)







chr1_snps <- as.data.frame(getFIX(chr1_vcfR))
chr1_snps$ID <- paste(chr1_snps$CHROM, chr1_snps$POS, sep = "_")

n_inds <- 100
n_snps <- nrow(chr1_freq)/3


outdir <- "simulated_source"
mclapply(1:ncol(chr1_freq), function(v) {
  grp_elai_input <- file.path(outdir, paste0("group_", colnames(chr1_freq)[v]))
  ind_names <- paste0("group_", colnames(chr1_freq)[v], "-", 1:n_inds)
  writeLines(c(n_inds, n_snps, paste(c("ID", ind_names), collapse = ",")), grp_elai_input)
  freq <- chr1_freq[,v] # get allele freq of the group
  geno <- matrix(nrow=n_inds, ncol=n_snps) # genotype matrix of the group
  for (i in 1:n_snps) { 
    # for locus i, frequencies of 0, 1, 2 genotype are in row 3i-2, 3i-1, 3i of freq matrix
    sml_geno <-  sample(c(0,1,2), n_inds, prob = c(freq[3*i-2], freq[3*i-1], freq[3*i]), replace = T)
    snp_id <- chr1_gl@loc.names[i]
    sml_gt <- geno2elai_gt(sml_geno, chr1_snps, snp_id)
    geno_li <- paste(c(snp_id, sml_gt), collapse = ",")
    cat(paste0(geno_li, "\n"), file = grp_elai_input, append = T)
  }
  return(v)
}, mc.cores = 6, mc.preschedule = F)