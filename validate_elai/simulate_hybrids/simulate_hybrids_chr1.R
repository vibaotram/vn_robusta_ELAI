#!/usr/bin/env Rscript

### simulate hybrids

outdir <- "simulated_hybrids"

## get parental genotypes
parents_vcf_file <- "vietcaf_final_chr01_C013_C034.recode.vcf"
parents_vcfR <-read.vcfR(parents_vcf_file)
parents_gt <- extract.gt(parents_vcfR, return.alleles = T) %>% 
  as.data.frame(stringsAsFactors = F)
parents_gt <- parents_gt %>% 
  tibble::rownames_to_column("rownames") %>% 
  dplyr::mutate(`C013-MERGED` = gsub(".", "?/?", `C013-MERGED`, fixed = T),
                `C034-MERGED` = gsub(".", "?/?", `C034-MERGED`, fixed = T)) %>% 
  tibble::column_to_rownames("rownames")
samp_factor <- 1:nrow(parents_gt)
# samp_factor <- 1:10
parents_gt <- parents_gt[samp_factor,]
rownames(parents_gt) <- gsub("\\.", "_", rownames(parents_gt))


## generate hybrid 3

hybrid <- data.frame("hybrid_3" = parents_gt$`C034-MERGED`,
                     stringsAsFactors = F,
                     row.names = rownames(parents_gt))
hybrid$pos <- gsub(".+_", "", rownames(hybrid)) %>% as.numeric()


### 1 allele from ER + 1 allele from AG
intgr_range_1 <- list(1e7:2e7, 2.003e7:2.795e7, 2.8e7:3e7, 4e7:4.4e7, 4.405e7:5e7)
intgr_1_id <- rownames(hybrid)[hybrid$pos %in% intgr_range_1]
hybrid[intgr_1_id,]$hybrid_3 <- unlist(mclapply(intgr_1_id, function(i) {
  cat(i, "-")
  parent_ER <- parents_gt[i,]$`C034-MERGED`
  parent_AG <- parents_gt[i,]$`C013-MERGED`
  if (any(is.na(c(parent_ER, parent_AG)))) return(NA)
  all1 <- sample(unlist(strsplit(parent_ER, "/")), 1)
  all2 <- sample(unlist(strsplit(parent_AG, "/")), 1)
  hb <- paste(all1, all2, sep = "/")
  return(hb)
}, mc.cores = 20, mc.preschedule = F))

### 2 alleles from AG
intgr_range_2 <- list(2e7:2.003e7, 4.4e7:4.405e7)
intgr_2_id <- rownames(hybrid)[hybrid$pos %in% intgr_range_2]
hybrid[intgr_2_id,]$hybrid_3 <- unlist(mclapply(intgr_2_id, function(i) {
  cat(i, "-")
  hb <- parents_gt[i,]$`C013-MERGED`
  return(hb)
}, mc.cores = 20, mc.preschedule = F))

hybrid$hybrid_3 <- gsub("/", "", hybrid$hybrid_3, fixed = T)
hb_gt <- hybrid %>% select(-pos)
hb3_file <- file.path(outdir, "test_hybrid3")
write_elai_geno(hb_gt, hb3_file)

### true inference of hybrid 3

true_infer <- data.frame(id = rownames(parents_gt), 
                         chr = gsub("\\_\\d{2,8}", "", rownames(parents_gt)), 
                         pos = as.numeric(gsub(".+_", "", rownames(parents_gt))), 
                         group_ER = 1, group_AG = 0)

true_infer <- true_infer %>% mutate(group_ER = case_when(pos %in% unlist(intgr_range_1) ~ 0.5,
                                                         pos %in% unlist(intgr_range_2) ~ 0,
                                                         TRUE ~ 1),
                                    group_AG = case_when(pos %in% unlist(intgr_range_1) ~ 0.5,
                                                         pos %in% unlist(intgr_range_2) ~ 1,
                                                         TRUE ~ 0))
true_inference <- reshape2::melt(true_infer, 
                                 variable.name = "ancestry", 
                                 value.name = "true_dosage", 
                                 measure.vars = c("group_ER", "group_AG"))
fwrite(true_inference, file.path(outdir, "true_inference_hybrid3.tsv"), sep = "\t")




## generate hybrid 4 with all large admixture lengths

hybrid <- data.frame("hybrid_4" = parents_gt$`C034-MERGED`,
                     stringsAsFactors = F,
                     row.names = rownames(parents_gt))
hybrid$pos <- gsub(".+_", "", rownames(hybrid)) %>% as.numeric()

# 1 allele from ER + 1 allele from AG
intgr_range_1 <- c(10e6:12e6, 17e6:19e6, 24e6:26e6, 27e6:30e6, 40e6:41e6, 43e6:44e6, 45e6:46e6, 48e6:50e6)
intgr_1_id <- rownames(hybrid)[hybrid$pos %in% intgr_range_1]
hybrid[intgr_1_id,]$hybrid_4 <- unlist(mclapply(intgr_1_id, function(i) {
  cat(i, "-")
  parent_ER <- parents_gt[i,]$`C034-MERGED`
  parent_AG <- parents_gt[i,]$`C013-MERGED`
  if (any(is.na(c(parent_ER, parent_AG)))) return(NA)
  all1 <- sample(unlist(strsplit(parent_ER, "/")), 1)
  all2 <- sample(unlist(strsplit(parent_AG, "/")), 1)
  hb <- paste(all1, all2, sep = "/")
  return(hb)
}, mc.cores = 20, mc.preschedule = F))

# 2 alleles from AG
intgr_range_2 <- c(12e6:17e6, 26e6:27e6, 46e6:48e6)
intgr_2_id <- rownames(hybrid)[hybrid$pos %in% intgr_range_2]
hybrid[intgr_2_id,]$hybrid_4 <- unlist(mclapply(intgr_2_id, function(i) {
  cat(i, "-")
  hb <- parents_gt[i,]$`C013-MERGED`
  return(hb)
}, mc.cores = 20, mc.preschedule = F))

hybrid$hybrid_4 <- gsub("/", "", hybrid$hybrid_4, fixed = T)
hb_gt <- hybrid %>% select(-pos)
hb4_file <- file.path(outdir, "test_hybrid4")
write_elai_geno(hb_gt, hb4_file)

true_infer <- data.frame(id = rownames(parents_gt), 
                         chr = gsub("\\_\\d{2,8}", "", rownames(parents_gt)), 
                         pos = as.numeric(gsub(".+_", "", rownames(parents_gt))), 
                         group_ER = 1, group_AG = 0)

true_infer <- true_infer %>% mutate(group_ER = case_when(pos %in% unlist(intgr_range_1) ~ 0.5,
                                                         pos %in% unlist(intgr_range_2) ~ 0,
                                                         TRUE ~ 1),
                                    group_AG = case_when(pos %in% unlist(intgr_range_1) ~ 0.5,
                                                         pos %in% unlist(intgr_range_2) ~ 1,
                                                         TRUE ~ 0))
true_inference <- reshape2::melt(true_infer, 
                                 variable.name = "ancestry", 
                                 value.name = "true_dosage", 
                                 measure.vars = c("group_ER", "group_AG"))
fwrite(true_inference, file.path(outdir, "true_inference_hybrid4.tsv"), sep = "\t")




## generate hybrid 5 with all small admixture lengths

hybrid <- data.frame("hybrid_5" = parents_gt$`C034-MERGED`,
                     stringsAsFactors = F,
                     row.names = rownames(parents_gt))
hybrid$pos <- gsub(".+_", "", rownames(hybrid)) %>% as.numeric()

### 1 allele from ER + 1 allele from AG
intgr_range_1 <- list(10e6:12e6, 12.05e6:13.2e6, 13.7e6:15e6, 15.5e6:16e6, 
                      16.05e6:24e6, 24.05e6:30e6, 30.5e6:33e6, 33.05e6:44e6, 
                      44.05e6:46e6, 46.5e6:48e6)
intgr_1_id <- rownames(hybrid)[hybrid$pos %in% intgr_range_1]
hybrid[intgr_1_id,]$hybrid_5 <- unlist(mclapply(intgr_1_id, function(i) {
  cat(i, "-")
  parent_ER <- parents_gt[i,]$`C034-MERGED`
  parent_AG <- parents_gt[i,]$`C013-MERGED`
  if (any(is.na(c(parent_ER, parent_AG)))) return(NA)
  all1 <- sample(unlist(strsplit(parent_ER, "/")), 1)
  all2 <- sample(unlist(strsplit(parent_AG, "/")), 1)
  hb <- paste(all1, all2, sep = "/")
  return(hb)
}, mc.cores = 20, mc.preschedule = F))

### 2 alleles from AG
intgr_range_2 <- list(12e6:12.05e6, 15e6:15.5e6, 24e6:24.05e6, 
                      30e6:30.5e6, 44e6:44.05e6)
intgr_2_id <- rownames(hybrid)[hybrid$pos %in% intgr_range_2]
hybrid[intgr_2_id,]$hybrid_5 <- unlist(mclapply(intgr_2_id, function(i) {
  cat(i, "-")
  hb <- parents_gt[i,]$`C013-MERGED`
  return(hb)
}, mc.cores = 20, mc.preschedule = F))

hybrid$hybrid_5 <- gsub("/", "", hybrid$hybrid_5, fixed = T)
# hybrid$hybrid_5 <- gsub(".", "??", hybrid$hybrid_5, fixed = T)
hb_gt <- hybrid %>% select(-pos)
hb5_file <- file.path(outdir, "test_hybrid5")
write_elai_geno(hb_gt, hb5_file)

true_infer <- data.frame(id = rownames(parents_gt), 
                         chr = gsub("\\_\\d{2,8}", "", rownames(parents_gt)), 
                         pos = as.numeric(gsub(".+_", "", rownames(parents_gt))), 
                         group_ER = 1, group_AG = 0)

true_infer <- true_infer %>% mutate(group_ER = case_when(pos %in% unlist(intgr_range_1) ~ 0.5,
                                                         pos %in% unlist(intgr_range_2) ~ 0,
                                                         TRUE ~ 1),
                                    group_AG = case_when(pos %in% unlist(intgr_range_1) ~ 0.5,
                                                         pos %in% unlist(intgr_range_2) ~ 1,
                                                         TRUE ~ 0))
true_inference <- reshape2::melt(true_infer, 
                                 variable.name = "ancestry", 
                                 value.name = "true_dosage", 
                                 measure.vars = c("group_ER", "group_AG"))
fwrite(true_inference, file.path(outdir, "true_inference_hybrid5.tsv"), sep = "\t")


write_elai_geno <- function(genotypes, file) {
  ## genotypes is a matrix/dataframe, with inds in columns and SNPs in rows, values are in ATGC format
  file.create(file)
  write(ncol(genotypes), file)
  write(nrow(genotypes), file, append = T)
  genotypes <- data.frame("ID" = rownames(genotypes), genotypes)
  data.table::fwrite(genotypes, file, append = T, sep = ",", quote = F, col.names = T, row.names = F)
}