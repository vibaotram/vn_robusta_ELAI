#!/usr/bin/env Rscript

### chromosome painting of the VN individuals
### using the results from the workflow

library(LEA)
library(vcfR)
library(tidyverse)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(karyoploteR)
library(Rsamtools)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(bstfun)
library(gt)
library(gridExtra)

group_col <- c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF", "#7A0403FF")
names(group_col) <- c("group_ER", "group_OB", "group_C", "group_AG", "group_D")

## functions to plot raw ELAI results
plot_elai <- function(snakedir, chrom, groups, plot_path = NULL) {
  snp <- c("all", "even")
  all_dosage <- do.call(rbind, lapply(snp, function(s) {
    result_data <- file.path(snakedir, "elai_results", chrom, s, "c5_mg20", "local_dosage.txt")
    print(result_data)
    all_dosage <- fread(result_data, sep = "\t", header = T, data.table = T)
    all_dosage1 <- data.frame(all_dosage, snp = s)
    return(all_dosage1)
  }))
  
  all_dosage$ancestry <- factor(all_dosage$ancestry)
  levels(all_dosage$ancestry) <- groups
  all_dosage$ancestry <- factor(all_dosage$ancestry, levels = c("group_ER", "group_OB", "group_C", "group_AG", "group_D"))
  
  if (!is.null(plot_path)) {
    cat(plot_path)
    # tiff(file = plot_path, 
         # width = 20, height = 12, units = "in", res = 800)
    p <- ggplot(all_dosage) +
      facet_grid(rows = vars(individual)) +
      # geom_point(aes(x = pos, y = dosage, color = ancestry), size = 1.5, shape = 16, alpha = .5) +
      geom_line(aes(x = pos, y = dosage, color = ancestry, linetype = snp, alpha = snp)) +
      scale_linetype_manual(values=c("dashed", "solid")) +
      # scale_color_viridis_d(end = 0.9) +
      scale_color_manual(values = group_col) +
      scale_alpha_manual(values = c(.8, 1)) +
      xlab(paste("Position on chromosome", unique(all_dosage$chr))) +
      ylab("Ancestral dosage") +
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
    ggsave(filename = plot_path, plot = p, width = 10, height = 6, dpi = 1000)
    # dev.off()
  }
  
  return(invisible(all_dosage))
}

## function to merge the raw ELAI results and convert to GR object, for each chromosome
final_elai <- function(snakedir, chrom, groups, plot_path) {
  dose <- plot_elai(snakedir, chrom, groups, plot_path)
  # gr <- elai_GR(dose)
  # return(gr)
}

plot_dir <- "../plots/elai_tr"
dir.create(plot_dir, showWarnings = F)

snakedir <- "./elai"

## in each chromosome, use snmf results to rename groups, then plot the raw ELAI results and create GR of the final dosage
chr <- "CC1.8.Chr01"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_C", "group_AG", "group_OB", "group_D", "group_ER")
dose_chr1 <- plot_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr1.tiff"))
elai_GR1 <- elai_GR(dose_chr1)

chr <- "CC1.8.Chr02"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_ER", "group_AG", "group_OB", "group_C", "group_D")
elai_GR2 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr2.tiff"))

chr <- "CC1.8.Chr03"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_D", "group_OB", "group_ER", "group_AG", "group_C")
elai_GR3 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr3.tiff"))

chr <- "CC1.8.Chr04"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_D", "group_ER", "group_C", "group_AG", "group_OB")
elai_GR4 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr4.tiff"))

chr <- "CC1.8.Chr05"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_D", "group_OB", "group_ER", "group_C", "group_AG")
elai_GR5 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr5.tiff"))

chr <- "CC1.8.Chr06"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_C", "group_ER", "group_D", "group_OB", "group_AG")
elai_GR6 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr6.tiff"))

chr <- "CC1.8.Chr07"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_AG", "group_D", "group_C", "group_OB", "group_ER")
elai_GR7 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr7.tiff"))

chr <- "CC1.8.Chr08"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_AG", "group_D", "group_OB", "group_ER", "group_C")
elai_GR8 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr8.tiff"))

chr <- "CC1.8.Chr09"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_ER", "group_C", "group_D", "group_OB", "group_AG")
elai_GR9 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr9.tiff"))

chr <- "CC1.8.Chr10"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_OB", "group_D", "group_C", "group_AG", "group_ER")
elai_GR10 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr10.tiff"))

chr <- "CC1.8.Chr11"
get_snmf(snakedir = snakedir, chrom = chr)
group_id <- c("group_OB", "group_ER", "group_AG", "group_D", "group_C")
elai_GR11 <- final_elai(snakedir, chr, group_id, plot_path = file.path(plot_dir, "raw_chr11.tiff"))

## combine all the GR
all_windows <- GRangesList(elai_GR1, elai_GR2, elai_GR3, elai_GR4, elai_GR5, elai_GR6, elai_GR7,
                           elai_GR8, elai_GR9, elai_GR10, elai_GR11) %>% unlist
seqlevels(all_windows) <- gsub("CC1.8.Chr", "Chr ", seqlevels(all_windows))

## plot
# genome
genome_fa <- "/data3/projects/vietcaf/reference/CC1.8_v2_pseudomolecule_cat.fa"
genome <- scanFa(genome_fa)
genome <- genome[grepl("Chr", names(genome))] # remove contigs
names(genome) <- gsub(" .+", "", names(genome))
names(genome) <- gsub("CC1.8.Chr", "Chr ", names(genome))
genome_GR <- GRanges(seqnames = names(genome), ranges = IRanges(start = 1, end = width(genome)))


plot_dir <- "../plots/elai_tr"



pp1 <- getDefaultPlotParams(plot.type=6)
pp1$data1height <- 5
pp1$data2height <- 5
pp1$topmargin <- 20
pp1$bottommargin <- 20
pp1$leftmargin <- .07
pp1$rightmargin <- .15

## plot TR6
id <- "TR6m"
all_windows_ind <- all_windows[all_windows$individual == id]
all_windows_large <- all_windows_ind[width(all_windows_ind) > 1e6,]
# all_windows_small <- all_windows_ind[width(all_windows_ind) < 1e6,]
groups <- colnames(mcols(all_windows))[colnames(mcols(all_windows)) != "individual"]

tiff(file.path(plot_dir, paste0(id, ".tiff")),
     width = 11, heigh = 6, unit = "in", res = 1000)
par(mar = c(2,2,2,2))
gkp <- plotKaryotype(genome = genome_GR, plot.type = 6,
                     plot.params = pp1) 
kpDataBackground(gkp, data.panel ="ideogram", color = "grey90")
kpAddBaseNumbers(gkp, tick.dist = 1e7)
for (g in 1:length(groups)) {
  if (g == 1) {
    ga <- 0
    gz <- mcols(all_windows_large)[,groups[g]]
  } else {
    ga <- gz
    gz <- gz+mcols(all_windows_large)[,groups[g]]
  }
  kpBars(gkp, data = all_windows_large, y0 = ga, y1 = gz, col = group_col[groups[g]], border = NA)
}
kpAxis(gkp, side = 2, data.panel=1, cex = 0.7)
legend("right", #x = 0.9, y = 0.4, #"bottomright", #  yjust = .5, xjust = 0, # 
       title = "ancestry", legend = c(gsub("group_", "", names(group_col)), "undetermined"), 
       fill = c(group_col, "grey90"), bty = "n", xjust = 0)
dev.off()

# plot all genomes together
pp2 <- getDefaultPlotParams(plot.type=7)
pp2$data1outmargin <- 0
# pp2$data1height <- 5
# pp2$data2height <- 5
pp2$topmargin <- 40
# pp2$bottommargin <- 20
pp2$leftmargin <- .07
pp2$rightmargin <- .2
pp2$ideogramlateralmargin <- 0.003

inds <- unique(all_windows$individual) %>% as.character() %>% str_sort(numeric = T, decreasing = T)

tiff(file.path(plot_dir, "all_TRs.tiff"),
     width = 20, heigh = 8, unit = "in", res = 1000)
par(mar = c(2,2,2,2))
gkp <- plotKaryotype(genome = genome_GR, plot.type = 7,
                     plot.params = pp2)
kpDataBackground(gkp, data.panel ="ideogram", color = "white")
kpAddBaseNumbers(gkp, tick.dist = 1e7)
for (id in 1:length(inds)) {
  if (id == 1) {
    ra <- 0
    rz <- 1/length(inds)-.02
  } else {
    ra <- rz +.02
    rz <- rz + 1/length(inds)
  }
  all_windows_ind <- all_windows[all_windows$individual == inds[id]]
  all_windows_large <- all_windows_ind[width(all_windows_ind) > 1e6,]
  all_windows_small <- all_windows_ind[width(all_windows_ind) < 1e6,]
  groups <- colnames(mcols(all_windows))[colnames(mcols(all_windows)) != "individual"]
  
  kpBars(gkp, data = genome_GR, y0 = rep(0, length(genome_GR)), y1 = rep(1, length(genome_GR)), col = "gray80", border = NA, r0 = ra, r1 = rz)
  kpAxis(gkp, side = 1, data.panel=1, cex = 0.7, r0 = ra, r1 = rz)
  kpText(gkp, chr=seqnames(genome_GR)[length(genome_GR)], x=width(genome_GR)[length(genome_GR)] + 2.6e7, y=0.5, labels = gsub("m", "", inds[id]), r0 = ra, r1 = rz)
  for (g in 1:length(groups)) {
    if (g == 1) {
      ga <- 0
      gz <- mcols(all_windows_large)[,groups[g]]
    } else {
      ga <- gz
      gz <- gz+mcols(all_windows_large)[,groups[g]]
    }
    kpBars(gkp, data = all_windows_large, y0 = ga, y1 = gz, col = group_col[groups[g]], border = NA, r0 = ra, r1 = rz)
  }
  # kpPlotRegions(gkp, data = all_windows_large, data.panel = "ideogram", col = "grey20", avoid.overlapping = T, r0 = ra, r1 = rz)
  # kpPlotRegions(gkp, data = all_windows_small, data.panel = "ideogram", col = "orange3", avoid.overlapping = T)
}
legend("right", #x = 0.9, y = 0.4, #"bottomright", #  yjust = .5, xjust = 0, # 
       title = "ancestry", legend = c(gsub("group_", "", names(group_col)), "undetermined"), 
       fill = c(group_col, "grey90"), bty = "n", xjust = 0, cex = 1.2)
dev.off()

## summary table for overall proportion

sum_TRs <- do.call(rbind, lapply(unique(all_windows$individual), function(i) {
  # cat("> individual:", i, "\n")
  all_windows_ind <- all_windows[all_windows$individual == i,]
  all_windows_ind <- all_windows_ind[width(all_windows_ind) > 1e6]
  grps <- mcols(all_windows_ind) %>% as.data.frame() %>% dplyr::select(contains("group")) %>% colnames()
  glob_adm <- do.call(rbind, lapply(grps, function(g) {
    pw <- sum(mcols(all_windows_ind)[,g]*width(all_windows_ind))/sum(width(genome_GR))
    return(pw)
  }))
  glob_adm <- as.data.frame(t(glob_adm))
  colnames(glob_adm) <- grps
  rownames(glob_adm) <- i
  return(glob_adm)
}))

sum_TRs <- sum_TRs[stringr::str_sort(rownames(sum_TRs), numeric = T), stringr::str_sort(colnames(sum_TRs))]
undetermined <- 1 - rowSums(sum_TRs)
sum_TRs$undetermined <- undetermined
colnames(sum_TRs) <- gsub("group_", "ancestry ", colnames(sum_TRs))
rownames(sum_TRs) <- gsub("m", "", rownames(sum_TRs))

t <- sum_TRs %>% 
  as.data.frame() %>% 
  gt(rownames_to_stub = T) %>% 
  fmt_number(columns = everything(), decimals = 3) %>% 
  data_color(colors = scales::col_numeric(palette = c("white", "#FB8022FF"), domain = c(0,1)),
             columns = "ancestry AG") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", "#30123BFF"), domain = c(0,1)),
             columns = "ancestry ER") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", "#28BBECFF"), domain = c(0,1)),
             columns = "ancestry OB") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", "black"), domain = c(0,1)),
             columns = "undetermined") 

gt <- as_ggplot(t, vwidth = 2400, vheight = 2000)

## mapping of groups
group_col <- c("#FB8022FF", "#A2FC3CFF", "#7A0403FF", "#30123BFF", "#28BBECFF")

africa <- ne_countries(scale = "medium", continent = "africa", returnclass = "sf")


africa <- africa %>% 
  mutate(ancestry = case_when(sovereignt %in% c("Ivory Coast", "Guinea", "Ghana") ~ "D",
                              sovereignt == "Cameroon" ~ "C",
                              sovereignt %in% c("Gabon", "Benin", "Angola") ~ "AG",
                              sovereignt %in% c("Uganda", "Central African Republic") ~ "OB",
                              sovereignt == "Democratic Republic of the Congo" ~ "ER"))

rob_countries <- africa %>% filter(!is.na(ancestry))

m <- ggplot(data = africa) +
  geom_sf(fill = "white", color = "grey") +
  geom_sf(data = rob_countries, aes(fill = ancestry), color = "grey") +
  geom_sf_text(data = rob_countries, aes(label = abbrev), color = "grey50") +
  scale_fill_manual(values = group_col, na.value = "grey50") +
  coord_sf(xlim = c(-15, 48), ylim = c(-40, 35), expand = TRUE) +
  theme_void() +
  theme(legend.box.margin = margin(l = 20), 
        legend.key.size = unit(14, "pt"),
        legend.key.height = unit(15, "pt"))

## plot the summary table with the map
plot_dir <- "/home/baotram/phd/robusta/writing/lai_paper/plots/elai_tr"
tiff(file.path(plot_dir, "table.tiff"),
     width = 8, heigh = 10, unit = "in", res = 1000)
gridExtra::grid.arrange(m, gt, heights = c( 1.2, 1))
dev.off()



## mapping of native population

ref_inds <- read_tsv("../genetic_structure/sequence_info_final.tsv")
ref_inds %<>% 
  filter(SSR_Group != "VN" | is.na(SSR_Group)) %>% 
  rowwise() %>% 
  mutate(ancestry = case_when(SSR_Group %in% c("O", "B") ~ "OB",
                              SSR_Group %in% c("E", "R") ~ "ER",
                              SSR_Group %in% c("A", "G") ~ "AG",
                              SSR_Group == "EB" ~ "ER",
                              TRUE ~ SSR_Group)) %>% 
  filter(!is.na(ancestry), 
         !(ancestry == "ER" && Country %in% c("Cameroon", "Guinea", "Gabon", "RCA")))

ggplot(data = africa) +
  geom_sf(fill = "white", color = "grey") +
  # geom_sf(data = rob_countries, aes(fill = ancestry), color = "grey") +
  # geom_point(data = ref_inds, aes(x = Long, y = Lat, color = ancestry), size = 1) +
  stat_ellipse(data = ref_inds, 
               aes(x = Long, y = Lat, group = ancestry, fill = ancestry), 
               type = "norm", level = 0.6, geom = "polygon") +
  geom_sf_text(data = rob_countries, aes(label = abbrev), color = "gray50") +
  scale_fill_manual(values = group_col, na.value = "grey50") +
  scale_color_manual(values = group_col, na.value = "grey50") +
  coord_sf(xlim = c(-15, 48), ylim = c(-40, 35), expand = TRUE) +
  theme_void() +
  theme(legend.box.margin = margin(l = 20), 
        legend.key.size = unit(14, "pt"),
        legend.key.height = unit(15, "pt"))
