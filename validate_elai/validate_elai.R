#!/usr/bin/env Rscript

library(data.table)
library(Metrics)
library(tidyverse)
library(huxtable)

true_files <- list.files(path = "./simulate_hybrids", pattern = "true",
                         full.names = T, recursive = T)
true_dosage <- do.call(rbind, lapply(true_files, function(f){
  ind <- gsub("(true_inference_|\\.tsv)", "", basename(f))
  df <- fread(f, header = T, sep = "\t", data.table = F)
  grp_df <- df %>% filter(ancestry == "group_ER")
  df <- rbind(df, 
        grp_df %>% mutate(ancestry = "group_OB", 
                    true_dosage = 0),
        grp_df %>% mutate(ancestry = "group_C", 
                          true_dosage = 0),
        grp_df %>% mutate(ancestry = "group_D", 
                          true_dosage = 0))
  df$individual <- ind
  return(df)
}))

ggplot(true_dosage) +
  facet_grid(rows = vars(individual)) +
  geom_line(aes(x = pos, y = true_dosage, color = ancestry))


true_window <- do.call(rbind, lapply(true_files, function(f){
  ind <- gsub("(true_inference_|\\.tsv)", "", basename(f))
  df <- fread(f, header = T, sep = "\t", data.table = F)
  grp_df <- df %>% filter(ancestry == "group_ER", 
                          true_dosage == 0)
  df <- rbind(grp_df, 
              grp_df %>% mutate(ancestry = "group_AG", 
                                true_dosage = 1),
              grp_df %>% mutate(ancestry = "group_OB", 
                                true_dosage = 0),
              grp_df %>% mutate(ancestry = "group_C", 
                                true_dosage = 0),
              grp_df %>% mutate(ancestry = "group_D", 
                                true_dosage = 0))
  df$individual <- ind
  return(df)
}))

elai_results <- list.files(path = "./all_tests", pattern = "local_dosage",
                           full.names = T, recursive = T)

elai_corr <- do.call(rbind, lapply(elai_results, function(f) {
  print(which(elai_results == f))
  param <- basename(dirname(f))
  c <- gsub("(c|_mg.+)", "", param)
  mg <- gsub(".+_mg", "", param)
  snp <- basename(dirname(dirname(f)))
  df <- fread(f, header = T, data.table = F)
  df <- df %>% 
    rename(id = rs) %>% 
    mutate(pos = as.numeric(gsub("CC1_8_Chr01_", "", id)))
  # colnames(df)[1] <- "id"
  # df <- merge(df, true_dosage, 
  #             by = c("id", "pos", "ancestry", "individual"), 
  #             all.x = T)
  df <- left_join(df, true_dosage, by = c("id", "ancestry", "individual", "pos"))
  df %>% 
    # filter(ancestry != "group_D") %>% 
    group_by(individual, ancestry) %>% 
    summarise(correlation = cor(dosage, true_dosage), .groups = "drop") %>% 
    group_by(individual) %>% 
    summarise(mean(correlation**2, na.rm = T))
  eva <- df %>% 
    group_by(individual, pos) %>% 
    summarise(correlation = cor(true_dosage, dosage), .groups = "drop_last") %>% 
    summarise(avg_cor = mean(correlation)) %>% 
    mutate(c = c, mg = mg, snp = snp)
  return(eva)
}))

elai_val <- lapply(elai_results, function(f) {
  print(which(elai_results == f))
  snp <- basename(dirname(dirname(f)))
  param <- basename(dirname(f))
  c <- gsub("c|_mg\\d+", "", param)
  mg <- gsub("c\\d+|_mg", "", param)
  dosage <- fread(f)
  val <- lapply(unique(dosage$individual), function(i) {
    ds <- dosage %>% 
      filter(individual == i) %>% 
      dplyr::rename(id = "rs") %>% 
      dplyr::select(id, ancestry, dosage, individual, pos)
    tds <- true_dosage %>% 
      filter(individual == i, id  %in% ds$id)
    all <- inner_join(tds, ds, by = c("pos", "ancestry", "individual")) %>% 
      # filter(!is.na(window)) %>%
      mutate(true_dosage = ifelse(is.na(true_dosage), 0, true_dosage))
    # cor <- all %>%  
    #   mutate(dif = abs(dosage - true_dosage)) %>% 
    #   group_by(id) %>% 
    #   summarise(acc = 1 - sum(dif), .group = "drop") %>% 
    #   summarise(mean_acc = mean(acc))
    cor <- all %>%
      filter(ancestry == "group_ER") %>%
      summarise(corr = cor(dosage, true_dosage), .groups = "drop") %>%
      summarise(cor = mean(corr, na.rm = T)) %>%
      mutate(snp = snp, c = c, mg = mg, individual = i)
    rmse <- all %>%
      # filter(individual != "hybrid_3") %>%
      mutate(window_length = ifelse(true_dosage == 0.5, NA, window_length)) %>%
      group_by(window_length) %>%
      summarise(rmse = rmse(dosage, true_dosage)) %>%
      mutate(snp = snp, c = c, mg = mg, individual = i)
    return(list(correlation = cor, rmse =rmse))
  })
  cor <- do.call(rbind, lapply(val, function(x) x$correlation))
  # cor %>% 
  #   group_by(id) %>% 
  #   summarise(cor_pos = mean(corr)) %>% 
  #   summarise(cor = mean(cor_pos))
  rmse <- do.call(rbind, lapply(val, function(x) x$rmse))
  return(list(correlation = cor, rmse = rmse))
})
elai_cor <- do.call(rbind, lapply(elai_val, function(x) x$correlation))


plot_dir <- "../plots/validate_elai"


tiff(file = file.path(plot_dir, "all_correlation_dodge.tiff"),
     # colormodel = "rgb",
     width = 7, height = 4, units = "in", res = 1200)
elai_cor %>% 
  group_by(snp, c, mg) %>% 
  summarise(correlation = mean(cor**2),
            sd = sd(cor**2, na.rm = T)) %>%
  mutate(mg = factor(mg, levels = c(5,10,20,30)), 
         c = as.numeric(c),
         snp = factor(snp, levels = c("10k", "100k", "even", "all"))) %>% 
  ggplot(aes(x = mg, y = correlation, col = snp, group = snp)) +
  facet_grid(cols = vars(c), labeller = function(lb) label_both(lb, sep = " = ")) +
  geom_point(size = 2, position = position_dodge(width = 0.6), shape = 18) +
  # geom_line(linewidth = 0.6) + 
  geom_errorbar(aes(ymin = correlation - sd, ymax = correlation + sd), 
                width = 0.6, linewidth = .5, 
                position = position_dodge2(width = 0.6)) + 
  scale_color_viridis_d(name = "SNP set", end = .9, option = "B", direction = -1) +
  # xlim(4, 31) +
  ylim(0.5, 1) +
  ylab(bquote("Squared Pearson's corrrelation (r"^2*")")) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "gray20", color = "gray20"),
        strip.text = element_text(color = "white", size = 12),
        panel.background = element_rect(colour = "gray20"),
        axis.text = element_text(size = 9))
dev.off()

elai_rmse <- do.call(rbind, lapply(elai_results, function(f) {
  print(which(elai_results == f))
  param <- basename(dirname(f))
  c <- gsub("(c|_mg.+)", "", param)
  mg <- gsub(".+_mg", "", param)
  snp <- basename(dirname(dirname(f)))
  df <- fread(f, header = T, data.table = F)
  colnames(df)[1] <- "id"
  df <- merge(df, true_window,
              by = c("id", "pos", "ancestry", "individual"))
  eva <- df %>%
    filter(!is.na(window), ancestry == "group_AG") %>%
    group_by(individual, window_length) %>%
    summarise(rmse_pos = rmse(true_dosage, dosage), .groups = "keep") %>%
    summarise(avg_rmse = mean(rmse_pos)) %>%
    mutate(c = c, mg = mg, snp = snp)
  return(eva)
}))

# elai_rmse <- do.call(rbind, lapply(elai_val, function(x) x$rmse))
tiff(file = file.path(plot_dir, "rmse_c5_rv.tiff"),
     # colormodel = "rgb",
     width = 9, height = 3, units = "in", res = 1200)
elai_rmse %>% 
  filter(window_length %in% c(5e4, 5e5, 1e6, 2e6, 5e6), #) %>% 
         individual != "hybrid_3") %>%
  dplyr::rename(error = "avg_rmse") %>% 
  group_by(window_length, snp, c, mg) %>% 
  summarise(rmse = mean(error),
            sd = sd(error, na.rm = T)) %>% 
  mutate(mg = factor(mg, levels = c(5,10,20,30)), 
         c = factor(c, levels = c(5,15,25)),
         snp = factor(snp, levels = c("10k", "100k", "even", "all"))) %>% 
  filter(c == 5) %>% 
  ggplot(aes(x = mg, y = rmse, fill = snp, color = snp)) +
  facet_grid(cols = vars(paste0(window_length/1e6, " Mb"))) +
  geom_col(position = position_dodge2(width = 0.3, padding = 0.15), width = 0.6) +
  scale_fill_viridis_d(name = "SNP set", end = .9, option = "B", direction = -1) +
  scale_color_viridis_d(name = "SNP set", end = .9, option = "B", direction = -1) +
  ylim(0, 0.52) +
  ylab("RMSE") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "gray20", color = "gray20"),
        strip.text = element_text(color = "white", size = 12),
        panel.background = element_rect(colour = "gray20"),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank()
        ) 
dev.off()


tiff(file = file.path(plot_dir, "rmse_allc_rv.tiff"),
     # colormodel = "rgb",
     width = 10, height = 9, units = "in", res = 1200)
elai_rmse %>% 
  filter(window_length %in% c(5e4, 5e5, 1e6, 2e6, 5e6), #) %>% 
         individual != "hybrid_3") %>%
  dplyr::rename(error = "avg_rmse") %>% 
  group_by(window_length, snp, c, mg) %>% 
  summarise(rmse = mean(error),
            sd = sd(error, na.rm = T)) %>% 
  mutate(mg = factor(mg, levels = c(5,10,20,30)), 
         c = factor(paste0("c = ", c), levels = c("c = 5","c = 15","c = 25")),
         snp = factor(snp, levels = c("10k", "100k", "even", "all"))) %>% 
  # filter(c == 5) %>% 
  ggplot(aes(x = mg, y = rmse, fill = snp, color = snp)) +
  facet_grid(cols = vars(paste0(window_length/1e6, " Mb")),
             rows = vars(c)) +
  geom_col(position = position_dodge2(width = 0.3, padding = 0.15), width = 0.6) +
  scale_fill_viridis_d(name = "SNP set", end = .9, option = "B", direction = -1) +
  scale_color_viridis_d(name = "SNP set", end = .9, option = "B", direction = -1) +
  ylim(0, 0.52) +
  ylab("RMSE") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "gray20", color = "gray20"),
        strip.text = element_text(color = "white", size = 12),
        panel.background = element_rect(colour = "gray20"),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank()) 
dev.off()



## benchmark
bm_files <- list.files(".", pattern = "benchmark.txt", 
                       full.name = T, recursive =  T)
timecount <- function(x) {
  d <- ifelse(str_detect(x, "-"), as.numeric(gsub("-.+", "", x)), 0)
  hms <- as.numeric(str_split(gsub(".+-", "", x), ":")[[1]])
  nb_hrs <- d*24+hms[1]+hms[2]/60+hms[3]/3600
  nb_hrs <- sprintf("%.1f", nb_hrs)
  return(nb_hrs)
}

bm <- do.call(rbind, lapply(bm_files, function(f) {
  print(f)
  df <- read.table(f, skip = 3, col.names = c("script", "time", "memory"))
  snp <- basename(dirname(dirname(dirname(f))))
  param <- basename(dirname(f))
  c <- gsub("(c|_mg.+)", "", param)
  mg <- gsub(".+_mg", "", param)
  df <- df %>% 
    mutate(c = c, mg = mg, snp = snp) %>% 
    mutate(memory = sprintf("%.1f", as.numeric(gsub("K", "", memory))/1e6),
           time = timecount(time)) %>% 
    select(-script) %>% 
    relocate(snp, c, mg, time, memory)
  return(df)
  }))

bmw <- bm %>% 
  mutate(mg = as.numeric(mg),
         c = as.numeric(c),
         snp = fct_relevel(snp, c("10k", "100k", "even"))) %>% 
  arrange(snp, c, mg) %>% 
  pivot_wider(names_from = c, values_from = c(time, memory))

bmw %>% 
  hux() %>% 
  set_contents(1, 3:8, gsub(".+_", "", colnames(bmw)[3:8])) %>% 
  insert_row("", "", "Time (hours)", "Time (hours)", "Time (hours)", 
             "Memory (Gb)", "Memory (Gb)", "Memory (Gb)", after = 0) %>%
  insert_row("", "", "c", "c", "c", "c", "c", "c", after = 1) %>% 
  merge_cells(2, 3:5) %>% 
  merge_cells(2, 6:8) %>% 
  merge_cells(4:7, 1) %>% 
  merge_cells(8:11,1) %>% 
  merge_cells(12:15,1) %>% 
  merge_cells(1, 3:5) %>%
  merge_cells(1, 6:8) %>%
  set_right_border(everywhere, everywhere, brdr(1, "solid", "grey")) %>% 
  set_bottom_border(everywhere, everywhere, brdr(1, "solid", "grey")) %>% 
  set_bold(row = 1:3, col = everywhere) %>%
  set_bold(row = everywhere, col = 1:2) %>% 
  set_align(everywhere, everywhere, "center") %>% 
  set_valign(everywhere, everywhere, "middle") %>% 
  quick_docx(file = "time_mem.docx")

