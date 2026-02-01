

########################################################
#############Fig_5
########################################################


## =========================================================
##  Fig_5_A: Codon frequency comparison (Others vs Chimera)
## =========================================================
rm(list = ls(all = TRUE))
## 0) Packages ------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(ggrepel)
  library(ggplot2)
})

## 1) Load data ----------------------------------------------
data_path <- "data/codon_data.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: codon_data


perform_chi_square <- function(f1, f2, total_count = 64) {
  observed <- c(f1 * total_count, f2 * total_count)
  expected <- rep(mean(observed), 2)
  chisq.test(observed, p = expected / sum(expected))$p.value
}

codon_data$p_value <- mapply(
  perform_chi_square,
  codon_data$Others ,
  codon_data$Chimera
)
codon_data$p_value=sprintf("%.2e", codon_data$p_value)
codon_data=codon_data[order(codon_data$Others),]
codon_data_long <- codon_data %>%
  select(Codon, Others, Chimera) %>%
  pivot_longer(cols = c(Others, Chimera),
               names_to = "Type",
               values_to = "Frequency")

codon_data$p_value <- as.numeric(codon_data$p_value)

codon_data_long$Codon  <- factor(codon_data_long$Codon,levels =codon_data$Codon)


p_points <- ggscatter(
  codon_data_long,
  x = "Codon",
  y = "Frequency",
  color = "Type",
  palette = c('#377EB8', '#FF7F00'),
  size = 2,
  alpha = 0.8
) +
  stat_summary(
    aes(group = Type, color = Type),
    fun = mean,
    geom = "point",
    size = 2.5,
    show.legend = FALSE
  ) +
  # mean ± SE
  stat_summary(
    aes(group = Type, color = Type),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.25,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  labs(
    x = "Codon",
    y = "Frequency (%)",
    color = "Group"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
  ) +
  geom_point(
    data = subset(codon_data, p_value < 0.05),
    aes(x = Codon, y = Chimera),
    color = "#E41A1C",
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = subset(codon_data, p_value < 0.05),
    aes(x = Codon, y = Chimera, label = AA),
    size = 3.5,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.4,
    max.overlaps = Inf,
    inherit.aes = FALSE
  )

print(p_points)

pdf(file='Fig_5_A.pdf',width =11,height = 5.5 )
print(p_points)
dev.off()



## =========================================================
##  Fig 5B: Abundance of vOTUs with chimera across environments
## =========================================================

## 0) Session setup ------------------------------------------
rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggpubr)
  library(ggplot2)
})

## 1) Load data ----------------------------------------------
data_path <- "data/all_type_Abundance.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: all_type_Abundance

if (!exists("all_type_Abundance")) {
  stop("Object `all_type_Abundance` was not found after loading: ", data_path)
}
stopifnot(is.data.frame(all_type_Abundance))

required_cols <- c("type", "rowsum")
missing_cols <- setdiff(required_cols, colnames(all_type_Abundance))
if (length(missing_cols) > 0) {
  stop("Missing required columns in `all_type_Abundance`: ", paste(missing_cols, collapse = ", "))
}

## 2) Order environments by median ----------------------------
type_order <- all_type_Abundance %>%
  group_by(type) %>%
  summarise(median_rowsum = median(rowsum, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_rowsum) %>%
  pull(type)

all_type_Abundance <- all_type_Abundance %>%
  mutate(type = factor(type, levels = type_order)) %>%
  filter(!is.na(type), type != "")

## 3) Pairwise Wilcoxon tests (optional) ----------------------
# Warning: many groups => many comparisons => cluttered annotations.
run_pairwise_tests <- TRUE
p_adjust_method <- "BH"   # more standard than raw p<0.05
alpha_cutoff <- 0.05

significant_pairs <- list()

if (run_pairwise_tests) {
  pairwise_df <- combn(levels(all_type_Abundance$type), 2, simplify = FALSE) %>%
    map_df(\(pair) {
      g1 <- pair[1]; g2 <- pair[2]
      sub_df <- all_type_Abundance %>% filter(type %in% c(g1, g2))
      wt <- wilcox.test(rowsum ~ type, data = sub_df)
      tibble(group1 = g1, group2 = g2, p_value = wt$p.value)
    }) %>%
    mutate(q_value = p.adjust(p_value, method = p_adjust_method))
  
  significant_pairs <- pairwise_df %>%
    filter(q_value < alpha_cutoff) %>%
    mutate(comparison = map2(group1, group2, ~ c(.x, .y))) %>%
    pull(comparison)
}

## 4) Color palette ------------------------------------------
type_colors <- c(
  "transportation or  hubs" = "#F47F60",
  "wwtp"                    = "#42A5F5",
  "air"                     = "#FFB74D",
  "biofilm"                 = "#6D7EC9",
  "bioreactor"              = "#D3D656",
  "coast"                   = "#9575CD",
  "Fresh water"             = "#26A69A",
  "hospital"                = "#DD6895",
  "indoor"                  = "#F48FB1",
  "indoor air"              = "#4CAF50",
  "outdoor"                 = "#26C6DA",
  "sludge"                  = "#7E57C2",
  "soil"                    = "#C5E1A5"
)

# Add fallback colors if there are levels not in the palette
missing_colors <- setdiff(levels(all_type_Abundance$type), names(type_colors))
if (length(missing_colors) > 0) {
  fallback <- rep("#BDBDBD", length(missing_colors))
  names(fallback) <- missing_colors
  type_colors <- c(type_colors, fallback)
}

## 5) Plot ----------------------------------------------------
p_box_type <- ggboxplot(
  all_type_Abundance,
  x = "type",
  y = "rowsum",
  fill = "type",
  bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.5,
  size = 0.5,
  outlier.shape = NA,
  short.panel.labs = TRUE
) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = NULL,
    y = "log2(number of vOTUs with Chimera)"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# OPTIONAL: add significance annotations (can be very crowded)
# If you really want to show only significant comparisons:
# if (length(significant_pairs) > 0) {
#   p_box_type <- p_box_type +
#     stat_compare_means(
#       comparisons = significant_pairs,
#       hide.ns = TRUE,
#       label = "p.signif"
#     )
# }

print(p_box_type)

## 6) Save ----------------------------------------------------
ggsave(
  filename = "Fig_5_B.pdf",
  plot = p_box_type,
  width = 12,
  height = 8,
  units = "in"
)

## =========================================================
##  Fig 5D: ARG category percentages (pie chart)
## =========================================================

## 0) Session setup ------------------------------------------
rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggpubr)   # ggpie()
  library(ggplot2)
})

## 1) Load data ----------------------------------------------
data_path <- "data/ARG_types.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: pieARG

if (!exists("pieARG")) stop("Object `pieARG` was not found after loading: ", data_path)
stopifnot(is.data.frame(pieARG))

required_cols <- c("ARG", "count")
missing_cols <- setdiff(required_cols, colnames(pieARG))
if (length(missing_cols) > 0) {
  stop("Missing required columns in `pieARG`: ", paste(missing_cols, collapse = ", "))
}

## 2) Color palette ------------------------------------------
# Keep your palette, but ensure it is long enough for all ARG categories.
base_colors <- c(
  "#90a955", "#00A6FB", "#0582CA", "#006494",
  "#c9e4ca", "#31572c", "#c9184a", "#ff4d6d",
  "#FF5733", "#66FFFF", "#9933CC", "#FF6699", "#339966", "#996633",
  "#FFFF33", "#663333", "#FF3366", "#99CCFF", "#33CC99", "#FF99CC",
  "#333399", "#FF9933", "#66CC66", "#a4133c",
  "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2"
)

n_cat <- n_distinct(pieARG$ARG)

# If palette is shorter than categories, extend it safely
if (length(base_colors) < n_cat) {
  palette_colors <- grDevices::colorRampPalette(base_colors)(n_cat)
} else {
  palette_colors <- base_colors[seq_len(n_cat)]
}

## 3) Compute percentages & labels ----------------------------
pieARG2 <- pieARG %>%
  mutate(
    percentage = count / sum(count, na.rm = TRUE) * 100,
    label = paste0(ARG, " (", round(percentage, 1), "%)")
  ) %>%
  arrange(desc(count)) %>%
  mutate(ARG = factor(ARG, levels = ARG))

# Show labels only for top N categories
top_n_labels <- 10
pieARG2 <- pieARG2 %>%
  mutate(show_label = if_else(row_number() <= top_n_labels, label, ""))

## 4) Plot ----------------------------------------------------
p_fig5d <- ggpie(
  data = pieARG2,
  x = "count",
  fill = "ARG",
  label = "show_label",
  color = "white",
  lab.pos = "out",
  lab.adjust = 1.2,
  lab.font = c(4, "plain", "black"),
  palette = palette_colors,
  legend = "right"
) +
  labs(title = "ARG category composition") +
  theme(plot.title = element_text(hjust = 0.5))

print(p_fig5d)

## 5) Save ----------------------------------------------------
ggsave(
  filename = "Fig_5_D.pdf",
  plot = p_fig5d,
  width = 10,
  height = 6,
  units = "in"
)



## =========================================================
##  Fig 5C: ARG host composition (Phylum + Family) circos rings
## =========================================================
## 0) Session setup ------------------------------------------
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(circlize)
  library(RColorBrewer)
})

## 1) Load data ----------------------------------------------
data_path <- "data/ARG_host.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected objects: df (at least)

if (!exists("df")) stop("Object `df` was not found after loading: ", data_path)
stopifnot(is.data.frame(df))

required_cols <- c("votu", "host_phylum", "host_family")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in `df`: ", paste(missing_cols, collapse = ", "))
}

## 2) Helper: keep top-N in a column, others -> "Other" -------
keep_top_n_levels <- function(data, col, n = 10) {
  stopifnot(col %in% names(data))
  data2 <- data %>%
    filter(!is.na(.data[[col]]), .data[[col]] != "")
  
  top_levels <- data2 %>%
    count(.data[[col]], sort = TRUE) %>%
    slice_head(n = n) %>%
    pull(1)
  
  data %>%
    mutate("{col}" := if_else(.data[[col]] %in% top_levels, .data[[col]], "Other"))
}

## 3) Prepare phylum/family table -----------------------------
# Keep only top phyla; others -> "Other"
top_phylum_n <- 10
df2 <- keep_top_n_levels(df, "host_phylum", n = top_phylum_n)

pf_df <- df2 %>%
  transmute(
    votu = votu,
    phylum = host_phylum,
    family = host_family
  ) %>%
  distinct() %>%
  select(-votu) %>%
  mutate(
    phylum = if_else(is.na(phylum) | phylum == "", "Other", phylum),
    family = if_else(is.na(family) | family == "", "Other", family)
  )

## 4) Rank phyla and collapse families within phyla -----------
# Phylum frequency (Other forced last later)
phylum_freq <- pf_df %>%
  count(phylum, name = "freq") %>%
  arrange(desc(freq))

# Define "top 1" and "rank 2" groups (your original intent)
top1_phylum <- phylum_freq$phylum[1]
top2to5_phyla <- phylum_freq$phylum[2:5]

pf_df2 <- pf_df %>%
  mutate(
    phylum_rank = case_when(
      phylum == top1_phylum ~ 1L,
      phylum %in% top2to5_phyla ~ 2L,
      TRUE ~ 3L
    )
  )

# Family counts within each phylum
family_counts <- pf_df2 %>%
  count(phylum, family, name = "family_freq") %>%
  group_by(phylum) %>%
  arrange(desc(family_freq), .by_group = TRUE) %>%
  mutate(family_rank = row_number()) %>%
  ungroup()

# Collapse families depending on phylum rank:
# - top1 phylum keeps top 5 families
# - phylum rank 2 (2-5) keeps top 3 families
pf_df3 <- pf_df2 %>%
  left_join(family_counts, by = c("phylum", "family")) %>%
  mutate(
    family = case_when(
      phylum_rank == 1L & family_rank > 5L ~ "Other",
      phylum_rank == 2L & family_rank > 3L ~ "Other",
      phylum == "Other" ~ "Other",
      TRUE ~ family
    )
  ) %>%
  select(phylum, family)

## 5) Build counts for circos rings ---------------------------
arg_counts_pf <- pf_df3 %>%
  count(phylum, family, name = "freq") %>%
  filter(freq > 0)

# Phylum totals; "Other" last
phylum_totals <- arg_counts_pf %>%
  group_by(phylum) %>%
  summarise(freq = sum(freq), .groups = "drop") %>%
  arrange(if_else(phylum == "Other", Inf, -freq)) %>%
  mutate(
    cum_freq = cumsum(freq),
    xmin = lag(cum_freq, default = 0),
    xmax = cum_freq
  )

# Family ranges inside each phylum sector
family_ranges <- arg_counts_pf %>%
  left_join(phylum_totals %>% select(phylum, xmin, xmax, freq_phylum = freq), by = "phylum") %>%
  arrange(desc(freq_phylum), phylum, desc(freq)) %>%
  group_by(phylum) %>%
  mutate(
    cum_family = cumsum(freq),
    family_xmin = xmin + lag(cum_family, default = 0),
    family_xmax = xmin + cum_family
  ) %>%
  ungroup() %>%
  select(phylum, family, freq, family_xmin, family_xmax)

## 6) Colors ---------------------------------------------------
# Phylum colors
base_phylum_colors <- c(
  "#245D90", "#00827C", "#72B33D", "#997300",
  "#8E4585", "#FCB887", "#EB051F", "#5E4886",
  "#F05C80", "#FCCC00", "#636363"
)

phylum_levels <- phylum_totals$phylum
phylum_palette <- if (length(base_phylum_colors) < length(phylum_levels)) {
  grDevices::colorRampPalette(base_phylum_colors)(length(phylum_levels))
} else {
  base_phylum_colors[seq_along(phylum_levels)]
}
phylum_colors <- setNames(phylum_palette, phylum_levels)
phylum_colors["Other"] <- "#BFBCCC"

# Family colors
family_levels <- unique(family_ranges$family)
family_palette <- grDevices::colorRampPalette(brewer.pal(12, "Paired"))(length(family_levels))
family_colors <- setNames(family_palette, family_levels)
family_colors["Other"] <- "#BFBCCC"

## 7) Plot (draw INSIDE the pdf device) -----------------------
top5_phyla <- phylum_totals$phylum[1:min(5, nrow(phylum_totals))]

pdf(file = "Fig_5_C_1.pdf", width = 6, height = 6)

circos.clear()
circos.par(
  start.degree = 180,
  gap.after = 0.01,
  track.margin = c(0.02, 0.02),
  cell.padding = c(0.02, 0, 0.02, 0)
)

circos.initialize(
  factors = factor(phylum_totals$phylum, levels = phylum_totals$phylum),
  xlim = cbind(phylum_totals$xmin, phylum_totals$xmax)
)

## Outer ring: families
label_threshold_frac <- 0.10  # label if freq > max(freq within phylum) * threshold

circos.track(
  ylim = c(0, 1),
  track.height = 0.15,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    fam_sub <- family_ranges[family_ranges$phylum == sector, , drop = FALSE]
    if (nrow(fam_sub) == 0) return()
    
    max_f <- max(fam_sub$freq, na.rm = TRUE)
    
    for (i in seq_len(nrow(fam_sub))) {
      circos.rect(
        fam_sub$family_xmin[i], 0,
        fam_sub$family_xmax[i], 1,
        col = family_colors[fam_sub$family[i]],
        border = "white"
      )
      
      # Only label families in top5 phyla and above threshold
      if (sector %in% top5_phyla && fam_sub$freq[i] > max_f * label_threshold_frac) {
        xpos <- (fam_sub$family_xmin[i] + fam_sub$family_xmax[i]) / 2
        text_y <- 1.8
        
        circos.text(
          xpos, text_y,
          labels = fam_sub$family[i],
          facing = "downward",
          adj = c(0.5, 0),
          col = "black",
          cex = 0.7
        )
        
        circos.lines(
          c(xpos, xpos),
          c(1, text_y - 0.05),
          col = "black",
          lwd = 0.5
        )
      }
    }
  }
)

## Inner ring: phyla
circos.track(
  ylim = c(0, 1),
  track.height = 0.20,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    circos.rect(
      CELL_META$xlim[1], 0,
      CELL_META$xlim[2], 1,
      col = phylum_colors[sector],
      border = "white"
    )
    
    # Only label top5 phyla
    if (sector %in% top5_phyla) {
      circos.text(
        CELL_META$xcenter, 0.5,
        sector,
        facing = "inside",
        niceFacing = TRUE,
        cex = 0.9,
        col = "white",
        font = 2
      )
    }
  }
)

dev.off()
circos.clear()



## =========================================================
##  Fig 5C-2: Host-known vs host-unknown proportions
##           (Phage with ARG vs All phage)
## =========================================================
## 0) Session setup ------------------------------------------
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## 1) Load data ------------------------------------------------
data_path <- "data/ARG_host_2.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: withhost

if (!exists("withhost")) stop("Object `withhost` was not found after loading: ", data_path)

# Ensure matrix-like object with required dimnames
stopifnot(!is.null(dim(withhost)))
stopifnot(!is.null(rownames(withhost)), !is.null(colnames(withhost)))

required_rows <- c("yes", "no")
required_cols <- c("arg", "all")

missing_rows <- setdiff(required_rows, rownames(withhost))
missing_cols <- setdiff(required_cols, colnames(withhost))
if (length(missing_rows) > 0 || length(missing_cols) > 0) {
  stop(
    "withhost must contain rows {yes,no} and cols {arg,all}. Missing: ",
    paste(c(paste0("rows:", missing_rows), paste0("cols:", missing_cols)), collapse = "; ")
  )
}

## 2) Build plotting dataframe --------------------------------
plot_df <- tibble::tibble(
  phage_group = rep(c("Phage\nwith ARG", "All\nphage"), each = 2),
  host_status = rep(c("Host-Known", "Host-Unknown"), times = 2),
  count = c(
    withhost["yes", "arg"],
    withhost["no",  "arg"],
    withhost["yes", "all"],
    withhost["no",  "all"]
  )
) %>%
  group_by(phage_group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

## 3) Plot -----------------------------------------------------
status_colors <- c("Host-Known" = "#E41A1C", "Host-Unknown" = "#377EB8")

p_fig5c2 <- ggplot(plot_df, aes(x = phage_group, y = percentage, fill = host_status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = status_colors) +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4,
    fontface = "bold"
  ) +
  labs(x = NULL, y = "Percentage (%)", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )

print(p_fig5c2)

## 4) Save -----------------------------------------------------
ggsave(
  filename = "Fig_5_C_2.pdf",
  plot = p_fig5c2,
  width = 3.5,
  height = 5,
  units = "in"
)




## =========================================================
##  Fig 5E: Type ARG heatmap
## =========================================================
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(pheatmap)
})

## 1) Load data ----------------------------------------------
data_path <- "data/type_ARG.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: type_ARG

if (!exists("type_ARG")) stop("Object `type_ARG` was not found after loading: ", data_path)

## 2) Prepare matrix -----------------------------------------
# Transpose so that rows = cities (or vice versa, depending on your original intent)
type_ARG_mat <- as.matrix(t(type_ARG))

# Keep rows with at least one non-zero entry
type_ARG_mat <- type_ARG_mat[rowSums(type_ARG_mat != 0, na.rm = TRUE) > 0, , drop = FALSE]

if (nrow(type_ARG_mat) == 0 || ncol(type_ARG_mat) == 0) {
  stop("No non-zero rows/columns left after filtering; cannot draw heatmap.")
}


## 3) Breaks & colors -----------------------------------------
breaks <- seq(-3, 3, by = 1)
hm_colors <- colorRampPalette(c("lightblue", "white", "darkorange"))(length(breaks) - 1)
## 4) Label formatting ----------------------------------------

type_ARG_mat=type_ARG_mat[,c("sludge","wwtp","bioreactor",
                             "biofilm","Fresh water","soil",
                             "indoor","transportation or  hubs",
                             "outdoor","hospital","indoor air","air")]

col_labels <- gsub("\\.", " ", colnames(type_ARG_mat))

heatmap_type=pheatmap(type_ARG_mat,scale = 'row',
                      color = hm_colors,
                      breaks = breaks,
                      treeheight_row = 0,  # 隐藏行的聚类树状图
                      treeheight_col = 0,
                      cluster_rows = TRUE,
                      cluster_cols = FALSE,
                      angle_col = '45',            # X轴标签旋转45度
                      labels_col = col_labels
) 

## 5) Draw & save ---------------------------------------------
pdf(file='Fig_5_E.pdf',width = 5,height = 5)
print(heatmap_type)
dev.off()



## =========================================================
##  Fig 5F: City ARG heatmap
## =========================================================
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(pheatmap)
})

## 1) Load data ----------------------------------------------
data_path <- "data/city_ARG.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: city_ARG

if (!exists("city_ARG")) stop("Object `city_ARG` was not found after loading: ", data_path)

## 2) Prepare matrix -----------------------------------------
# Transpose so that rows = cities (or vice versa, depending on your original intent)
city_arg_mat <- as.matrix(t(city_ARG))

# Keep rows with at least one non-zero entry
city_arg_mat <- city_arg_mat[rowSums(city_arg_mat != 0, na.rm = TRUE) > 0, , drop = FALSE]

if (nrow(city_arg_mat) == 0 || ncol(city_arg_mat) == 0) {
  stop("No non-zero rows/columns left after filtering; cannot draw heatmap.")
}

## 3) Label formatting ----------------------------------------
col_labels <- gsub("\\.", " ", colnames(city_arg_mat))

## 4) Breaks & colors -----------------------------------------
breaks <- seq(-6, 6, by = 1)
hm_colors <- colorRampPalette(c("lightblue", "white", "darkorange"))(length(breaks) - 1)

## 5) Draw & save ---------------------------------------------
heatmap_city<-pheatmap(
  mat = city_arg_mat,
  scale = "row",
  color = hm_colors,
  breaks = breaks,
  treeheight_row = 0,
  treeheight_col = 0,
  angle_col = '315',               
  labels_col = col_labels
)

pdf(file='Fig_5_F.pdf',width = 10,height = 5)
print(heatmap_city)
dev.off()
