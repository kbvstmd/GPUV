
## =========================================================
##  Fig 6A (combined): 
##   - Fig 6A-1: Human-pathogenic virus families (horizontal bar)
##   - Fig 6A-2: Nested donut (outer: phage vs other DNA viruses;
##                            inner: human pathogens vs others)
##  Output: Fig_6_A_combined.pdf
## =========================================================

## 0) Session setup ------------------------------------------
rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggpubr)
  library(ggplot2)
  library(patchwork)
})

## 1) Load data ----------------------------------------------
data_path <- "data/nophage.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: nophage

if (!exists("nophage")) stop("Object `nophage` was not found after loading: ", data_path)
stopifnot(is.data.frame(nophage))
stopifnot("V11" %in% colnames(nophage))

## 2) Summarise family counts ---------------------------------
family_df <- as.data.frame(table(nophage$V11), stringsAsFactors = FALSE) %>%
  rename(family = Var1, count = Freq) %>%
  mutate(family = if_else(is.na(family) | family == "", "Unknown", family))

## 3) Define human pathogen families --------------------------
human_pathogen_families <- c(
  "Adenoviridae", "Herpesviridae", "Papillomaviridae",
  "Poxviridae", "Polyomaviridae", "Parvoviridae"
)

family_df <- family_df %>%
  mutate(group = if_else(family %in% human_pathogen_families, "human pathogens", "others"))

## 4) Fig 6A-1: Horizontal bar (human pathogens only) ----------
bar_df <- family_df %>%
  filter(group == "human pathogens") %>%
  arrange(desc(count)) %>%
  mutate(family = factor(family, levels = rev(family)))

p_fig6a1 <- ggbarplot(
  bar_df,
  x = "family",
  y = "count",
  fill = "group",
  orientation = "horiz",
  xlab = "Family",
  ylab = "Count",
  palette = c("human pathogens" = "#990000", "others" = "#1F4E79"),
  title = ""
) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

## 5) Fig 6A-2: Nested donut ----------------------------------
total_contigs <- 23800

n_non_phage <- nrow(nophage)
n_phage <- total_contigs - n_non_phage
if (n_phage < 0) stop("total_contigs is smaller than nrow(nophage); check `total_contigs`.")

# IMPORTANT: use ONLY human-pathogen families for the inner pie
n_human_pathogen <- sum(family_df$count[family_df$group == "human pathogens"], na.rm = TRUE)
n_other_non_phage <- n_non_phage - n_human_pathogen
if (n_other_non_phage < 0) stop("Human-pathogen counts exceed nrow(nophage); check filtering.")

outer_df <- tibble(
  category = c("Bacteriophage", "Other DNA viruses"),
  value = c(n_phage, n_non_phage)
) %>%
  mutate(
    percent = value / sum(value) * 100,
    label = paste0(category, "\n", sprintf("%.2f%%", percent))
  )

inner_df <- tibble(
  subcategory = c("Human pathogens", "Others"),
  value = c(n_human_pathogen, n_other_non_phage)
) %>%
  mutate(
    percent = value / sum(value) * 100,
    label = paste0(subcategory, "\n", sprintf("%.2f%%", percent), " (", value, ")")
  )

outer_palette <- c(
  "Bacteriophage" = "#4ECDC4",
  "Other DNA viruses" = "#FFA07A"
)
inner_palette <- c(
  "Human pathogens" = "#FF6B6B",
  "Others" = "#45B7D1"
)

outer_ring <- ggdonutchart(
  outer_df,
  x = "value",
  label = "label",
  fill = "category",
  color = "white",
  lab.pos = "in",
  lab.font = c(4, "plain", "black"),
  palette = outer_palette
) +
  theme(legend.position = "none") +
  xlim(0.5, 2.5)

inner_pie <- ggpie(
  inner_df,
  x = "value",
  label = "label",
  fill = "subcategory",
  color = "white",
  lab.pos = "in",
  lab.font = c(3, "plain", "black"),
  palette = inner_palette
) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

p_fig6a2 <- outer_ring +
  inset_element(
    inner_pie,
    left = 0.35, right = 0.65,
    bottom = 0.35, top = 0.65
  )

## 6) Combine into one figure ---------------------------------
# Layout: bar plot on the left, donut on the right
combined_fig6a <- p_fig6a1 + p_fig6a2 +
  plot_layout(widths = c(1.4, 1))

print(combined_fig6a)

## 7) Save ----------------------------------------------------
ggsave(
  filename = "Fig_6_A.pdf",
  plot = combined_fig6a,
  width = 11,
  height = 6,
  units = "in",
  dpi = 300
)


## =========================================================
##  Fig 6B: Prevalence of pathogen families (binary presence)
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggpubr)
  library(ggplot2)
})

## 0) Load data ----------------------------------------------
data_path <- "data/all_type_nophage.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: res_type
res_type_bin <- res_type

## 1) Columns to binarize / count -----------------------------
# 你这里建议把列名写死（更稳），或至少统一范围
cols_bin   <- 5:11   # binarize these columns
cols_count <- 5:10   # compute prevalence from these columns (按你原逻辑)

stopifnot(max(cols_bin) <= ncol(res_type_bin))
stopifnot(max(cols_count) <= ncol(res_type_bin))

## 2) Binarize (presence/absence) -----------------------------
threshold <- 1e-4
res_type_bin <- res_type_bin %>%
  mutate(across(all_of(cols_bin), ~ if_else(.x > threshold, 1L, 0L)))

## 3) Column sums (prevalence) --------------------------------
prev_vec <- colSums(res_type_bin[, cols_count, drop = FALSE], na.rm = TRUE)

prev_df <- enframe(prev_vec, name = "pathogen", value = "count") %>%
  arrange(count) %>%
  mutate(pathogen = factor(pathogen, levels = pathogen))  # control plotting order

## 4) Plot -----------------------------------------------------
p_fig6b <- ggbarplot(
  prev_df,
  x = "pathogen",
  y = "count",
  orientation = "horiz",
  xlab = "Pathogen family",
  ylab = "Count",
  fill = "pathogen",
  palette = "jama",
  title = ""
) +
  theme_pubr() +
  theme(legend.position = "none")

print(p_fig6b)

## 5) Save -----------------------------------------------------
ggsave(
  filename = "Fig_6_B.pdf",
  plot = p_fig6b,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)


## =========================================================
##  Fig 6C: Violin plot of pathogen relative abundance by environment
## =========================================================
rm(list = ls(all.names = TRUE))


suppressPackageStartupMessages({
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

## 0) Load data ----------------------------------------------
data_path <- "data/all_type_nophage.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: res_type
plot_df <- res_type

required_cols <- c("type", "pathogens")
missing_cols <- setdiff(required_cols, colnames(plot_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in `res_type`: ", paste(missing_cols, collapse = ", "))
}

plot_df <- plot_df %>%
  filter(!is.na(type), type != "", !is.na(pathogens))

## 1) Transform ------------------------------------------------
# Recommended, interpretable transform for abundances:
# log10(x + 1) avoids -Inf and doesn't require arbitrary eps.
plot_df <- plot_df %>% mutate(pathogens_log10 = log10(pathogens + 1e-12))

x_col <- "type"
y_col <- "pathogens_log10"

## 2) Order groups by median ----------------------------------
type_levels <- plot_df %>%
  group_by(.data[[x_col]]) %>%
  summarise(med = median(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(.data[[x_col]])

plot_df[[x_col]] <- factor(plot_df[[x_col]], levels = type_levels)

## 3) Colors ---------------------------------------------------
# If you already have `colors` defined globally as a named vector, it will be used.
# Otherwise, define a default palette here.
colors <- c(
  "air"                     = "#0B3C5D",
  "hospital"                = "#0E6F6D",
  "transportation or  hubs" = "#8B0000",
  "indoor"                  = "#D4A017",
  "outdoor"                 = "#1B7F2A",
  "indoor air"              = "#8A9AA5",
  "Fresh water"             = "#5A189A",
  "sludge"                  = "#D35400",
  "biofilm"                 = "#6CA0DC",
  "bioreactor"              = "#7BA23F",
  "coast"                   = "#7A001F",
  "soil"                    = "#2F3E46",
  "wwtp"                    = "#3DD5C2"
)
# Ensure all levels have a color; add gray fallback if missing
missing_cols <- setdiff(levels(plot_df[[x_col]]), names(colors))
if (length(missing_cols) > 0) {
  fallback <- rep("#BDBDBD", length(missing_cols))
  names(fallback) <- missing_cols
  colors <- c(colors, fallback)
}

## 4) Plot -----------------------------------------------------
p_fig6c <- ggviolin(
  plot_df,
  x = x_col, y = y_col,
  fill = x_col, color = x_col,
  add = "boxplot",
  add.params = list(fill = "white", color = "black", width = 0.1),
  xlab = "Environment",
  ylab = "log10(pathogens relative abundance + 1)"
) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  stat_compare_means(
    method = "anova",
    label.x.npc = "left",
    label.y.npc = "top",
    aes(label = paste0("ANOVA p = ", ..p.format..))
  )

print(p_fig6c)

## 5) Save -----------------------------------------------------
ggsave(
  filename = "Fig_6_C.pdf",
  plot = p_fig6c,
  device = "pdf",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)


## =========================================================
##  Fig 6D: Heatmap of human pathogen families across types
## =========================================================
rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
})

## 1) Load data ----------------------------------------------
data_path <- "data/all_type_nophage.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: res_type

if (!exists("res_type")) stop("Object `res_type` was not found after loading: ", data_path)
stopifnot(is.data.frame(res_type))
stopifnot("type" %in% colnames(res_type))

## 2) Summarise by type --------------------------------------
# Ensure the target columns exist
target_cols <- c("Adenoviridae", "Herpesviridae", "Papillomaviridae",
                 "Poxviridae", "Polyomaviridae", "Parvoviridae")
missing_cols <- setdiff(target_cols, colnames(res_type))
if (length(missing_cols) > 0) {
  stop("Missing pathogen-family columns in `res_type`: ", paste(missing_cols, collapse = ", "))
}

type_means <- res_type %>%
  group_by(type) %>%
  summarise(across(all_of(target_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# type_means: rows = type, cols = families
type_means_df <- as.data.frame(type_means)
rownames(type_means_df) <- type_means_df$type
type_means_df$type <- NULL

# Heatmap matrix wants rows = families, cols = type (your original transpose)
heat_mat <- t(as.matrix(type_means_df))

## 3) Row annotation -----------------------------------------
row_annotation <- data.frame(
  Transmission = c(
    "Respiratory Secretion",
    "Direct Contact",
    "Skin Contact",
    "Direct Contact",
    "Direct Contact",
    "Respiratory Secretion"
  ),
  row.names = rownames(heat_mat)
)

# Define annotation colors explicitly (do NOT rely on `colors` existing)
transmission_colors <- c(
  "Respiratory Secretion" = "#1F4E79",
  "Skin Contact"          = "#005B5C",
  "Direct Contact"        = "#990000"
)


ann_colors <- list(Transmission = transmission_colors)

## 4) Draw & save heatmap ------------------------------------
# For z-scored rows, a symmetric break range is usually best
breaks <- seq(-3, 3, length.out = 101)
hm_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(length(breaks) - 1)

pathogens_heatmap<-pheatmap(
  mat = heat_mat,
  scale = "row",
  color = hm_colors,
  breaks = breaks,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  fontsize_col = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  cutree_rows = 3,
  treeheight_row = 0,
  treeheight_col = 0,
)


# 保存图表为PDF文件
pdf(file="Fig_6_D.pdf",width = 8,height = 6)
print(pathogens_heatmap)
dev.off()

