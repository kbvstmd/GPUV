
########################################################
#############Fig_2_A-F
########################################################

## =========================================================
##  Alpha diversity plots (box/violin) - FULL SCRIPT
##  - Load data from: all_type_alpha.RData
##  - Output: PDF files for Shannon/Simpson
## =========================================================
rm(list = ls(all.names = TRUE))
## 0) Packages ------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(rlang)
  library(ggplot2)
  library(ggpubr)
})

## 1) Load data ----------------------------------------------
data_path <- "data/all_type_alpha.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: plot_data

if (!exists("plot_data")) {
  stop("Object `plot_data` was not found after loading: ", data_path)
}
stopifnot(is.data.frame(plot_data))

required_cols <- c("type", "shannon", "simpson")
missing_cols <- setdiff(required_cols, colnames(plot_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in `plot_data`: ", paste(missing_cols, collapse = ", "))
}

## 2) Color palette (define once) -----------------------------
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

## 3) Helper: save PDF via ggsave -----------------------------
save_pdf <- function(plot, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
}

## 4) Plot function: violin + inner box + ANOVA ---------------
plot_alpha_violin_anova <- function(data, x_col, y_col, colors) {
  stopifnot(is.data.frame(data))
  stopifnot(x_col %in% names(data), y_col %in% names(data))
  
  # Order x levels by median(y) desc
  medians <- data %>%
    group_by(.data[[x_col]]) %>%
    summarise(median_y = median(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(median_y))
  
  data[[x_col]] <- factor(data[[x_col]], levels = medians[[x_col]])
  
  ggpubr::ggviolin(
    data,
    x = x_col, y = y_col,
    fill = x_col, color = x_col,
    add = "boxplot",
    add.params = list(fill = "white", color = "black", width = 0.1)
  ) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggpubr::stat_compare_means(
      method = "anova",
      label.x.npc = "left",
      label.y.npc = "top",
      aes(label = paste0("ANOVA p = ", ..p.format..))
    )
}


## 6) Batch plotting ------------------------------------------
alpha_metrics <- c("shannon")
x_col <- "type"

for (metric in alpha_metrics) {
  y_col <- metric
  message("Plotting metric: ", y_col)
  
  if (!y_col %in% names(plot_data)) {
    warning("Skipping metric (column not found): ", y_col)
    next
  }
  
  # Violin + ANOVA label
  p_violin_nop <- plot_alpha_violin_anova(plot_data, x_col, y_col, type_colors)
  save_pdf(p_violin_nop,  "Fig_2_A.pdf", width = 12, height = 6)
}


########################################################
#############Fig_2_B
########################################################


## =========================================================
##  Fig 2B1: PCoA scatter (Bray–Curtis) colored by sample type
## =========================================================

## 0) Session setup ------------------------------------------

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## 1) Load data ----------------------------------------------
dist_path <- "data/all_type_PCoA.RData"

stopifnot(file.exists(dist_path))
load(dist_path)  # expected: PCoA

## 2) Build PCoA coordinates ---------------------------------
pcoa_scores <- as.data.frame(PCoA$points[, 1:2, drop = FALSE])
pcoa_scores$samples <- rownames(pcoa_scores)

colnames(pcoa_scores)[1:2] <- c("PCoA1", "PCoA2")

# Explained variance (optional; not used in your original plot, but kept)
pcoa_var_pct <- round(PCoA$eig / sum(PCoA$eig) * 100, digits = 2)

## 3) Join with metadata -------------------------------------
# Keep only samples present in both PCoA and metadata
common_samples <- intersect(pcoa_scores$samples, metadata$Run)

plot_df <- pcoa_scores %>%
  filter(samples %in% common_samples) %>%
  inner_join(
    metadata %>%
      select(Run, type) %>%
      rename(samples = Run, group = type),
    by = "samples"
  ) %>%
  filter(!is.na(group), group != "")

## 4) Filter groups with enough samples ----------------------
min_n_per_group <- 30

group_counts <- plot_df %>%
  count(group, name = "n") %>%
  arrange(desc(n)) %>%
  filter(n > min_n_per_group)

plot_df <- plot_df %>%
  filter(group %in% group_counts$group) %>%
  mutate(group = factor(group, levels = group_counts$group))

## 5) Color palette ------------------------------------------
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

# Ensure all groups in the plot have a color; otherwise add a neutral fallback
missing_colors <- setdiff(levels(plot_df$group), names(type_colors))
if (length(missing_colors) > 0) {
  fallback <- rep("#BDBDBD", length(missing_colors))
  names(fallback) <- missing_colors
  type_colors <- c(type_colors, fallback)
}

## 6) Plot ----------------------------------------------------
p <- ggplot(plot_df, aes(x = PCoA1, y = PCoA2, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = type_colors) +
  labs(
    x = "PCoA1",
    y = "PCoA2"
    # If you want variance in axis labels, use:
    # x = paste0("PCoA1 (", pcoa_var_pct[1], "%)"),
    # y = paste0("PCoA2 (", pcoa_var_pct[2], "%)")
  ) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, angle = 90),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 16)
  )

## 7) Save ----------------------------------------------------
ggsave(
  filename = "Fig_2_B1.pdf",
  plot = p,
  width = 18,
  height = 15,
  units = "in"
)



## =========================================================
##  Fig 2B2: ANOSIM (Bray–Curtis) distance-rank boxplots
## =========================================================

## 0) Session setup ------------------------------------------
# NOTE: Avoid clearing the whole environment in scripts.
# rm(list = ls(all.names = TRUE))
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## 1) Load data ----------------------------------------------
data_path <- "data/all_type_anosim.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: df_anosim

if (!exists("df_anosim")) {
  stop("Object `df_anosim` was not found after loading: ", data_path)
}

required_cols <- c("class.vec", "dis.rank", "statistic", "signif")
missing_cols <- setdiff(required_cols, names(df_anosim))
if (length(missing_cols) > 0) {
  stop("Missing required fields in `df_anosim`: ", paste(missing_cols, collapse = ", "))
}

## 2) Build plotting dataframe --------------------------------
plot_df <- tibble::tibble(
  class = as.character(df_anosim$class.vec),
  distance_rank = df_anosim$dis.rank
)


## 3) Order classes by median distance rank -------------------
median_df <- plot_df %>%
  group_by(class) %>%
  summarise(median_rank = median(distance_rank, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_rank))

ordered_levels <- median_df$class

# Put "Between" first if it exists (robust, no index math)
if ("Between" %in% ordered_levels) {
  ordered_levels <- c("Between", setdiff(ordered_levels, "Between"))
}

plot_df$class <- factor(plot_df$class, levels = ordered_levels)

## 4) Color palette -------------------------------------------
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

# If there are classes not present in the palette (e.g., "Between"),
# ggplot will drop fill for those unless we add a default.
# Here we add a neutral color for "Between" if needed.
if ("Between" %in% levels(plot_df$class) && !"Between" %in% names(type_colors)) {
  type_colors <- c(type_colors, "Between" = "#BDBDBD")
}

## 5) Labels ---------------------------------------------------
x_label <- sprintf(
  "R = %s, p = %s",
  signif(df_anosim$statistic, 3),
  signif(df_anosim$signif, 3)
)

## 6) Plot -----------------------------------------------------
p <- ggplot(plot_df, aes(x = class, y = distance_rank)) +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.8) +
  geom_boxplot(aes(fill = class),
               outlier.colour = "white",
               linewidth = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(
    title = "Bray–Curtis ANOSIM",
    x = x_label,
    y = "Rank of Distance (Bray–Curtis)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

## 7) Save -----------------------------------------------------
ggsave(
  filename = "Fig_2_B2.pdf",
  plot = p,
  width = 12,
  height = 8,
  units = "in"
)




########################################################
#############Fig_2_C
########################################################



## =========================================================
##  Pairwise overlap test (hypergeometric) + network plot
##  Clean / standardized naming (FULL SCRIPT)
## =========================================================
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(reshape2)
  library(igraph)
  library(RColorBrewer)
})

## 0) Load data ------------------------------------------
data_path <- "data/network_type.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: data,metadata


## 1) Pairwise hypergeometric test function -------------------
pair_overlap_test <- function(group_a, group_b, pa_table, verbose = FALSE) {
  stopifnot(group_a %in% names(pa_table), group_b %in% names(pa_table))
  
  a <- pa_table[[group_a]]
  b <- pa_table[[group_b]]
  
  # remove rows with NA in either
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]; b <- b[ok]
  
  both_present   <- sum(a == 1 & b == 1)
  a_only_present <- sum(a == 1 & b == 0)
  b_only_present <- sum(a == 0 & b == 1)
  neither        <- sum(a == 0 & b == 0)
  
  # Standard 2x2 contingency table (rows: A present/absent; cols: B present/absent)
  contingency <- matrix(
    c(both_present, b_only_present,
      a_only_present, neither),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      A = c("Present", "Absent"),
      B = c("Present", "Absent")
    )
  )
  
  # Hypergeometric test:
  # N = total
  # K = number of successes in population (= B present)
  # n = draws (= A present)
  # x = observed successes in draws (= both present)
  N <- both_present + a_only_present + b_only_present + neither
  K <- both_present + b_only_present     # B present
  n <- both_present + a_only_present     # A present
  x <- both_present
  
  p_value <- phyper(q = x - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
  
  # Jaccard index
  denom <- (both_present + a_only_present + b_only_present)
  jaccard <- if (denom == 0) 0 else both_present / denom
  
  # Keep your original "normalized_JI" logic (as-is, but guarded)
  sample_a_size <- n + neither
  sample_b_size <- K - both_present + b_only_present + neither  # equals (B present + neither)
  # Safer: B present + neither = (both_present + b_only_present) + neither
  sample_b_size <- (both_present + b_only_present) + neither
  
  normalized_jaccard <- if (N == 0) 0 else jaccard * (sample_a_size * sample_b_size) / (N^2)
  
  if (verbose) {
    message("Pair: ", group_a, " vs ", group_b)
    print(contingency)
    message("Jaccard = ", signif(jaccard, 4),
            "; normalized = ", signif(normalized_jaccard, 4),
            "; p = ", signif(p_value, 4))
  }
  
  data.frame(
    A = group_a,
    B = group_b,
    both = both_present,
    a_only = a_only_present,
    b_only = b_only_present,
    neither = neither,
    jaccard = jaccard,
    normalized_jaccard = normalized_jaccard,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

## 2) Compute pairwise matrices -------------------------------
groups <- colnames(data)

score_mat <- matrix(0, nrow = length(groups), ncol = length(groups),
                    dimnames = list(groups, groups))
pvalue_mat <- matrix(1, nrow = length(groups), ncol = length(groups),
                     dimnames = list(groups, groups))

pair_results <- list()

idx <- 1L
for (i in seq_len(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    res <- pair_overlap_test(groups[i], groups[j], data, verbose = FALSE)
    
    # store only upper triangle
    score_mat[i, j]  <- max(0, res$normalized_jaccard)
    pvalue_mat[i, j] <- min(1, res$p_value)
    
    pair_results[[idx]] <- res
    idx <- idx + 1L
  }
}

pair_df <- bind_rows(pair_results)

## 3) Long table + filtering ----------------------------------
# Convert upper triangle to long format
library(tidyr)
library(dplyr)

score_long <- score_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("A") %>%
  pivot_longer(
    cols = -A,
    names_to = "B",
    values_to = "score"
  )

pval_long <- pvalue_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("A") %>%
  pivot_longer(
    cols = -A,
    names_to = "B",
    values_to = "p_value"
  )
result_long <- score_long %>%
  inner_join(pval_long, by = c("A", "B")) %>%
  filter(A != B) %>%
  # keep only upper triangle (avoid duplicate edges)
  mutate(is_upper = match(A, groups) < match(B, groups)) %>%
  filter(is_upper) %>%
  select(A, B, score, p_value)

# Thresholds (your original)
score_cutoff <- 0.2
p_cutoff <- 0.05

edges_df <- result_long %>%
  filter(score >= score_cutoff, p_value < p_cutoff)

## 4) Build node table (from metadata$type) -------------------
node_df <- as.data.table(table(metadata$type))
setnames(node_df, c("node", "sample_size"))

# Assign meta-groups (your original mapping, cleaned)
node_df[, super_group := "wwtp"]
node_df[node == "air", super_group := "air"]
node_df[node == "indoor air", super_group := "air"]
node_df[node %in% c("coast", "Fresh water"), super_group := "water"]
node_df[node %in% c("indoor", "outdoor", "transportation or  hubs", "hospital"), super_group := "BE"]
node_df[node == "soil", super_group := "soil"]

## 5) Network plot --------------------------------------------
# igraph expects columns named: from/to (or first two columns)
edges_igraph <- edges_df %>%
  rename(from = A, to = B) %>%
  as.data.frame()

nodes_igraph <- node_df %>%
  rename(name = node, group = super_group) %>%
  as.data.frame()

g <- graph_from_data_frame(d = edges_igraph, vertices = nodes_igraph, directed = FALSE)

# Vertex size: log2(sample_size + 1)
V(g)$size <- log2(V(g)$sample_size + 1) * 2

# Vertex colors by group
grp_levels <- sort(unique(V(g)$group))
pal <- brewer.pal(n = max(3, length(grp_levels)), name = "Set2")[seq_along(grp_levels)]
group_colors <- setNames(pal, grp_levels)
V(g)$color <- group_colors[V(g)$group]

# Edge width by score
E(g)$width <- E(g)$score * 5
E(g)$color <- adjustcolor("grey40", alpha.f = 0.6)

set.seed(123)
layout_fr <- layout_with_fr(g)

pdf("Fig_2_C.pdf", width = 8, height = 8)
plot(
  g,
  layout = layout_fr,
  vertex.label = V(g)$name,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  main = "Network of shared viruses"
)
legend(
  "topleft",
  legend = names(group_colors),
  col = group_colors,
  pch = 19,
  pt.cex = 1.5,
  bty = "n"
)
dev.off()


########################################################
#############Fig_2_D
########################################################
## =========================================================
##  Fig 2D: City-level sharing network on world map
##  - Hypergeometric overlap test + (optional) normalized Jaccard
##  - BH-FDR correction across all pairs
## =========================================================

rm(list = ls(all.names = TRUE))


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(maps)        # for map_data("world")
  library(readr)       # for read_csv
  library(scales)      # for rescale()
})

## 0) Inputs --------------------------------------------------
# Required objects:
# - data: presence/absence table (rows: features/contigs; cols: cities/groups), values 0/1
# - city_coordinates2.csv: columns City, lat, lon
data_path <- "data/city_shared_viral.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected: data, metadata

geolocations=geolocations2
## function  hypergeom_test
#
hypergeom_test <- function(group1, group2, data) {
  both <- sum(data[[group1]] == 1 & data[[group2]] == 1)
  only_group1 <- sum(data[[group1]] == 1 & data[[group2]] == 0)
  only_group2 <- sum(data[[group1]] == 0 & data[[group2]] == 1)
  neither <- sum(data[[group1]] == 0 & data[[group2]] == 0)
  
  
  contingency_table <- matrix(c(both, only_group1, only_group2, neither), nrow = 2)
  colnames(contingency_table) <- c(paste0(group1, "_only"), paste0(group2, "_only"))
  rownames(contingency_table) <- c("Present", "Absent")
  print(contingency_table)
  total <- both + only_group1 + only_group2 + neither
  m <- both + only_group1  
  n <- total - m           
  k <- both + only_group2  
  x <- both                
  
  # phyper p
  p_value <- phyper(x - 1, m, n, k, lower.tail = FALSE)
  JI=x/(m+k-x)
  print(JI)
  #  Jaccard
  sample1_size <- m + neither
  sample2_size <- k + neither
  
  # normalize Jaccard 
  normalized_JI <- JI * ( as.numeric(sample1_size) * as.numeric(sample2_size)) / total^2
  result <- data.frame(
    JI = JI,
    normalized_JI=normalized_JI,
    p_value = p_value
  )
  
  return(result)
}

groups <- colnames(data)

result_matrix <- matrix(NA, ncol = length(groups), nrow = length(groups))
colnames(result_matrix) <- groups
rownames(result_matrix) <- groups
result_matrix_pvalue=result_matrix

for (i in 1:(length(groups) - 1)) {
  print(i)
  for (j in (i + 1):length(groups)) {
    result <- hypergeom_test(groups[i], groups[j], data)
    
    if(result$p_value<0.05&result$normalized_JI>=0.2){
      result_matrix[i, j]<-result$normalized_JI
      result_matrix_pvalue[i, j]<-result$p_value
    }
    
  }
}
library(reshape2)

long_result_matrix <- reshape2::melt(as.matrix(result_matrix), 
                                     varnames = c("from", "to"), 
                                     value.name = "value")

long_result_matrix=long_result_matrix[!is.na(long_result_matrix$value),]


cts=data.frame(table(c(long_result_matrix$from,long_result_matrix$to)))

long_result_matrix_pvalue <- reshape2::melt(
  as.matrix(result_matrix_pvalue),
  varnames = c("from", "to"),
  value.name = "value"
)
long_result_matrix_pvalue=long_result_matrix_pvalue[!is.na(long_result_matrix_pvalue$value),]

long_result_matrix_pvalue["qval"]=p.adjust(long_result_matrix_pvalue$value)

long_result_matrix2=long_result_matrix[which(long_result_matrix$value>0.4),]
cts=data.frame(table(c(long_result_matrix2$from,long_result_matrix2$to)))

edges=long_result_matrix2
colnames(edges)[3]='weight'

# 
edges <- edges %>%
  left_join(geolocations, by = c("from" = "City")) %>%
  rename(lat_from = lat, lon_from = lon) %>%
  left_join(geolocations, by = c("to" = "City")) %>%
  rename(lat_to = lat, lon_to = lon)

# world map
world_map <- map_data("world")
p_city_map<-ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = geolocations, aes(x = lon, y = lat), color = "red", size = 2) +
  geom_curve(data = edges, aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to, size = weight*100, color = weight, alpha = 1/weight),
             curvature = 0.2) +
  scale_size(range = c(0.5, 2)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +  
  scale_alpha(range = c(0.1, 1)) + 
  theme_minimal() +
  labs(title = "Distribution of Viral Contigs",
       subtitle = "Shared Viral Contigs Between Cities",
       x = "Longitude", y = "Latitude")

p_city_map


ggsave(
  filename = "Fig_2_D.pdf",
  plot = p_city_map,
  width = 15,
  height = 10,
  units = "in"
)


########################################################
#############Fig_2_E
########################################################

## =========================================================
##  Fig 2E: Heatmap of mean abundance by environment
## =========================================================
rm(list = ls(all.names = TRUE))
suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
})

## 1) Load data ------------------------------------------------
data_path <- "data/Core_Virome_Abundance.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected: Core_Virome_Abundance, metadata

if (!exists("Core_Virome_Abundance")) stop("Missing object: Core_Virome_Abundance")
if (!exists("metadata")) stop("Missing object: metadata")
stopifnot("type" %in% colnames(metadata))

## 2) Compute mean abundance per environment -------------------
# Core_Virome_Abundance: samples x features (assumed)
# metadata$type: environment label per sample

env_means <- aggregate(Core_Virome_Abundance, by=list(Environment=metadata$type), mean)
rownames(env_means) <- env_means$Environment
env_means=env_means[-c(1,3),]
env_means=env_means[,-1]
normalized_env_means <- apply(env_means, 2, function(x) x / sum(x))


# pheatmap
colors <- colorRampPalette(c("#003366", "#F0F0F0", "#990000"))(10)
p<-pheatmap(env_means ,show_rownames = T,show_colnames = F,
            color = colors,
            scale = 'column',
            clustering_method='ward.D2',
            cluster_rows = T,
            cluster_cols = T)

# save pdf
pdf("Fig_2_E.pdf", width = 10, height = 6)
print(p)
dev.off()

########################################################
#############Fig_2_F
########################################################

## =========================================================
## Fig 2F: XGBoost multi-class model with Boruta feature selection
## =========================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(xgboost)
  library(Boruta)
  library(dplyr)
  library(pROC)
  library(ggplot2)
})

set.seed(42)

## 1) Load data ------------------------------------------------
data_path <- "data/city_model.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected objects: otu, metadata

## 2) Prepare dataset ------------------------------------------
label_key <- "city"

# Align OTU table rows with metadata$Run
x <- otu[metadata$Run, , drop = FALSE]

# Combine predictors and label
df <- cbind(x, metadata[, label_key, drop = FALSE])
colnames(df)[ncol(df)] <- "label"

# Convert label to factor, then to 0-based integer class index
df$label <- factor(df$label)
df$label <- as.integer(df$label) - 1L

# Create label mapping (original factor levels -> class id)
id_map <- data.frame(
  label = levels(factor(metadata[, label_key])),
  class = seq_along(levels(factor(metadata[, label_key]))) - 1L,
  stringsAsFactors = FALSE
)

cat("Class distribution (0-based):\n")
print(table(df$label))

## 3) Boruta feature selection --------------------------------
# Note: Boruta expects a data.frame; label should be a factor for classification
df_boruta <- df
df_boruta$label <- factor(df_boruta$label)
#load('data/boruta_fit.RData')
boruta_fit <- Boruta(label ~ ., data = df_boruta, doTrace = 2)

confirmed_features <- getSelectedAttributes(boruta_fit, withTentative = FALSE)
feature_stats <- attStats(boruta_fit)

cat("Number of confirmed features:", length(confirmed_features), "\n")

df_sel <- df[, confirmed_features, drop = FALSE]
df_sel$label <- df$label

## 4) Train/test split -----------------------------------------
train_idx <- sample(seq_len(nrow(df_sel)), size = floor(0.8 * nrow(df_sel)))
train_df <- df_sel[train_idx, , drop = FALSE]
test_df  <- df_sel[-train_idx, , drop = FALSE]

train_y <- train_df$label
test_y  <- test_df$label

train_x <- as.matrix(train_df[, setdiff(colnames(train_df), "label"), drop = FALSE])
test_x  <- as.matrix(test_df[, setdiff(colnames(test_df), "label"), drop = FALSE])

dtrain <- xgb.DMatrix(data = train_x, label = train_y)
dtest  <- xgb.DMatrix(data = test_x,  label = test_y)

n_classes <- length(unique(train_y))

## 5) XGBoost model --------------------------------------------
params <- list(
  objective = "multi:softprob",
  num_class = n_classes,
  eval_metric = "mlogloss",
  eta = 0.1,
  max_depth = 3
)

watchlist <- list(train = dtrain, test = dtest)

xgb_fit <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 50,
  watchlist = watchlist,
  verbose = 1
)

## 6) Predict probabilities ------------------------------------
pred_probs <- predict(xgb_fit, dtest)

# Convert the prediction vector to a probability matrix (n_test x n_classes)
prob_mat <- matrix(pred_probs, ncol = n_classes, byrow = TRUE)
colnames(prob_mat) <- as.character(sort(unique(train_y)))  # "0", "1", "2", ...

# Ensure test labels are factor with matching levels
test_y_factor <- factor(test_y, levels = sort(unique(train_y)))

## 7) ROC per class + mean AUC ---------------------------------
# Colors (keep yours)
colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#a4133c", "#90a955", "#00A6FB",
  "#0582CA", "#006494", "#c9e4ca", "#31572c", "#1F78B4",
  "#c9184a", "#ff4d6d", "#FF5733", "#9933CC", "#FF6699",
  "#339966", "#996633", "#663333", "#FF3366", "#99CCFF",
  "#33CC99", "#FF99CC", "#333399", "#FF9933", "#66CC66",
  "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2",
  "#F39B7FB2", "#8491B4B2", "#0073C2CC", "#EFC000CC",
  "#868686CC", "#CD534CCC", "#7AA6DCCC", "#003C67CC",
  "#8F7700CC", "#3B3B3BCC", "#A73030CC", "#4A6990CC"
)

# One-vs-rest ROC for each class
roc_list <- lapply(levels(test_y_factor), function(cls) {
  roc(response = (test_y_factor == cls), predictor = prob_mat[, cls], quiet = TRUE)
})
names(roc_list) <- levels(test_y_factor)

auc_values <- sapply(roc_list, function(r) as.numeric(r$auc))
mean_auc <- mean(auc_values)

# 95% CI using normal approximation across class AUCs (same as your logic)
se_auc <- sd(auc_values) / sqrt(length(auc_values))
ci_lower <- mean_auc - 1.96 * se_auc
ci_upper <- mean_auc + 1.96 * se_auc

cat(sprintf("Mean AUC = %.3f (95%% CI: %.3f–%.3f)\n", mean_auc, ci_lower, ci_upper))

## 8) Plot ROC curves + mean ROC -------------------------------
roc_plot <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    title = "",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    plot.title   = element_text(size = 16)
  )

# Add per-class ROC curves
for (i in seq_along(roc_list)) {
  r <- roc_list[[i]]
  roc_df <- data.frame(
    fpr = 1 - r$specificities,
    tpr = r$sensitivities
  )
  roc_plot <- roc_plot +
    geom_line(data = roc_df, aes(x = fpr, y = tpr), color = colors[i], size = 1) +
    ylim(0, 1)
}

# Interpolate TPRs onto a common FPR grid to compute mean ROC curve
fpr_grid <- seq(0, 1, length.out = 100)
tpr_interp <- sapply(roc_list, function(r) {
  approx(x = 1 - r$specificities, y = r$sensitivities, xout = fpr_grid)$y
})
mean_tpr <- rowMeans(tpr_interp, na.rm = TRUE)

mean_roc_df <- data.frame(fpr = fpr_grid, tpr = mean_tpr)

roc_plot <- roc_plot +
  geom_line(data = mean_roc_df, aes(x = fpr, y = tpr), color = "black", size = 1.3) +
  annotate(
    "text", x = 0.6, y = 0.1,
    label = sprintf("Mean AUC = %.3f (95%% CI: %.3f–%.3f)", mean_auc, ci_lower, ci_upper),
    color = "black", size = 5, hjust = 0
  )

## 9) Save figure ----------------------------------------------
pdf(file = "Fig_2_F.pdf", width = 10, height = 10)
print(roc_plot)
dev.off()
