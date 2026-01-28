## =========================================================
##  ComplexHeatmap multi-track annotation figure
## =========================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
})

## ---- Load data ----
data_path <- "data/complexheatmap.RData"
stopifnot(file.exists(data_path))
load(data_path)

# expected: dfnew, class_dist, specialist_dist, lifestyle_dist, uniqueless_dist
stopifnot(exists("dfnew"))
df <- dfnew

## ---- Helpers ----
normalize_matrix_to_percent <- function(mat) {
  mat <- as.matrix(mat)
  # keep dimnames
  rn <- rownames(mat); cn <- colnames(mat)
  out <- t(apply(mat, 1, function(x) {
    total <- sum(x, na.rm = TRUE)
    if (is.na(total) || total == 0) rep(0, length(x)) else x / total * 100
  }))
  rownames(out) <- rn
  colnames(out) <- cn
  out
}

## ---- Normalize distributions ----
class_dist       <- normalize_matrix_to_percent(class_dist)
specialist_dist  <- normalize_matrix_to_percent(specialist_dist)
lifestyle_dist   <- normalize_matrix_to_percent(lifestyle_dist)
uniqueless_dist  <- normalize_matrix_to_percent(uniqueless_dist)

## ---- Ensure specialist columns (optional; only if 3 cols) ---
if (ncol(specialist_dist) == 3) {
  colnames(specialist_dist) <- c("Generalist", "Specialist", "Unknown")
}

## ---- Colors ----
phylum_colors <- c(
  Acidobacteriota   = "#4DBCD7",
  Actinomycetota    = "#F29B7F",
  Bacillota         = "#1CA088",
  Bacillota_A       = "#8491B4",
  Bacteroidota      = "#8C5B87",
  Chloroflexota     = "#D08263",
  Halobacteriota    = "#45609B",
  Methanobacteriota = "#D0423F",
  Nanoarchaeota     = "#2A3E61",
  Patescibacteria   = "#2A8675",
  Planctomycetota   = "#9A60A7",
  Pseudomonadota    = "#247E94",
  Verrucomicrobiota = "#E27C55",
  other__Archaea    = "#BDBDBD",
  other__Bacteria   = "#BDBDBD",
  Unhosted          = "#BDBDBD"
)

class_colors <- c(
  Caudoviricetes     = "#2FB79A",
  Faserviricetes     = "#EF5E5F",
  Malgrandaviricetes = "#3A538A"
)

specialist_colors <- c(
  Generalist = "#AAD9E9",
  Specialist = "#4DBCD7",
  Unknown    = "#BDBDBD"
)

lifestyle_colors <- c(
  Temperate = "#1CA088",
  Virulent  = "#99C9BD",
  Unknown   = "#BDBDBD"
)

uniqueless_colors <- c(
  Yes = "#3B5588",
  No  = "#9295B8"
)

## ---- Validate phylum levels in df ----
stopifnot(all(c("Host_Class", "Host_Phylum", "vOTU", "genome_size") %in% colnames(df)))

# Use ONLY phyla that appear in df, keep order defined by phylum_colors, and drop missing
phylum_in_df <- unique(as.character(df$Host_Phylum))
phylum_levels <- intersect(names(phylum_colors), phylum_in_df)

# If there are phyla not in palette, put them in "Unhosted"/"other__Bacteria" or gray
missing_phyla <- setdiff(phylum_in_df, names(phylum_colors))
if (length(missing_phyla) > 0) {
  # append gray for missing
  phylum_colors <- c(phylum_colors, setNames(rep("#BDBDBD", length(missing_phyla)), missing_phyla))
  phylum_levels <- c(phylum_levels, missing_phyla)
}

df <- df %>%
  mutate(
    Host_Phylum = factor(as.character(Host_Phylum), levels = phylum_levels),
    Host_Class  = as.character(Host_Class)
  ) %>%
  arrange(Host_Phylum, Host_Class)

## ---- Placeholder heatmap matrix ----
mat <- matrix(1, nrow = nrow(df), ncol = 1)
rownames(mat) <- df$Host_Class
colnames(mat) <- "Phage"

## ---- Left annotation: Host class labels ----
left_ha <- rowAnnotation(
  Host_Class = anno_text(df$Host_Class, just = "right", location = 0.5, gp = gpar(fontsize = 7)),
  width = unit(0.1, "cm")
)

## ---- Right annotations (multi-tracks) ----
# Safer genome_size list -> list of numeric vectors
genome_list <- lapply(df$genome_size, function(x) as.numeric(x))
genome_list <- lapply(genome_list, function(x) log10(x))

right_ha <- rowAnnotation(
  Phylum = anno_empty(width = unit(3.5, "cm"), border = TRUE),
  
  gap1 = anno_empty(width = unit(0.1, "cm"), border = FALSE),
  
  vOTU = anno_barplot(
    log10(as.numeric(df$vOTU)),
    gp = gpar(fill = phylum_colors[as.character(df$Host_Phylum)], col = NA),
    width = unit(2, "cm"),
    border = FALSE,
    axis_param = list(
      at = 0:4,
      labels = expression(10^0, 10^1, 10^2, 10^3, 10^4),
      facing = "outside"
    )
  ),
  
  gap2 = anno_empty(width = unit(0.1, "cm"), border = FALSE),
  
  Genome = anno_boxplot(
    genome_list,
    outline = FALSE,
    box_width = 0.6,
    gp = gpar(
      fill = phylum_colors[as.character(df$Host_Phylum)],
      col  = phylum_colors[as.character(df$Host_Phylum)]
    ),
    width = unit(3, "cm"),
    border = FALSE,
    axis_param = list(
      at = 4:6,
      labels = expression(10^4, 10^5, 10^6),
      facing = "outside"
    )
  ),
  
  gap3 = anno_empty(width = unit(0.1, "cm"), border = FALSE),
  
  Class = anno_barplot(
    class_dist,
    gp = gpar(fill = class_colors[colnames(class_dist)], col = NA),
    border = FALSE,
    bar_width = 0.8,
    width = unit(3, "cm"),
    axis_param = list(at = c(0, 50, 100), labels = c("0%", "50%", "100%"), facing = "outside")
  ),
  
  gap4 = anno_empty(width = unit(0.1, "cm"), border = FALSE),
  
  Specialist = anno_barplot(
    specialist_dist,
    gp = gpar(fill = specialist_colors[colnames(specialist_dist)], col = NA),
    border = FALSE,
    bar_width = 0.8,
    width = unit(2, "cm"),
    axis_param = list(at = c(0, 50, 100), labels = c("0%", "50%", "100%"), facing = "outside")
  ),
  
  gap5 = anno_empty(width = unit(0.1, "cm"), border = FALSE),
  
  Lifestyle = anno_barplot(
    lifestyle_dist,
    gp = gpar(fill = lifestyle_colors[colnames(lifestyle_dist)], col = NA),
    border = FALSE,
    bar_width = 0.8,
    width = unit(2, "cm"),
    axis_param = list(at = c(0, 50, 100), labels = c("0%", "50%", "100%"), facing = "outside")
  ),
  
  gap6 = anno_empty(width = unit(0.1, "cm"), border = FALSE),
  
  Uniqueless = anno_barplot(
    uniqueless_dist,
    gp = gpar(fill = uniqueless_colors[colnames(uniqueless_dist)], col = NA),
    border = FALSE,
    bar_width = 0.8,
    width = unit(2, "cm"),
    axis_param = list(at = c(0, 50, 100), labels = c("0%", "50%", "100%"), facing = "outside")
  ),
  
  annotation_name_side = "bottom",
  border = FALSE
)

## ---- Heatmap object ----
ht <- Heatmap(
  mat,
  name = "placeholder",
  col = "white",
  width = unit(0.01, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  show_row_names = FALSE,
  left_annotation = left_ha,
  right_annotation = right_ha,
  row_split = df$Host_Phylum,
  row_title = NULL,
  gap = unit(2, "mm")
)

## ---- Legends ----
lgd_list <- list(
  Legend(title = "Phylum",
         labels = phylum_levels,
         legend_gp = gpar(fill = unname(phylum_colors[phylum_levels])),
         ncol = 3),
  Legend(title = "Class",
         labels = colnames(class_dist),
         legend_gp = gpar(fill = unname(class_colors[colnames(class_dist)]))),
  Legend(title = "Specialist",
         labels = colnames(specialist_dist),
         legend_gp = gpar(fill = unname(specialist_colors[colnames(specialist_dist)]))),
  Legend(title = "Lifestyle",
         labels = colnames(lifestyle_dist),
         legend_gp = gpar(fill = unname(lifestyle_colors[colnames(lifestyle_dist)]))),
  Legend(title = "Uniqueless",
         labels = colnames(uniqueless_dist),
         legend_gp = gpar(fill = unname(uniqueless_colors[colnames(uniqueless_dist)])))
)

## ---- Draw & Save ----
pdf("Fig3.pdf", width = 14, height = 10)
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_list = packLegend(lgd_list, direction = "horizontal"),
  annotation_legend_side = "bottom"
)

## ---- Add left border lines for each split (robust slices) ---
n_groups <- length(unique(df$Host_Phylum))

for (i in seq_len(n_groups)) {
  for (nm in c("Phylum", "vOTU", "Genome", "Class", "Specialist", "Lifestyle", "Uniqueless")) {
    decorate_annotation(nm, slice = i, {
      grid.lines(
        x = unit(0, "npc"), y = unit(c(0, 1), "npc"),
        gp = gpar(col = "black", lwd = 1)
      )
    })
  }
}

## ---- Add phylum text inside Phylum block (only existing splits) ----
split_levels <- levels(df$Host_Phylum)

for (i in seq_along(split_levels)) {
  decorate_annotation("Phylum", slice = i, {
    grid.text(
      label = split_levels[i],
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      just = "center",
      gp = gpar(fontsize = 10, fontface = "bold", col = phylum_colors[split_levels[i]])
    )
  })
}

dev.off()
