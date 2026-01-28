
########################################################
#############Fig_1_A-E
########################################################


########################################################
#############Fig_1_A
########################################################
# Clean environment
rm(list = ls(all = TRUE))

# Set working directory
setwd("D:/R/1008/code")

# Load required packages
library(ggplot2)
library(ggforce)       # For drawing pie charts
library(dplyr)
library(tidyr)
library(sf)            # For spatial data handling
library(rnaturalearth) # For world map data

# Load geographic data
load("data/geolocations.RData")

# Convert data to long format
data_long <- data %>%
  pivot_longer(
    cols = -c(City, lon, lat),
    names_to = "Category",
    values_to = "Value"
  ) %>%
  group_by(City) %>%
  mutate(
    Total = sum(Value),             # Total value per city
    Fraction = Value / Total,       # Fraction for each category
    SizeCategory = case_when(       # Assign size class based on total
      Total <= 10 ~ "1-10",
      Total <= 100 ~ "10-100",
      Total <= 300 ~ "100-300",
      Total > 300 ~ "300+",
      TRUE ~ "Other"
    )
  )

# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define pie radius by total value range
data_long <- data_long %>%
  mutate(
    r = case_when(
      SizeCategory == "1-10" ~ 1,
      SizeCategory == "10-100" ~ 3,
      SizeCategory == "100-300" ~ 5,
      SizeCategory == "300+" ~ 7,
      TRUE ~ 1
    )
  )

# Plot global city-level pie map
p_all_world_city <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "white") +  # Base world map
  geom_arc_bar(
    data = data_long,
    aes(
      x0 = lon, y0 = lat,          # Center of each pie
      r0 = 0, r = r,               # Pie radius
      amount = Fraction,           # Fraction for each slice
      fill = Category
    ),
    stat = "pie"
  ) +
  scale_fill_manual(values = c(
    "transportation.or..hubs" = "#F47F60",
    "wwtp" = "#42A5F5",
    "air" = "#FFB74D",
    "biofilm" = "#6D7EC9",
    "bioreactor" = "#D3D656",
    "coast" = "#9575CD",
    "Fresh.water" = "#26A69A",
    "hospital" = "#DD6895",
    "indoor" = "#F48FB1",
    "indoor.air" = "#4CAF50",
    "outdoor" = "#26C6DA",
    "sludge" = "#7E57C2",
    "soil" = "#C5E1A5",
    "hot spring" = "#556B2F",
    "Game animal" = "#A9A9A9",
    "Green land" = "#D4A017",
    "built env" = "#CD5C5C"
  )) +
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.spacing.x = unit(0.5, "cm")
  ) +
  # Custom legend for city size categories
  geom_point(
    data = data.frame(
      SizeCategory = c("1-10", "10-100", "100-300", "300+"),
      x = c(0, 1, 2, 3),
      y = rep(90, 4),
      size = c(1, 3, 5, 7)
    ),
    aes(x = x, y = y, size = size, color = SizeCategory),
    shape = 16
  ) +
  scale_size_identity() +
  scale_color_manual(values = c(
    "1-10" = "black",
    "10-100" = "gray",
    "100-300" = "darkgray",
    "300+" = "lightgray"
  )) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

# Save to PDF
pdf(file = "Fig_1_A.pdf", width = 15, height = 10)
print(p_all_world_city)
dev.off()



########################################################
#############Fig_1_B
########################################################

rm(list = ls(all = TRUE))

## =========================================================
##  Genome length analysis + plots (EN, full script)
## =========================================================

## 0) Global options -----------------------------------------
# NOTE: Avoid rm(list = ...) in scripts; it can break interactive debugging
set.seed(123)

## 1) Load data ----------------------------------------------
data_path <- "data/density_all.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected objects: density_ref, density_urban, dt1, dt2, dt3

required_objs <- c("density_ref", "density_urban", "dt1", "dt2", "dt3")
missing_objs <- setdiff(required_objs, ls())
if (length(missing_objs) > 0) {
  stop("Missing objects after load(): ", paste(missing_objs, collapse = ", "))
}

## 2) Helper functions ---------------------------------------
sample_rows <- function(df, n) {
  stopifnot(is.data.frame(df))
  if (n > nrow(df)) stop("Requested n (", n, ") exceeds nrow(df) (", nrow(df), ").")
  df[sample.int(nrow(df), n), , drop = FALSE]
}

make_length_group <- function(df, group_name) {
  # Your original dt* objects needed column swap dt[, c(2,1)]
  stopifnot(is.data.frame(df), ncol(df) >= 2)
  out <- df[, c(2, 1), drop = FALSE]
  names(out) <- c("Length", "group")
  out$group <- group_name
  out
}

## 3) Subsample and harmonize data ---------------------------
n_target <- nrow(density_ref)

density_urban_sub <- sample_rows(density_urban, n_target)

dt1_sub <- make_length_group(sample_rows(dt1, n_target), "MGV")
dt2_sub <- make_length_group(sample_rows(dt2, n_target), "GOV2")
dt3_sub <- make_length_group(sample_rows(dt3, n_target), "IMGVR4")

## Ensure reference has the expected columns
stopifnot(all(c("Length", "group") %in% names(density_ref)))
stopifnot(all(c("Length", "group") %in% names(density_urban_sub)))

## 4) Statistical tests --------------------------------------
# Wilcoxon test: Reference vs Urban environment
wilcox_res <- wilcox.test(density_urban_sub$Length, density_ref$Length)
print(wilcox_res)

## 5) Combine and convert units ------------------------------
density_all <- rbind(density_ref, density_urban_sub, dt1_sub, dt2_sub, dt3_sub)
density_all <- as.data.frame(density_all)

# Convert to kb (assuming original Length is in bp)
density_all$Length <- density_all$Length / 1000

## Split by threshold (<=200 kb vs >200 kb)
density_all_200  <- subset(density_all, Length <= 200)
density_all_1000 <- subset(density_all, Length > 200)

## 6) Colors --------------------------------------------------
custom_colors <- c(
  "Reference"         = "#514e4e",
  "Urban environment" = "#cc0000",
  "MGV"               = "#3b6ea8",
  "GOV2"              = "#2f8f9d",
  "IMGVR4"            = "#e69f00"
)

## 7) Plot A: ECDF (<=200 kb) --------------------------------
library(ggplot2)

pA_ecdf <- ggplot(density_all_200, aes(x = Length, color = group)) +
  stat_ecdf(linewidth = 1.2) +
  geom_vline(xintercept = 200, linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(5, 200, 50)) +
  labs(
    x = "Genome size (kb)",
    y = "Cumulative fraction"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black")
  )

print(pA_ecdf)

## 8) Dunn test (<=200 kb) -----------------------------------
library(FSA)      # alternative: dunn.test
dunn_200 <- dunnTest(Length ~ group, data = density_all_200, method = "bonferroni")
print(dunn_200)

## 9) Summary stats (>200 kb) --------------------------------
library(dplyr)

summary_df <- density_all_1000 %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_len = median(Length, na.rm = TRUE),
    .groups = "drop"
  )

# Scale factor for overlaying median points on count bars
scale_factor <- max(summary_df$n) / max(summary_df$median_len)

## 10) Dunn test (>200 kb) -----------------------------------
dunn_1000 <- dunnTest(Length ~ group, data = density_all_1000, method = "bonferroni")
print(dunn_1000)

## 11) Plot inset: counts + median (>200 kb) ------------------
library(ggsignif)

# Annotation for Reference vs Urban environment
sig_df <- data.frame(
  xmin = "Reference",
  xmax = "Urban environment",
  y_position = max(summary_df$n) * 1.08,
  annotation = "p = 2.2 × 10\u207B\u2077"  # 2.2 × 10⁻⁷
)

p_inset <- ggplot(summary_df, aes(x = group)) +
  geom_col(aes(y = n, fill = group), alpha = 0.6, width = 0.7) +
  geom_point(aes(y = median_len * scale_factor, color = group), size = 3) +
  scale_fill_manual(values = custom_colors, name = "Group") +
  scale_color_manual(values = custom_colors, name = "Group") +
  scale_y_continuous(
    name = "Number of genomes",
    sec.axis = sec_axis(~ . / scale_factor, name = "Median genome length (kb)")
  ) +
  geom_signif(
    data = sig_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotation,
        y_position = y_position),
    manual = TRUE,
    textsize = 3.5,
    tip_length = 0.01
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank()
  )

print(p_inset)

## 12) Combine using patchwork -------------------------------
library(patchwork)

combined_plot <- pA_ecdf +
  inset_element(
    p_inset,
    left = 0.43, bottom = 0.12, right = 0.93, top = 0.75,
    align_to = "full"
  )

print(combined_plot)

## 13) Save figure -------------------------------------------
ggsave(
  filename = "Fig_1_B.pdf",
  plot = combined_plot,
  width = 5.8,
  height = 4.7357,
  units = "in",
  dpi = 300
)



########################################################
#############Fig_1_C
########################################################

## =========================================================
##  Fig 1C: N Samples vs N vOTUs by ecosystem (EN, full script)
## =========================================================

rm(list = ls(all.names = TRUE))

set.seed(123)

## 1) Load packages ------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

## 2) Load data ----------------------------------------------
data_path <- "data/all_type_Curve.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: datacombine

if (!exists("datacombine")) {
  stop("Object `datacombine` was not found after loading: ", data_path)
}

## 3) Basic checks & derived values --------------------------
# Make sure required columns exist
required_cols <- c("Eco", "N_sample", "Min", "Max", "N_ANI")
missing_cols <- setdiff(required_cols, colnames(datacombine))
if (length(missing_cols) > 0) {
  stop("Missing required columns in `datacombine`: ", paste(missing_cols, collapse = ", "))
}

# Frequency table of the first column (as in your original code)
freq_df <- as.data.frame(table(datacombine[[1]]))
colnames(freq_df) <- c("Category", "Freq")

x_max <- max(freq_df$Freq, na.rm = TRUE) * 1.2

## 4) Color palette ------------------------------------------
category_colors_type <- c(
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

## 5) Theme ---------------------------------------------------
nice_theme <- theme(
  text = element_text(color = "black", size = 8),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(color = "black", size = 8),
  panel.background = element_blank(),
  plot.background = element_blank(),
  panel.grid.major = element_line(color = "gray"),
  panel.grid.minor = element_line(color = "gray"),
  panel.border = element_rect(color = "black", fill = NA)
)

## 6) Label data (right-most point per ecosystem) ------------
label_df <- datacombine %>%
  group_by(Eco) %>%
  filter(N_sample == max(N_sample, na.rm = TRUE)) %>%
  ungroup()

## 7) Plot ----------------------------------------------------
# Parameters you may want to tune in one place:
label_nudge_x <- 50
label_size    <- 5
line_width    <- 1
ribbon_alpha  <- 0.3

p <- ggplot(datacombine) +
  geom_ribbon(
    aes(x = N_sample, ymin = Min, ymax = Max, fill = Eco),
    alpha = ribbon_alpha
  ) +
  geom_line(
    aes(x = N_sample, y = N_ANI, color = Eco),
    linewidth = line_width
  ) +
  geom_text_repel(
    data = label_df,
    aes(x = N_sample, y = N_ANI, label = Eco, color = Eco),
    nudge_x = label_nudge_x,
    direction = "y",
    hjust = 0,
    segment.size = 0.3,
    size = label_size,
    show.legend = FALSE
  ) +
  labs(x = "N Samples", y = "N vOTUs") +
  nice_theme +
  scale_x_continuous(
    labels = comma,
    limits = c(1, x_max),
    expand = expansion(mult = c(0, 0.02))  # slight right padding
  ) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = category_colors_type) +
  scale_fill_manual(values = category_colors_type) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y  = element_text(size = 16)
  )

print(p)

## 8) Save figure --------------------------------------------
ggsave(
  filename = "Fig_1_C.pdf",
  plot = p,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)


########################################################
#############Fig_1_D/E
########################################################


## =========================================================
##  Fig 1D (UpSet) + Fig 1E (Match vs Unmatched fractions)
##  Clean / standardized naming (EN)
## =========================================================
rm(list = ls(all.names = TRUE))
## 0) Packages ------------------------------------------------
suppressPackageStartupMessages({
  library(UpSetR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(reshape2)
})

## 1) Load data ----------------------------------------------
data_path <- "data/genus.RData"
stopifnot(file.exists(data_path))
load(data_path)  # expected object: datacombine


## 2) Colors ---------------------------------------------------

set_colors <- c("Urban virus" = "#FF0000",  # 红色
                "IMG/VR 4.0" = "#56B4E9",  # 绿色
                "GOV2.0" = "#56B4E9",  # 蓝色
                "MGV" = "#56B4E9",  # 蓝色
                "RefSeq" = "#56B4E9")  # 紫色



upset_plot=upset(dt01,nsets = 60,sets.bar.color =set_colors,
                 order.by = 'freq',keep.order = TRUE,
                 main.bar.color = c("red",rep("black",11)),  # 设置柱子颜色
                 mainbar.y.label = "Genus Size",
                 text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5))  # 分别调整标题、刻度、标签和交集数量的文字大小)

upset_plot

pdf(file = "Fig_1_D.pdf", width = 9, height = 6)
print(upset_plot)
dev.off()

## 4) Fig 1E: Match vs Unmatched fractions --------------------
# Totals provided by you (use clear names)

dt=data.frame(type=c("MGV","IMG/VR 4.0","GOV2.0","RefSeq"),
              all=c(189680,3661704,488131,1849),
              match=0)
dt$match=data.frame(colSums(res1))[3:6,1]

dt[,'no']=dt$all-dt$match
dt$match=dt$match/dt$all
dt$no=dt$no/dt$all
dt=dt[,-2]
dt2=data.frame(t(dt))
colnames(dt2)=dt2[1,]
dt2=dt2[-1,]
dt2[,"genus"]=row.names(dt2)

data_melt <- reshape2::melt(dt2, id.vars = "genus")

colnames(data_melt)[1]='Datasource'
data_melt$value=as.numeric(data_melt$value)

stacked_plot<-ggbarplot(data_melt, x = "variable", y = "value", fill = "Datasource", 
                        color =NA, 
                        position = position_stack(), 
                        legend = "right",palette =c('#42A5F5','#999999'),ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))+   
  labs(x = "Datasource", y = "Frequency", fill = "Datasource") +
  theme_pubr()+theme(
    axis.ticks.x = element_blank()) 

stacked_plot
pdf(file = "Fig_1_E.pdf", width = 4, height = 6)
print(stacked_plot)
dev.off()
