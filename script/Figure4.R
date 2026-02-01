

########################################################
#############Fig_4_B
########################################################

rm(list = ls(all = TRUE))

library(data.table)
library(stringr)
library(tidyr)
library(ggpubr)


load('data/data_hugephage.RData')
density_ref=data_hugephage$length[which(data_hugephage$source=='urban environment')]
density_ref=data.frame(density_ref)
density_ref[,'group']='Reference'
colnames(density_ref)[1]="Length"

density_urban=data_hugephage$length[which(data_hugephage$source=='Reference')]
density_urban=data.frame(density_urban)
density_urban[,'group']='Urban environment'
colnames(density_urban)[1]="Length"

p_Length=wilcox.test(density_urban$Length,density_ref$Length)
density_all=rbind(density_ref,density_urban)
density_all$Length=density_all$Length/1000
density_all=data.frame(density_all)

comparisons <- list(c("Reference", "Urban environment"))
p_boxplot <- ggboxplot(density_all, x = "group", y = "Length",
                       color = "group", palette = c("#514e4e", "#cc0000"),
                       add = "jitter",  
                       ylab = "Length") +
  theme(legend.position = "none",  
        axis.text.x = element_blank(),  
        axis.title.x = element_blank())+  
  stat_compare_means(comparisons = comparisons,method = "wilcox.test", label = "p.signif")



total_count <- nrow(density_all)


bw_value <- density(density_all$Length)$bw

custom_colors <- c("Reference" = "#514e4e", "Urban environment" = "#cc0000")


p_density<-ggplot(density_all, aes(x = Length, color = group, fill = group)) +
  geom_density(alpha = 0.6, adjust = 0.5, aes(y = ..density.. * total_count * bw_value)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = c(seq(0, 200, by = 50),seq(200,500,by=100),seq(500,1500,by=300)))+
  labs(x = "Genome size (kb)", y = "Genomes (count)", title = "") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    legend.position = c(0.95, 0.3),             
    legend.justification = c(1, 1)     
  )
p_density

comparisons <- list(c("Reference", "urban environment"))

p_boxplot_gc <- ggboxplot(data_hugephage, x = "source", y = "GC",
                          color = "source", palette = c("#514e4e", "#cc0000"),
                          add = "jitter",  
                          ylab = "GC%",
                          xlab = "") +
  theme(legend.position = "none",  
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())+  
  stat_compare_means(comparisons = comparisons,method = "wilcox.test", label = "p.signif")+
  coord_flip()

bar_life=data.frame(table(data_hugephage[,c("lifestyle","source")]))

tab <- xtabs(Freq ~ lifestyle + source, data = bar_life)
print(tab)


chi_res <- chisq.test(tab)
pval <- chi_res$p.value
print(chi_res)

colnames(bar_life)[2]="group"
colnames(bar_life)[3]="value"
gpl=unique(bar_life$group)
for (i in gpl) {
  bar_life[which(bar_life$group==i),"value"]= bar_life[which(bar_life$group==i),"value"]/sum( bar_life[which(bar_life$group==i),"value"])
}


library(ggpubr)


height_summary <- aggregate(value ~ group, bar_life, sum)
mid_y <- mean(height_summary$value)

pval_df <- data.frame(
  group1 = "Reference",
  group2 = "urban environment",
  y.position = mid_y * 1.05,   
  p = pval,
  p.signif = ifelse(pval < 0.0001, "****",
                    ifelse(pval < 0.001, "***",
                           ifelse(pval < 0.01, "**",
                                  ifelse(pval < 0.05, "*", "ns"))))
)


category_colors <- c("Temperate" = "#1f77b4",   
                     "Virulent"  = "#ff7f0e")   

lifestyle_stacked_bar_plot <- ggbarplot(
  bar_life, x = "group", y = "value",
  fill = "lifestyle", color = "white",
  position = position_stack(reverse = TRUE)
) +
  scale_fill_manual(values = category_colors) +
  labs(x = "Group", y = "Value") +
  theme_minimal() +
  stat_pvalue_manual(
    pval_df,
    label = "p.signif",
    tip.length = 0.02,
    angle = 90    
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  coord_flip()


library(patchwork)

combined_plot <- p_density + 
  inset_element(p_boxplot_gc, left = 0.28, bottom = 0.65, right = 0.85, top = 0.95, 
                align_to = "full")+
  inset_element(lifestyle_stacked_bar_plot, left = 0.3, bottom = 0.4, right = 0.95, top = 0.65, 
                align_to = "full")

print(combined_plot)
ggsave("Fig_4_B.pdf", combined_plot, width = 8, height = 6, units = "in", dpi = 300)

