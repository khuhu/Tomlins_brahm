#install.packages("ggplotify")
library(grid)
library(ggplotify)
library(pheatmap)
library(ggplot2)
load("/home/kevhu/data/20181218heatmapIMG.RData")


colors.breaksNorm <- seq(-2,2,4/1000)

a <- pheatmap((heatmap.tcga), color = heatMapCol,
         cluster_cols = FALSE, cluster_rows = TRUE,fontsize = 7,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         breaks = colors.breaksNorm, border_color = NA,
         cellwidth = 5, cellheight = .5, treeheight_row = 10,
         treeheight_col = 10, show_rownames = FALSE)

grobHeat <- as.grob(a)
grid.newpage()
grid.draw(grobHeat)


grobHeat.ggplot <- as.ggplot(a)
grobHeat.ggplot

ggplot_build(grobHeat.ggplot)[[2]][1]
