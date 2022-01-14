# Title     : Heatmap generation
# Objective : To generate heatmaps comparing tools' annotation results
# Created by: iquasere
# Created on: 06/01/2022

packages <- c("DESeq2", "pheatmap", "RColorBrewer")

for (package in packages){eval(bquote(library(.(package))))}

generate.heatmap <- function(fide) {
  df <- read.table(paste('ann_paper/heatmap_df_', fide, '.tsv', sep=''), h=T, row.names=1, sep='\t')
  cs <- c('c1', 'c1', 'c1', 'c2', 'c2', 'c2')
  condition <- factor(cs)
  cd <- data.frame(cs)
  colnames(cd)[1] <- "condition"
  rownames(cd) <- colnames(df)
  dds <- DESeqDataSetFromMatrix(countData = df, colData = cd, design = ~condition)
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(cogs)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  jpeg(paste('ann_paper/', fide, '_results_differences.jpeg', sep=''))
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  dev.off()
}

generate.heatmap('COG')
generate.heatmap('EC number')
