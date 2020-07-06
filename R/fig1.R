
# code for FIGURE 1
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('sceMNN')
tsne <- readRDS('sceMNN_tsne.rds')
require(Seurat)
ccrcc <- readRDS('Data/ccrcc.rds')
counts <- ccrcc@assays$RNA@counts
rm(ccrcc)

#######
# 1 C #
#######
####----------------------------------------------------------------------------
shuffle <- sample(nrow(tsne))
clusmeans <- t(sapply(sort(unique(sce$seurat_clusters)),function(cn){
    colMeans(tsne[which(sce$seurat_clusters==cn), ])
}))

#png(filename = '~/Desktop/tSNE_clusterLabels.png', width=5, height=5, units = 'in', res = 300)
clr <- colorby(sce$seurat_clusters, colors = brewer.pal(11,'Spectral'))
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
points(clusmeans, pch=16, col=rgb(1,1,1,.6), cex=2)
text(clusmeans, labels=sort(unique(sce$seurat_clusters)), font = 2, cex = .7,
     col = colorby(sort(unique(sce$seurat_clusters)),
                   colors = brewer.pal(11,'Spectral')))
#dev.off()
#png(filename = '~/Desktop/FIGURES/tSNE_allCells_cluster_legend.png', width=3, height=5, units = 'in', res = 300)
lgd <- unique(cbind(sce$seurat_clusters, clr))
lgd <- lgd[order(as.numeric(lgd[,1])), ]
plot.new()
legend('left','(x,y)',bty='n',legend=lgd[1:13,1], col=lgd[1:13,2], pch=16)
legend('center','(x,y)',bty='n',legend=lgd[14:26,1], col=lgd[14:26,2], pch=16)
legend('right','(x,y)',bty='n',legend=lgd[27:39,1], col=lgd[27:39,2], pch=16)
#dev.off()

    



#######
# 1 D #
#######
####----------------------------------------------------------------------------
# tSNE colored by disease stage (normal vs. early, vs. locally advanced vs. metastatic)
#png(filename = '~/Desktop/tSNE_stage.png', width=5, height=5, units = 'in', res = 300)
clr <- c(brewer.pal(8,'Set1')[2], brewer.pal(9,'Reds')[c(4,6,8)])[factor(sce$type)]
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
#dev.off()
#png(filename = '~/Desktop/tSNE_allCells_stage_legend.png', width=3, height=5, units = 'in', res = 300)
lgd <- unique(cbind(sce$type, clr))
lgd <- lgd[order(lgd[,1]), ]
plot.new()
legend('left','(x,y)',bty='n',legend=lgd[,1], col=lgd[,2], pch=16)
#dev.off()

    
#######
# 1 E #
#######
####----------------------------------------------------------------------------
# highlight different broad classes of cell types in the tSNE

# tSNE with immune cells colored
png(filename = '~/Desktop/tSNE_allCells_immune.png', width=5, height=5, units = 'in', res = 300)
clr <- rep('grey65', ncol(sce))
clr[which(sce$clusIMM != -1)] <- brewer.pal(9,'Set1')[2]
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
dev.off()
png(filename = '~/Desktop/tSNE_allCells_immune_legend.png', width=3, height=5, units = 'in', res = 300)
plot.new()
legend('left','(x,y)',bty='n',legend=c('Immune Cells','Other'), col=c(brewer.pal(9,'Set1')[2], 'grey65'), pch=16)
dev.off()


# tSNE with T cells colored
png(filename = '~/Desktop/tSNE_allCells_tcell.png', width=5, height=5, units = 'in', res = 300)
clr <- rep('grey65', ncol(sce))
clr[which(sce$clusTCELL != -1)] <- brewer.pal(9,'Set1')[2]
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
dev.off()
png(filename = '~/Desktop/tSNE_allCells_tcell_legend.png', width=3, height=5, units = 'in', res = 300)
plot.new()
legend('left','(x,y)',bty='n',legend=c('T Cells','Other'), col=c(brewer.pal(9,'Set1')[2], 'grey65'), pch=16)
dev.off()


# tSNE with myeloid cells colored
png(filename = '~/Desktop/tSNE_allCells_myeloid.png', width=5, height=5, units = 'in', res = 300)
clr <- rep('grey65', ncol(sce))
clr[which(sce$clusMYEL != -1)] <- brewer.pal(9,'Set1')[2]
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
dev.off()
png(filename = '~/Desktop/tSNE_allCells_myeloid_legend.png', width=3, height=5, units = 'in', res = 300)
plot.new()
legend('left','(x,y)',bty='n',legend=c('Myeloid Cells','Other'), col=c(brewer.pal(9,'Set1')[2], 'grey65'), pch=16)
dev.off()


# tSNE with tumor cells colored
#png(filename = '~/Desktop/tSNE_tumor.png', width=5, height=5, units = 'in', res = 300)
pdf(file='~/Desktop/tSNE_tumor.pdf', width=5,height = 5)
clr <- rep('grey65', ncol(sce))
clr[which(sce$clusALL %in% c(8,17,20))] <- brewer.pal(9,'Set1')[2]
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
dev.off()
png(filename = '~/Desktop/tSNE_tumor_legend.png', width=3, height=5, units = 'in', res = 300)
plot.new()
legend('left','(x,y)',bty='n',legend=c('Tumor Cells','Other'), col=c(brewer.pal(9,'Set1')[2], 'grey65'), pch=16)
dev.off()


#######
# 1 ? #
#######
####----------------------------------------------------------------------------
# heatmap of cluster marker genes (just need an Rdata object for David)
genelist <- c('PTPRC', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'FOXP3', 'IL2RA', 'NCAM1', 'FCGR3A', 'NCR1', 'KLRB1', 'CD19', 'MS4A1', 'CD38', 'SDC1', 'TNFRSF17', 'ITGAM', 'ITGAX', 'CSF1R', 'CD68', 'CD163', 'THBD', 'CLEC9A', 'CLEC4C', 'TPSAB1', 'KIT', 'EPCAM', 'CA9', 'ALDOB', 'UMOD')
mat <- counts[genelist, ]

clusMeans <- sapply(unique(sce$clusALL), function(cl){
    Matrix::rowMeans(mat[, which(sce$clusALL==cl)])
})
colnames(clusMeans) <- unique(sce$clusALL)
clusMeans <- clusMeans[, order(colnames(clusMeans))]

save(clusMeans, file='~/Desktop/clusMeans_ALL.Rdata')


#######
# 1 ? # Supplement
#######
####----------------------------------------------------------------------------
# Number of cells in each cluster:
# All Cells clustering
# T cells clustering
# Myeloid clustering
tab <- data.frame(table(sce$clusALL))
colnames(tab) <- c('cluster','cells')
write.csv(tab, file='~/Desktop/ncells_ALLclus.csv', row.names = FALSE)


tab <- data.frame(table(sce$clusTCELL))
colnames(tab) <- c('cluster','cells')
write.csv(tab, file='~/Desktop/ncells_TCELLclus.csv', row.names = FALSE)

tab <- data.frame(table(sce$clusMYEL))
colnames(tab) <- c('cluster','cells')
write.csv(tab, file='~/Desktop/ncells_MYELclus.csv', row.names = FALSE)




#######
# 1 ? # Supplement
#######
####----------------------------------------------------------------------------
# stacked barplot of which samples are represented in each cluster, using correct cluster names:

tab <- as.matrix(table(sce$sample, sce$clusALL))

samp.ord <- c("S1_N","S1_T","S2_N","S2_T",
              "S3_N","S3_T","S5_N","S5_T","S6_N","S6_T","S7_N","S7_T","S8_N",
              "S8_T","S10_N","S10_T2","S11_M","S11_N","S11_T","S12_N","S12_T",
              "S12_V","S14_N","S14_T","S15_N","S15_T","S16_N","S16_T")

tab <- tab[match(samp.ord, rownames(tab)), ]


lgd <- data.frame(unique(cbind(as.character(sce$sample), sce$type)))
names(lgd) <- c('sample','type')
lgd <- lgd[match(samp.ord, lgd$sample), ]
lgd$color <- brewer.pal(8,'Set1')[2]
lgd$color[lgd$type=='T_early'] <- brewer.pal(9,'Reds')[4]
lgd$color[lgd$type=='T_locAdv'] <- brewer.pal(9,'Reds')[6]
lgd$color[lgd$type=='T_metast'] <- brewer.pal(9,'Reds')[8]



pdf(file='~/Desktop/barplot_clusALL_counts.pdf',width=10,height=5)
barplot(tab, col = lgd$color, xlab='Cluster', ylab='Count',las=2)
dev.off()

pct <- t(t(tab)/colSums(tab))

pdf(file='~/Desktop/barplot_clusALL_percent.pdf',width=10,height=5)
barplot(pct, col = lgd$color, xlab='Cluster', ylab='Percent',las=2)
dev.off()







