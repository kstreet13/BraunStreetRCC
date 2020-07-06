
# code for FIGURE 2
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('sceMNN')
require(Seurat)
ccrcc <- readRDS('Data/ccrcc.rds')
counts <- ccrcc@assays$RNA@counts
rm(ccrcc)
sceTCELL <- sce[, which(sce$clusTCELL != -1)]

{
    save_pheatmap_png <- function(x, filename, ...) {
        png(filename, ...)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
    }
} # helper functions

#######
# 2 A #
#######
####----------------------------------------------------------------------------
# tSNE with T cell clusters
tsne <- readRDS('sceMNN_tsne.rds')
shuffle <- sample(1:nrow(tsne))
clusmeans <- t(sapply(sort(unique(sce$clusTCELL))[-1],function(cn){
    colMeans(tsne[which(sce$clusTCELL==cn), ])
}))

png(filename = '~/Desktop/tSNE_tcellClus.png', width=5, height=5, units = 'in', res = 300)
clr <- rep('grey65', ncol(sce))
clr[which(sce$clusTCELL != -1)] <- colorby(sce$clusTCELL[which(sce$clusTCELL != -1)], colors = brewer.pal(11,'Spectral'))
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
points(clusmeans, pch=16, col=rgb(1,1,1,.6), cex=2)
text(clusmeans, labels=sort(unique(sce$clusTCELL))[-1], font = 2, cex = .7,
     col = colorby(sort(unique(sce$clusTCELL))[-1],
                   colors = brewer.pal(11,'Spectral')))
dev.off()
png(filename = '~/Desktop/tSNE_tcellClus_legend.png', width=3, height=5, units = 'in', res = 300)
lgd <- unique(cbind(Tclus, clr[which(colnames(sce) %in% t_cells)]))
lgd <- lgd[order(as.numeric(lgd[,1])), ]
lgd <- rbind(lgd, c('Other','grey65'))
plot.new()
legend('left','(x,y)',bty='n',legend=lgd[1:13,1], col=lgd[1:13,2], pch=16)
legend('center','(x,y)',bty='n',legend=lgd[14:20,1], col=lgd[14:20,2], pch=16)
dev.off()



#######
# 2 B #
#######
####----------------------------------------------------------------------------
require(pheatmap)
clus.num <- c(0,1,14,5,17,18,4,8,16,2,10,15,6,7,3,11,9,13,12)
clus.name <- c("PMCH+ terminally Exhausted CD8.1","Terminally exhausted CD8.2","Terminally exhausted CD8.3","Exhausted CD8.1","Exhausted CD8.2","Exhausted CD8.3","Cytotoxic CD8.1",'Cytotoxic CD8.2',"Cytotoxic CD8.3","CD8 TRM.1","CD8 TRM.2","CD8 TRM.3","Activated CD4.1","Activated CD4.2","CD4 TCM.1","CD4 TCM.2","Treg","Proliferating T cells","Mito-high T cells")
clus.cat <- rep(c('Terminally exhausted CD8','Exhausted CD8','Cytotoxic CD8','Tissue resident memory CD8','CD4','Proliferating','Mito-high'),
                times = c(3,3,3,3,5,1,1))
clus <- data.frame(num = clus.num, name = clus.name, cat = clus.cat)
gene <- c('CD3E', 'CD4', 'CD8A',
          'SELL', 'CCR7', 'IL7R', 'FAS', 'CD27', 'CD28', 'ITGAE', 'KLRB1',
          'PDCD1', 'TIGIT', 'HAVCR2', 'LAG3', 'CTLA4', 'KLRG1', 'TNFRSF14', 'BTLA', 'CD244', 'CD160', 'VSIR',
          'CD38', 'ENTPD1', 'NT5E',
          'CD69', 'IL2RA', 'ICOS', 'TNFRSF4', 'TNFRSF9', 'HLA-DRA', 'FASLG', 'CD40LG',
          'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GZMM', 'PRF1', 'NKG7', 'GNLY', 'LAMP1',
          'IL2', 'IFNG', 'TNF', 
          'TOX', 'TCF7', 'EOMES', 'TBX21', 'PRDM1', 'ZNF683', 'FOXP3', 'GATA3', 'LEF1', 'ID2', 'ID3', 
          'MKI67', 'TOP2A',
          'NCAM1', 'NCR1', 'FCGR3A', 
          'PMCH', 'CD74', 'CD200R1', 'TOX2', 'CCL3', 'CCL4', 'CX3CR1', 'CXCL13', 'XCL1', 'XCL2', 'FOS', 'FOSB', 'JUN')
gene.cat <- rep(c('Lineage','Naive/Memory','Inhibitory','Ectonuclease','Activation','Cytotoxic','Cytokines','Transcription Factors','Cell Cycle','NK-like','Cluster Markers'),
                times = c(3,8,11,3,8,9,3,11,2,3,13))
genes <- data.frame(gene, cat = gene.cat)
genes <- genes[genes$gene %in% rownames(counts), ]

# each cell is a column
mat <- NULL
for(cn in clus.num){
    mat <- cbind(mat, counts[genes$gene, which(sce$clusTCELL==cn)])
}
# checks
all(colnames(mat) %in% colnames(sceTCELL))
all(colnames(sceTCELL) %in% colnames(mat))
sceTCELL <- sceTCELL[, colnames(mat)]


annCol <- clus.colData[1:1000, c('clusTCELL','clus.name','clus.cat')]
annRow <- data.frame(gene.cat)
gapRow <- which(gene.cat[-length(gene.cat)] != gene.cat[-1])
gapCol <- which(clus.colData$clus.cat[-nrow(clus.colData)] != clus.colData$clus.cat[-1])

pheatmap(as.matrix(log1p(mat[,1:1000])), border_color = NA,
         gaps_row = gapRow, gaps_col = gapCol,
         annotation_row = annRow, annotation_names_row = FALSE,
         annotation_col = annCol,
         color = viridis_pal()(100))

save_pheatmap_png(ph, '~/Desktop/f_2c.png', width=20, height=10, units = 'in', res = 300)


# each cluster is a column
# average expression within cluster
avMat <- sapply(clus$num, function(cn){
    rowMeans(mat[, sceTCELL$clusTCELL == cn])
})
rownames(avMat) <- 1:nrow(avMat)
colnames(avMat) <- clus$num

annCol <- clus[,c('name','cat')]
rownames(annCol) <- clus$num
annRow <- data.frame(genes$cat)
gapRow <- which(genes$cat[-length(genes$cat)] != genes$cat[-1])
gapCol <- which(clus$cat[-nrow(clus)] != clus$cat[-1])

ph <- pheatmap(as.matrix(log1p(avMat)), border_color = NA, scale = 'row',
               cluster_rows = FALSE, cluster_cols = FALSE,
               color = viridis_pal()(100),
               gaps_row = gapRow, gaps_col = gapCol,
               annotation_row = annRow, annotation_names_row = FALSE,
               labels_row = genes$gene, fontsize_row = 7,
               annotation_col = annCol)

save_pheatmap_png(ph, '~/Desktop/f_2c.png', width=15, height=10, units = 'in', res = 300)

#











genesigs <- read_excel('GeneSets/Markers_TcellSignature_20200505.xlsx')
allgenes <- do.call('c',genesigs)
allgenes <- unique(allgenes[!is.na(allgenes)])
mat <- counts[allgenes[allgenes%in%rownames(counts)], colnames(sce)]
cs <- Matrix::colSums(counts[,colnames(sce)])






#######
# 2 D #
#######
####----------------------------------------------------------------------------
# feature plots on full tSNE
require(readxl)
genes <- read_excel('FIGURES/Genes_for_tSNE.xlsx')$Gene

for(gene in genes){
    clr <- colorby(log1p(counts[gene,]), 
                   colors = brewer.pal(9,'Greens')[-c(2:4)],
                   alpha = .5)
    #png(filename = paste0('~/Desktop/tSNE_genes/',gene,'.png'), width=5, height=5, units = 'in', res = 300)
    pdf(file=paste0('~/Desktop/tSNE_',gene,'.pdf'), width=5,height = 5)
    plot(tsne[shuffle,], asp=1, pch=16, cex=.4, main=gene,
         xlab = 'tSNE-1', ylab='tSNE-2',
         col = clr[shuffle])
    dev.off()
}






#######
# 2 E #
#######
####----------------------------------------------------------------------------

#     Dendrogram
# 2) Euclidean distance between avg. log exprn for each cluster (100 cell subsets)
avgs <- sapply(sort(unique(sceTCELL$clusTCELL)), function(clID){
    ind <- which(sceTCELL$clusTCELL==clID)
    if(length(ind) > 100){ ind <- sample(ind, 100) }
    Matrix::rowMeans(log1p(counts(sceTCELL)[, ind]))
})
hc <- hclust(dist(t(avgs)))
plot(hc)






#######
# 2 ? # Supplement
#######
####----------------------------------------------------------------------------
source('~/Projects/OLD/thingsandstuff/violinplot.R')
markers <- c('CD3E', 'CD4', 'CD8A', 'IL2RA', 'FOXP3', 'ITGAE', 'ZNF683', 'SELL', 'CCR7', 'IL7R', 'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'TOX', 'TCF7', 'MKI67', 'PMCH', 'ENTPD1', 'CD69', 'KLRG1', 'CD40LG')
mat <- log1p(counts[markers, colnames(sceTCELL)])

pdf(file='~/Desktop/markers_newTcell.pdf',width=8.5,height=11)
layout(matrix(1:5))
for(i in 1:nrow(mat)){
    violinplot( by(mat[i,], sceTCELL$clusTCELL, function(x){x}) , from=0,
                col = colorby(1:lenu(sceTCELL$clusTCELL), colors=brewer.pal(11,'Spectral')),
                main = rownames(mat)[i], ylab = 'log Counts+1', xlab = 'Cluster')
}
dev.off()
layout(1)


#######
# 2 ? # Supplement
#######
####----------------------------------------------------------------------------
# average expression in each cluster (for David to make heatmap)
genelist <- c('CD3E', 'CD4', 'CD8A',
              'SELL', 'CCR7', 'IL7R', 'FAS', 'CD27', 'CD28', 'ITGAE', 'KLRB1',
              'PDCD1', 'TIGIT', 'HAVCR2', 'LAG3', 'CTLA4', 'KLRG1', 'TNFRSF14', 'BTLA', 'CD244', 'CD160', 'VSIR',
              'CD38', 'ENTPD1', 'NT5E',
              'CD69', 'IL2RA', 'ICOS', 'TNFRSF4', 'TNFRSF9', 'HLA-DRA', 'FASLG', 'CD40LG',
              'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GZMM', 'PRF1', 'NKG7', 'GNLY', 'LAMP1',
              'IL2', 'IFNG', 'TNF', 
              'TOX', 'TCF7', 'EOMES', 'TBX21', 'PRDM1', 'ZNF683', 'FOXP3', 'GATA3', 'LEF1', 'ID2', 'ID3', 
              'MKI67', 'TOP2A',
              'NCAM1', 'NCR1', 'FCGR3A', 
              'PMCH', 'CD74', 'CD200R1', 'TOX2', 'CCL3', 'CCL4', 'CX3CR1', 'CXCL13', 'XCL1', 'XCL2', 'FOS', 'FOSB', 'JUN')
mat <- counts[genelist[genelist %in% rownames(counts)], ]

clusMeans <- sapply(unique(sce$clusALL), function(cl){
    Matrix::rowMeans(mat[, which(sce$clusALL==cl)])
})
colnames(clusMeans) <- unique(sce$clusALL)
clusMeans <- clusMeans[, order(colnames(clusMeans))]

save(clusMeans, file='~/Desktop/clusMeans_ALL.Rdata')



#######
# 2 ? # Supplement
#######
####----------------------------------------------------------------------------
# stacked barplot of which samples are represented in each cluster, using correct cluster names:

tab <- as.matrix(table(sce$sample, sce$clusTCELL))
tab <- tab[, colnames(tab) != -1]

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



pdf(file='~/Desktop/barplot_clusTCELL_counts.pdf',width=10,height=5)
barplot(tab, col = lgd$color, xlab='Cluster', ylab='Count',las=2)
dev.off()

pct <- t(t(tab)/colSums(tab))

pdf(file='~/Desktop/barplot_clusTCELL_percent.pdf',width=10,height=5)
barplot(pct, col = lgd$color, xlab='Cluster', ylab='Percent',las=2)
dev.off()



