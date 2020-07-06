
# code for FIGURE 5
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('sceMNN')
require(Seurat)
ccrcc <- readRDS('Data/ccrcc.rds')
counts <- ccrcc@assays$RNA@counts
rm(ccrcc)
sceMYEL <- sce[, which(sce$clusMYEL != -1)]
SCE <- loadHDF5SummarizedExperiment('Data/MYELtrajectory')



#######
# 5 A #
#######
####----------------------------------------------------------------------------
tsne <- readRDS('sceMNN_tsne.rds')
shuffle <- sample(1:nrow(tsne))
clusmeans <- t(sapply(sort(unique(sce$clusMYEL))[-1],function(cn){
    colMeans(tsne[which(sce$clusMYEL==cn), ])
}))

#png(filename = '~/Desktop/tSNE_myelClus.png', width=5, height=5, units = 'in', res = 300)
pdf(file = '~/Desktop/tSNE_myelClus.pdf', width=5, height=5)
clr <- rep('grey65', ncol(sce))
clr[which(sce$clusMYEL != -1)] <- colorby(sce$clusMYEL[which(sce$clusMYEL != -1)], colors = brewer.pal(11,'Spectral'))
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
points(clusmeans, pch=16, col=rgb(1,1,1,.6), cex=2)
text(clusmeans, labels=sort(unique(sce$clusMYEL))[-1], font = 2, cex = .7,
     col = colorby(sort(unique(sce$clusMYEL))[-1],
                   colors = brewer.pal(11,'Spectral')))
dev.off()
# png(filename = '~/Desktop/tSNE_myelClus_legend.png', width=3, height=5, units = 'in', res = 300)
# lgd <- unique(cbind(sce$clusMYEL[which(sce$clusMYEL != -1)], clr[which(sce$clusMYEL != -1)]))
# lgd <- lgd[order(as.numeric(lgd[,1])), ]
# lgd <- rbind(lgd, c('Other','grey65'))
# plot.new()
# legend('left','(x,y)',bty='n',legend=lgd[1:9,1], col=lgd[1:9,2], pch=16)
# legend('center','(x,y)',bty='n',legend=lgd[10:17,1], col=lgd[10:17,2], pch=16)
# dev.off()





#######
# 5 C #
#######
####----------------------------------------------------------------------------
# gene feature plots, slingshot clusters

genelist <- c('CD14','FCGR3A','MSR1','IL1B','IL6','IL8','TNF','LGALS3','LGALS9','CD163','FOLR2')
mat <- counts[genelist, colnames(SCE)]

for(gene in genelist){
    if(gene %in% rownames(mat)){
        #png(filename = paste0('~/Desktop/MYELtraj_',gene,'.png'), width=6, height=6, units = 'in', res = 100)
        pdf(file=paste0('~/Desktop/MYELtraj_',gene,'.pdf'), width=6,height = 6)
        plot(UMAPalt2D, asp=1, pch=16, main=gene,
             xlab = 'UMAP-1', ylab='UMAP-2',
             col = colorby(log1p(mat[gene,]), 
                           colors = brewer.pal(9,'Greens')[-c(2:4)],
                           alpha = .5))
        for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
        dev.off()
    }
}


#png(filename = '~/Desktop/clusters.png', width=6, height=6, units = 'in', res = 200)
pdf(file = '~/Desktop/clusters.pdf', width=6, height=6)
plot(UMAPalt2D, asp=1, pch=16, col = colorby(as.character(SCE$clusMYEL), alpha=.5),
     xlab='UMAP-1', ylab='UMAP-3', main = "Myeloid Trajectory")
for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
legendby(as.character(SCE$clusMYEL), pos = 'topright', cex=.7)
dev.off()




#######
# 5 D #
#######
####----------------------------------------------------------------------------
# umaps by stage
# stage as feature plots on umap
pdf(file = '~/Desktop/MYELtraj_N.pdf', width=6, height=6)
plot(UMAPalt2D, asp=1, col = 'grey80', main='Normal', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(UMAPalt2D[SCE$type=='N',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Set1')[2], alpha=.5), pch=16)
for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
dev.off()

pdf(file = '~/Desktop/MYELtraj_Tearly.pdf', width=6, height=6)
plot(UMAPalt2D[,c(1,2)], asp=1, col = 'grey80', main='Tumor - Early', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(UMAPalt2D[SCE$type=='T_early',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Reds')[5], alpha=.5), pch=16)
for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
dev.off()

pdf(file = '~/Desktop/MYELtraj_TlocAdv.pdf', width=6, height=6)
plot(UMAPalt2D[,c(1,2)], asp=1, col = 'grey80', main='Tumor - Locally Advanced', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(UMAPalt2D[SCE$type=='T_locAdv',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Reds')[7], alpha=.5), pch=16)
for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
dev.off()

pdf(file = '~/Desktop/MYELtraj_Tmetast.pdf', width=6, height=6)
plot(UMAPalt2D[,c(1,2)], asp=1, col = 'grey80', main='Tumor - Metastatic', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(UMAPalt2D[SCE$type=='T_metast',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Reds')[9], alpha=.5), pch=16)
for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
dev.off()


#######
# 5 E # density by stage
#######
####----------------------------------------------------------------------------
# density plot of different stages along pseudotime
SCE <- loadHDF5SummarizedExperiment('Data/MYELtrajectory')
UMAP <- reducedDim(SCE,'slingReducedDim')
UMAPalt2D <- cbind(UMAP[,1], cos(pi/6)*UMAP[,3]-sin(pi/6)*UMAP[,2])
require(slingshot)
pst <- slingPseudotime(SCE)[,1]

pdf(file='~/Desktop/MYELtraj_density.pdf', width=5,height = 5)
plot(c(0,max(pst,na.rm=TRUE)), c(0,4), col='white', axes=FALSE, 
     xlab='Pseudotime (Right Lineage)', ylab='Density', main='Pseudotime by Type')
axis(2, at=0:3, labels = c("Normal","T_early","T_locAdv","T_metast"),las=3)
abline(h=0:3)
for(i in 1:4){
    ty <- c("N","T_early","T_locAdv","T_metast")[i]
    cc <- c(brewer.pal(8,'Set1')[2], brewer.pal(9,'Reds')[c(4,6,8)])[i]
    d <- density(pst[SCE$type==ty], na.rm=TRUE)
    xx <- c(min(d$x), d$x, max(d$x))
    yy <- 4-i+c(0, 4*d$y, 0)
    polygon(xx,yy, col = alpha(cc))
    lines(xx,yy, lwd=2)
}
dev.off()


#######
# 5 E # signature feature plots - umap
#######
####----------------------------------------------------------------------------
# signature score plots
require(readxl)
genesigs <- read_excel('GeneSets/Markers_MyeloidSignature_20200505.xlsx')
cs <- Matrix::colSums(counts[,colnames(SCE)])
sigCounts <- lapply(genesigs, function(x){
    Matrix::colSums(counts[which(rownames(counts) %in% x), colnames(SCE)])
})
sigCounts <- data.frame(sigCounts)
sigScores <- log1p(1e4 * sigCounts / cs)


for(i in 1:ncol(sigScores)){
    #png(filename = paste0('~/Desktop/myeloidTrajectoryPlots/signatureScores/',names(sigScores)[i],'.png'), width=6, height=6, units = 'in', res = 150)
    pdf(file = paste0('~/Desktop/MYELtraj_',names(sigScores)[i],'.pdf'), width=6, height=6)
    plot(UMAPalt2D, asp=1, pch=16, main = names(sigScores)[i], xlab='UMAP-1',ylab='cos(30)*UMAP-3 - sin(30)*UMAP-2',
         col = colorby(sigScores[[i]], colors = brewer.pal(9,'YlGnBu'), alpha=.5))
    for(c in slingCurves(SCE)){ lines(c$s[,1], cos(pi/6)*c$s[,3]-sin(pi/6)*c$s[,2], lwd=2) }
    # li <- as.raster(matrix(colorby(seq(0,1,length.out=50), colors = brewer.pal(9,'YlGnBu')), nrow=1))
    # rasterImage(li, 1, -4.3, 4, -4)
    # text(c(1,4), c(-4,-4), c(0, round(max(sigScores[[i]]),digits = 1)), pos=3)
    # text(2.5,-3.7, 'Scores', pos=3)
    dev.off()
}







#######
# 5 E # signature tradeSeq plots
#######
####----------------------------------------------------------------------------
require(readxl)
genesigs <- read_excel('GeneSets/Markers_MyeloidSignature_20200505.xlsx')
cs <- Matrix::colSums(counts[,colnames(SCE)])
sigCounts <- lapply(genesigs, function(x){
    Matrix::colSums(counts[which(rownames(counts) %in% x), colnames(SCE)])
})
sigCounts <- data.frame(sigCounts)

mat <- t(as.matrix(sigCounts))
require(tradeSeq)
pst <- cbind(SCE$slingPseudotime_1, SCE$slingPseudotime_2)
pst[is.na(pst)] <- mean(pst, na.rm=TRUE) # necessary for tradeSeq, shouldn't make a difference (check this)
scw <- slingCurveWeights(SCE)
tot <- SCE$nUMI

ind1 <- !is.na(SCE$slingPseudotime_1)
l1sig <- fitGAM(counts = mat[,ind1],
                pseudotime = pst[ind1,1],
                cellWeights = scw[ind1,1],
                offset = log1p(tot[ind1]),
                nknots = 5)

###
require(ggplot2)
for(sig in rownames(l1sig)){
    p <- plotSmoothers(l1sig, counts(l1sig), gene = sig)
    ggsave(filename = paste0('~/Desktop/',sig,'_tradeSeq.pdf'), plot = p, width=6, height=4, units = 'in', device = 'pdf')
}

require(ggplot2)
for(sig in rownames(l1sig)){
    #png(filename = paste0('~/Desktop/tradeSeq_lines/',sig,'.png'), width=6, height=4, units = 'in', res = 150)
    pdf(file=paste0('~/Desktop/MYELtraj_',sig,'.pdf'), width=6,height = 4)
    p <- plotSmoothers(l1sig, counts(l1sig), gene = sig)
    pb <- ggplot_build(p)
    l1 <- pb$data[[2]]
    col1 <- pb$data[[3]]$colour[1]
    dat <- pb$data
    plot(dat$x, dat$y, col=col1, pch=16, cex=.5,
         xlab = 'Pseudotime', ylab = 'Avg Expression', main = sig)
    lines(l1, col=col1, lwd=6)
    #lines(l2, col=col2, lwd=6)
    #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
    dev.off()
}




#######
# 5 ? # cluster average expression for heatmap
#######
####----------------------------------------------------------------------------
genelist <- c('S100A8', 'S100A9', 'CD14', 'FCGR3A', 'SELL', 'CD68', 'CD163', 'CSF1R', 'CD69', 'CD80', 'CD86', 'HLA-DRA', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2', 'TNF', 'IL1B', 'IFNA1', 'IFNA2', 'IFNA4', 'IFNA5', 'IFNA6', 'IFNA7', 'IFNA8', 'IFNA10', 'IFNA14', 'IFNA16', 'IFNA17', 'IFNA21', 'IFNB1', 'IFNB3', 'IFNG', 'IL6', 'IL8', 'APOC1', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'FLT3', 'THBD', 'BATF', 'CLEC9A', 'CD1C', 'CLEC4C', 'TPSAB1', 'KIT', 'MKI67', 'TOP2A')
mat <- counts[genelist[genelist %in% rownames(counts)], ]

clusMeans <- sapply(unique(sce$clusMYEL), function(cl){
    Matrix::rowMeans(mat[, which(sce$clusMYEL==cl)])
})
colnames(clusMeans) <- unique(sce$clusMYEL)

save(clusMeans, file='~/Desktop/clusMeans_MYEL.Rdata')





#######
# 5 ? # Supplement
#######
####----------------------------------------------------------------------------
# Monocyte/Macrophage TradeSeq plots (full plots) for:
genelist <- c('CD14', 'FCGR3A', 'MSR1', 'IL1B', 'IL6', 'IL8', 'TNF', 'LGALS3', 'LGALS9', 'CD163', 'FOLR2')
require(tradeSeq); require(slingshot)
mat <- counts[genelist, colnames(SCE)]
mat <- mat[rowSums(mat) >= 10, ]
pst <- cbind(SCE$slingPseudotime_1, SCE$slingPseudotime_2)
pst[is.na(pst)] <- mean(pst, na.rm=TRUE) # necessary for tradeSeq, shouldn't make a difference (check this)
scw <- slingCurveWeights(SCE)
tot <- colSums(counts[,colnames(SCE)])

ind1 <- !is.na(SCE$slingPseudotime_1)
l1 <- fitGAM(counts = mat[,ind1],
             pseudotime = pst[ind1,1],
             cellWeights = scw[ind1,1],
             offset = log1p(tot[ind1]),
             nknots = 6)
ind2 <- !is.na(SCE$slingPseudotime_2)
l2 <- fitGAM(counts = mat[,ind2],
             pseudotime = pst[ind2,2],
             cellWeights = scw[ind2,2],
             offset = log1p(tot[ind2]),
             nknots = 5)


for(i in 1:nrow(l1)){
    gene <- rownames(l1)[i]
    pdf(file=paste0('~/Desktop/MYELtraj_',gene,'_tradeSeq.pdf'), width=6,height = 4)
    p <- plotSmoothers(l1, counts(l1), gene = gene)
    pb <- ggplot_build(p)
    l1.i <- pb$data[[2]]
    col1 <- pb$data[[3]]$colour[1]
    plot(range(l1.i$x), range(l1.i$y), col='white',
         xlab = 'Pseudotime', ylab = 'Avg Expression', main = gene)
    lines(l1.i, col='white', lwd=5)
    lines(l1.i, col=col1, lwd=3)
    #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
    dev.off()
}


for(i in 1:nrow(l1)){
    gene <- rownames(l1)[i]
    pdf(file=paste0('~/Desktop/MYELtraj_',gene,'_tradeSeqFull.pdf'), width=6,height = 4)
    p <- plotSmoothers(l1, counts(l1), gene = gene)
    pb <- ggplot_build(p)
    l1.i <- pb$data[[2]]
    col1 <- pb$data[[3]]$colour[1]
    p12 <- pb$data[[1]]
    plot(p12$x,p12$y, col=alpha(c(col1,col2)[p12$group], alpha=.5), pch=16,
         xlab = 'Pseudotime', ylab = 'log Counts +1 ', main = gene)
    
    lines(l1.i, col='white', lwd=5)
    lines(l1.i, col=col1, lwd=3)
    #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
    dev.off()
}


#######
# 5 ? # Supplement
#######
####----------------------------------------------------------------------------
# stacked barplot of which samples are represented in each cluster, using correct cluster names:

tab <- as.matrix(table(sce$sample, sce$clusMYEL))
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



pdf(file='~/Desktop/barplot_clusMYEL_counts.pdf',width=10,height=5)
barplot(tab, col = lgd$color, xlab='Cluster', ylab='Count',las=2)
dev.off()

pct <- t(t(tab)/colSums(tab))

pdf(file='~/Desktop/barplot_clusMYEL_percent.pdf',width=10,height=5)
barplot(pct, col = lgd$color, xlab='Cluster', ylab='Percent',las=2)
dev.off()







