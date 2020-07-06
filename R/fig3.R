
# code for FIGURE 3
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('Data/CD8trajectory')
umap <- reducedDim(sce, 'UMAP')
sds <- readRDS('Data/CD8trajectory_sds.rds')
require(Seurat)
ccrcc <- readRDS('Data/ccrcc.rds')
counts <- ccrcc@assays$RNA@counts
rm(ccrcc)

#######
# 3 A #
#######
####----------------------------------------------------------------------------

# stage as feature plots on umap
plot(umap[,c(1,2)], asp=1, col = 'grey80', main='Normal', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(umap[sce$type=='N',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Set1')[2], alpha=.5), pch=16)
lines(sds, dims = c(1,2))

plot(umap[,c(1,2)], asp=1, col = 'grey80', main='Tumor - Early', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(umap[sce$type=='T_early',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Reds')[5], alpha=.5), pch=16)
lines(sds, dims = c(1,2))

plot(umap[,c(1,2)], asp=1, col = 'grey80', main='Tumor - Locally Advanced', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(umap[sce$type=='T_locAdv',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Reds')[7], alpha=.5), pch=16)
lines(sds, dims = c(1,2))

plot(umap[,c(1,2)], asp=1, col = 'grey80', main='Tumor - Metastatic', xlab='UMAP-1', ylab='UMAP-2', pch=16)
points(umap[sce$type=='T_metast',c(1,2)], asp=1, col = alpha(brewer.pal(9,'Reds')[9], alpha=.5), pch=16)
lines(sds, dims = c(1,2))


# density by stage
#
pst <- sce$slingPseudotime_1
plot(c(0,max(pst,na.rm=TRUE)), c(0,4), col='white', axes=FALSE, 
     xlab='Pseudotime (Lower Lineage)', ylab='Density', main='Pseudotime by Type')
axis(2, at=0:3, labels = c("Normal","T_early","T_locAdv","T_metast"),las=3)
abline(h=0:3)
for(i in 1:4){
    ty <- c("N","T_early","T_locAdv","T_metast")[i]
    cc <- c(brewer.pal(8,'Set1')[2], brewer.pal(9,'Reds')[c(4,6,8)])[i]
    d <- density(pst[sce$type==ty], na.rm=TRUE)
    xx <- c(min(d$x), d$x, max(d$x))
    yy <- 4-i+c(0, 2.7*d$y, 0)
    polygon(xx,yy, col = alpha(cc))
    lines(xx,yy, lwd=2)
}

pst <- sce$slingPseudotime_2
plot(c(0,max(pst,na.rm=TRUE)), c(0,4), col='white', axes=FALSE, 
     xlab='Pseudotime (Upper Lineage)', ylab='Density', main='Pseudotime by Type')
axis(2, at=0:3, labels = c("Normal","T_early","T_locAdv","T_metast"),las=3)
abline(h=0:3)
for(i in 1:4){
    ty <- c("N","T_early","T_locAdv","T_metast")[i]
    cc <- c(brewer.pal(8,'Set1')[2], brewer.pal(9,'Reds')[c(4,6,8)])[i]
    d <- density(pst[sce$type==ty], na.rm=TRUE)
    xx <- c(min(d$x), d$x, max(d$x))
    yy <- 4-i+c(0, 2.7*d$y, 0)
    polygon(xx,yy, col = alpha(cc))
    lines(xx,yy, lwd=2)
}








#######
# 3 B #
#######
####----------------------------------------------------------------------------
# tradeSeq plots
require(tradeSeq)
genelist <- readxl::read_excel('Marker Genes/Markers_TcellIndividual_20200505.xlsx')$TcellGenes
genelist <- unique(genelist)
mat <- counts[genelist[genelist %in% rownames(counts)], colnames(sce)]
pst <- cbind(sce$slingPseudotime_2, sce$slingPseudotime_1)
pst[is.na(pst)] <- mean(pst, na.rm=TRUE) # necessary for tradeSeq, shouldn't make a difference (check this)
scw <- cbind(sce$slingCurveWeights_2, sce$slingCurveWeights_1)
tot <- sce$nUMI
#icMat <- evaluateK(counts = as.matrix(counts(sce)), k=3:10, nGenes = 200, pseudotime = pst, cellWeights = scw, offset = log1p(tot))
# 6 seems like a good number (maybe try up to 20 next time)
g <- fitGAM(counts = mat,
            pseudotime = pst,
            cellWeights = scw,
            offset = log1p(tot),
            nknots = 6)

# ASSOCIATION TEST
at_global <- associationTest(g, global = TRUE, lineages = FALSE)
at_lineages <- associationTest(g, global = FALSE, lineages = TRUE)
# PATTERN TEST
pt <- patternTest(g)

for(gene in genelist){
    if(gene %in% rownames(g)){
        #png(filename = paste0('~/Desktop/tradeSeq_lines/',gene,'.png'), width=6, height=4, units = 'in', res = 150)
        pdf(file=paste0('~/Desktop/CD8traj_',gene,'_tradeSeq.pdf'), width=6,height = 4)
        p <- plotSmoothers(g, counts(g), gene = gene)
        pb <- ggplot_build(p)
        l1 <- pb$data[[2]]
        col1 <- pb$data[[3]]$colour[1]
        l2 <- pb$data[[4]]  
        col2 <- '#35B779FF' # pb$data[[5]]$colour[1]
        plot(range(c(l1$x,l2$x)), range(c(l1$y,l2$y)), col='white',
             xlab = 'Pseudotime', ylab = 'Avg Expression', main = gene)
        lines(l1, col=col1, lwd=3)
        lines(l2, col=col2, lwd=3)
        #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
        dev.off()
    }
}




#######
# 3 E #
#######
####----------------------------------------------------------------------------
# gene feature plots

genelist <- c('PDCD1','TCF7','HAVCR2','TBX21','LAG3','TOX','ENTPD1','EOMES','JUN','FO')
#gene <- 'JUN'
mat <- counts[genelist, colnames(sce)]

for(gene in genelist){
    if(gene %in% rownames(mat)){
        #png(filename = paste0('~/Desktop/newCD8TrajectoryPlots/genes/',gene,'.png'), width=6, height=6, units = 'in', res = 100)
        pdf(file=paste0('~/Desktop/CD8traj_',gene,'.pdf'), width=6,height = 6)
        plot(umap, asp=1, pch=16, main=gene,
             xlab = 'UMAP-1', ylab='UMAP-2',
             col = colorby(log1p(mat[gene,]), 
                           colors = brewer.pal(9,'Greens')[-c(2:4)],
                           alpha = .5))
        lines(sds, dims=c(1,2), lwd=2)
        dev.off()
    }
}


png(filename = '~/Desktop/JUN.png', width=6, height=6, units = 'in', res = 150)
plot(umap, asp=1, pch=16, main=gene,
     xlab = 'UMAP-1', ylab='UMAP-2',
     col = colorby(log1p(mat[gene,]), 
                   colors = brewer.pal(9,'Greens')[-c(2:4)],
                   alpha = .5))
lines(sds, dims=c(1,2), lwd=2)
dev.off()









#######
# 3 F #
#######
####----------------------------------------------------------------------------
# gene set signature feature plots
require(readxl)
genesigs <- read_excel('GeneSets/Markers_TcellSignature_20200505.xlsx')
cs <- Matrix::colSums(counts[,colnames(sce)])
sigCounts <- lapply(genesigs, function(x){
    Matrix::colSums(counts[which(rownames(counts) %in% x), rownames(umap),drop=FALSE])
})
sigCounts <- data.frame(sigCounts)
sigScores <- log1p(1e4 * sigCounts / cs)

# signature scores as feature plots

for(i in 1:ncol(sigScores)){
    #png(filename = paste0('~/Desktop/newCD8TrajectoryPlots/signatures/',names(sigScores)[i],'.png'), width=6, height=6, units = 'in', res = 100)
    pdf(file=paste0('~/Desktop/CD8traj_signatures/CD8traj_',names(sigScores)[i],'.pdf'), width=6,height = 6)
    plot(umap[,c(1,2)], asp=1,  main = names(sigScores)[i], xlab='UMAP-1', ylab='UMAP-2', pch=16,
         col = colorby(sigScores[[i]], 
                       colors = brewer.pal(9,'YlGnBu'),
                       alpha = .5))
    lines(sds, dims = c(1,2))
    li <- as.raster(matrix(colorby(seq(0,1,length.out=50), colors = brewer.pal(9,'YlGnBu')), nrow=1))
    rasterImage(li, -6, 4, -4, 4.3)
    text(c(-6,-4), c(4.3,4.3), c(0, round(max(sigScores[[i]]),digits = 1)), pos=3)
    text(-5,4.5, 'Scores', pos=3)
    dev.off()
}




#######
# 3 F #
#######
####----------------------------------------------------------------------------
# gene set signature tradeSeq plots
require(readxl)
genesigs <- read_excel('GeneSets/Markers_TcellSignature_20200505.xlsx')
cs <- Matrix::colSums(counts[,colnames(sce)])
sigCounts <- lapply(genesigs, function(x){
    Matrix::colSums(counts[which(rownames(counts) %in% x), rownames(umap),drop=FALSE])
})
sigCounts <- data.frame(sigCounts)
sigCounts <- t(as.matrix(sigCounts))
require(tradeSeq)
pst <- cbind(sce$slingPseudotime_2, sce$slingPseudotime_1)
pst[is.na(pst)] <- mean(pst, na.rm=TRUE) # necessary for tradeSeq, shouldn't make a difference (check this)
scw <- cbind(sce$slingCurveWeights_2, sce$slingCurveWeights_1)
tot <- sce$nUMI
#icMat <- evaluateK(counts = as.matrix(counts(sce)), k=3:10, nGenes = 200, pseudotime = pst, cellWeights = scw, offset = log1p(tot))
# 6 seems like a good number (maybe try up to 20 next time)
g <- fitGAM(counts = sigCounts,
            pseudotime = pst,
            cellWeights = scw,
            offset = log1p(tot),
            nknots = 6)

# ASSOCIATION TEST
at_global <- associationTest(g, global = TRUE, lineages = FALSE)
at_lineages <- associationTest(g, global = FALSE, lineages = TRUE)
# PATTERN TEST
pt <- patternTest(g)
#pt <- patternTest(g, l2fc = 0.25)

# plot
sig <- 'Wherry_Exhaustion'
plotSmoothers(g, counts(g), gene = sig)

require(ggplot2)
for(sig in rownames(sigCounts)){
    #png(filename = paste0('~/Desktop/tradeSeq_lines/',sig,'.png'), width=6, height=4, units = 'in', res = 150)
    pdf(file=paste0('~/Desktop/CD8traj_signatures_tradeSeq/CD8traj_',sig,'.pdf'), width=6,height = 4)
    p <- plotSmoothers(g, counts(g), gene = sig)
    pb <- ggplot_build(p)
    l1 <- pb$data[[2]]
    col1 <- pb$data[[3]]$colour[1]
    l2 <- pb$data[[4]]
    col2 <- '#35B779FF' # pb$data[[5]]$colour[1]
    plot(range(c(l1$x,l2$x)), range(c(l1$y,l2$y)), col='white',
         xlab = 'Pseudotime', ylab = 'Avg Expression', main = sig)
    lines(l1, col=col1, lwd=6)
    lines(l2, col=col2, lwd=6)
    #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
    dev.off()
}



#######
# 3 ? # Supplement
#######
####----------------------------------------------------------------------------
# CD8 TradeSeq plots (full plots, with individual cell expression) for:
genelist <- c('PDCD1', 'HAVCR2', 'LAG3', 'ENTPD1', 'TCF7', 'TBX21', 'TOX', 'EOMES', 'JUN', 'FOS')
require(tradeSeq)
mat <- counts[genelist, colnames(sce)]
pst <- cbind(sce$slingPseudotime_2, sce$slingPseudotime_1)
pst[is.na(pst)] <- mean(pst, na.rm=TRUE) # necessary for tradeSeq, shouldn't make a difference (check this)
scw <- cbind(sce$slingCurveWeights_2, sce$slingCurveWeights_1)
tot <- sce$nUMI
#icMat <- evaluateK(counts = as.matrix(counts(sce)), k=3:10, nGenes = 200, pseudotime = pst, cellWeights = scw, offset = log1p(tot))
# 6 seems like a good number (maybe try up to 20 next time)
g <- fitGAM(counts = mat,
            pseudotime = pst,
            cellWeights = scw,
            offset = log1p(tot),
            nknots = 6)


for(gene in genelist){
    #png(filename = paste0('~/Desktop/tradeSeq_lines/',gene,'.png'), width=6, height=4, units = 'in', res = 150)
    pdf(file=paste0('~/Desktop/CD8traj_',gene,'_tradeSeqFull.pdf'), width=6,height = 4)
    p <- plotSmoothers(g, counts(g), gene = gene)
    pb <- ggplot_build(p)
    l1 <- pb$data[[2]]
    col1 <- pb$data[[3]]$colour[1]
    l2 <- pb$data[[4]]  
    col2 <- '#35B779FF' # pb$data[[5]]$colour[1]
    p12 <- pb$data[[1]]
    plot(p12$x,p12$y, col=alpha(c(col1,col2)[p12$group], alpha=.5), pch=16,
         xlab = 'Pseudotime', ylab = 'log Counts +1 ', main = gene)
    
    lines(l1, col='white', lwd=5)
    lines(l1, col=col1, lwd=3)
    lines(l2, col='white', lwd=5)
    lines(l2, col=col2, lwd=3)
    #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
    dev.off()
}

#

require(readxl)
genesigs <- read_excel('GeneSets/Markers_TcellSignature_20200505.xlsx')
siglist <- c('Nick_Inhibition', 'Azizi_Tcell-Terminal', 'SadeFeldman_Exhausted_T_cells', 'Singer_Activation_Dysfunction_module_100', 'IdoAmit_Stress')
genesigs <- genesigs[, siglist]
sigCounts <- lapply(genesigs, function(x){
    Matrix::colSums(counts[which(rownames(counts) %in% x), rownames(umap),drop=FALSE])
})
sigCounts <- data.frame(sigCounts)
sigCounts <- t(as.matrix(sigCounts))
require(tradeSeq)
pst <- cbind(sce$slingPseudotime_2, sce$slingPseudotime_1)
pst[is.na(pst)] <- mean(pst, na.rm=TRUE) # necessary for tradeSeq, shouldn't make a difference (check this)
scw <- cbind(sce$slingCurveWeights_2, sce$slingCurveWeights_1)
tot <- sce$nUMI
#icMat <- evaluateK(counts = as.matrix(counts(sce)), k=3:10, nGenes = 200, pseudotime = pst, cellWeights = scw, offset = log1p(tot))
# 6 seems like a good number (maybe try up to 20 next time)
g <- fitGAM(counts = sigCounts,
            pseudotime = pst,
            cellWeights = scw,
            offset = log1p(tot),
            nknots = 6)

for(i in 1:nrow(g)){
    sig <- rownames(g)[i]
    #png(filename = paste0('~/Desktop/tradeSeq_lines/',sig,'.png'), width=6, height=4, units = 'in', res = 150)
    pdf(file=paste0('~/Desktop/CD8traj_',sig,'_tradeSeqFull.pdf'), width=6,height = 4)
    p <- plotSmoothers(g, counts(g), gene = sig)
    pb <- ggplot_build(p)
    l1 <- pb$data[[2]]
    col1 <- pb$data[[3]]$colour[1]
    l2 <- pb$data[[4]]  
    col2 <- '#35B779FF' # pb$data[[5]]$colour[1]
    p12 <- pb$data[[1]]
    plot(p12$x,p12$y, col=alpha(c(col1,col2)[p12$group], alpha=.5), pch=16,
         xlab = 'Pseudotime', ylab = 'log Counts +1 ', main = sig)
    
    lines(l1, col='white', lwd=5)
    lines(l1, col=col1, lwd=3)
    lines(l2, col='white', lwd=5)
    lines(l2, col=col2, lwd=3)
    #legend('top','(x,y)',bty='n',legend=c('Lineage 1','Lineage 2'), col=c(col1,col2), lwd=3)
    dev.off()
}



