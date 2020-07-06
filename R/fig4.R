
# code for FIGURE 4
source('TCRsetup.R')
#rm(contigs)
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('sceMNN')
sceTCELL <- sce[, which(sce$clusTCELL != -1)]

{
    entropy <- function(p, base = exp(1)){
        p <- p[p > 0]
        -sum(p * log(p, base = base))
    }
    
    tricube <- function(x){
        out <- (1-abs(x)^3)^3
        out[abs(x) >= 1] <- 0
        return(out)
    }
    
    loess.entropy <- function(x,y, span = .5, weights = NULL, base = NULL, approx_points = 100){
        n <- length(x)
        stopifnot(length(y)==n)
        if(is.null(weights)){
            weights <- rep(1/n, n)
        }else{
            stopifnot(length(weights)==n)
            weights <- weights / sum(weights)
        }
        stopifnot(span > 0, span <= 1)
        if(is.null(base)){ base <- length(unique(y))}
        if(approx_points){
            x_int <- seq(min(x), max(x), length.out = approx_points)
        }else{
            x_int <- x
        }
        k <- max(2,ceiling(n * span))
        
        knn <- FNN::get.knnx(data = matrix(x), 
                             query = matrix(x_int), 
                             k=k)$nn.index
        
        fitted <- vapply(seq_along(x_int), function(i){
            idx <- knn[i,]
            dist <- x_int[i] - x[idx]
            dist <- dist/max(abs(dist))
            locwts <- tricube(dist)
            y.sub <- y[idx]
            tab <- unlist(by(locwts, y.sub, sum))
            entropy(tab/sum(tab), base = base)
        }, 0.0)
        
        fitted <- approx(x = x_int, y = fitted, xout = x)$y
        
        return(list(x=x, y=y, fitted=fitted))
    }
    
    # custom pie charts
    mypie <- function(countscol, ...){
        if(all(countscol==0)){
            plot(0,0,col=NA,xlab='',ylab='',axes=FALSE,
                 ...)
            return(invisible(NULL))
        }
        piecounts1 <- countscol
        piecounts1 <- piecounts1[piecounts1!=0]
        singletons <- sum(piecounts1==1)
        doublets <- sum(piecounts1==2)
        piecounts1 <- piecounts1[! piecounts1 %in% 1:2]
        clr <- hcl.colors(length(unique(piecounts1))+2)
        clr <- rev(rep(clr, times = c(1,1,table(piecounts1))))
        piecounts1 <- c(piecounts1, doublets, singletons)
        
        pie(piecounts1, labels = '', col = clr, ...)
    }
    
    # same color scheme and data, but in stacked bar plot form
    mybar <- function(countscol, ...){
        piecounts1 <- countscol
        piecounts1 <- piecounts1[piecounts1!=0]
        singletons <- sum(piecounts1==1)
        doublets <- sum(piecounts1==2)
        piecounts1 <- piecounts1[! piecounts1 %in% 1:2]
        clr <- hcl.colors(length(unique(piecounts1))+2)
        clr <- rev(rep(clr, times = c(1,1,table(piecounts1))))
        piecounts1 <- c(piecounts1, doublets, singletons)
        
        barplot(matrix(piecounts1), col = clr, border = NA, ...)
    }
    
    
} # helper functions


#######
# 4 A #
#######
####----------------------------------------------------------------------------
# overall tSNE showing TCR data, T cell clusters
tsne <- readRDS('sceMNN_tsne.rds')
shuffle <- sample(nrow(tsne))
red <- brewer.pal(9,'Set1')[1]
blue <- brewer.pal(9,'Set1')[2]
purple <- brewer.pal(11,'Set3')[10]
        
#png(filename = '~/Desktop/tSNE_TcellsTCRs.png', width=5, height=5, units = 'in', res = 300)
pdf(file='~/Desktop/tSNE_TcellsTCRs_onlyTcells.pdf', width=5,height = 5)
clr <- rep('grey75', ncol(sce))
clr[which(sce$clusTCELL != -1)] <- blue
clr[which(colnames(sce) %in% paired$barcode)] <- red
#clr[which(colnames(sce) %in% paired$barcode & sce$clusTCELL != -1)] <- purple
plot(tsne[shuffle,], col = alpha(clr[shuffle], alpha=.3), pch = 16, cex = .4, 
     asp=1, xlab='Dim 1', ylab='Dim 2', las=1,
     main='tSNE - All Cells')
dev.off()

#png(filename = '~/Desktop/tSNE_TcellsTCRs_legend.png', width=4, height=4, units = 'in', res = 300)
pdf(file='~/Desktop/tSNE_TcellsTCRs_legend.pdf', width=5,height = 5)
plot(c(0,0,1,1), c(0,1,0,1), xlim=c(-1,2), ylim=c(-1,2), asp=1, axes=FALSE, xlab='',ylab='', pch=16,
     cex = .00025 * c(sum(! colnames(sce) %in% paired$barcode & sce$clusTCELL == -1),
                     sum(colnames(sce) %in% paired$barcode & sce$clusTCELL == -1),
                     sum(! colnames(sce) %in% paired$barcode & sce$clusTCELL != -1),
                     sum(colnames(sce) %in% paired$barcode & sce$clusTCELL != -1)),
     col = c('grey75',red,blue,purple))
axis(1, at=0:1, labels = c('Other Cells','T Cells'), las=2)
axis(2, at=0:1, labels = c('No TCR','TCR'), las=2)
text(c(0,0,1,1), c(0,1,0,1), cex=.75, pos = 3,
     labels = format(c(sum(! colnames(sce) %in% paired$barcode & sce$clusTCELL == -1),
                sum(colnames(sce) %in% paired$barcode & sce$clusTCELL == -1),
                sum(! colnames(sce) %in% paired$barcode & sce$clusTCELL != -1),
                sum(colnames(sce) %in% paired$barcode & sce$clusTCELL != -1)), big.mark = ','))
dev.off()

# window <- par("usr")
# rect(window[2]-.02*(window[2]-window[1]),
#      2*window[3]/3+window[4]/3,
#      window[2],
#      window[3]/3+2*window[4]/3, col='grey75')
# rect(window[2]-.02*(window[2]-window[1]),
#      2*window[3]/3+window[4]/3,
#      window[2],
#      (2*window[3]/3+window[4]/3)+(ncol(sce)/nrow(tsne))*(window[4]-window[3])/3,
#      col='grey75')



#######
# 4 B #
#######
####----------------------------------------------------------------------------
tcrcounts <- lapply(sample_ids,function(n){
    sort(as.numeric(table(paired$clonotype_tcr[which(paired$sample==n)])),decreasing = TRUE)
})
m <- max(sapply(tcrcounts,length))
tcrcounts <- lapply(tcrcounts, function(x){
    x <- c(x,rep(0,m-length(x)))
})
names(tcrcounts) <- sample_ids
tcrcounts <- as.matrix(data.frame(tcrcounts))

pcts <- t(t(tcrcounts) / colSums(tcrcounts))


samp_plot_ord <- c("S1_N","S1_T","S2_N","S2_T","S6_N","S6_T","S8_N","S8_T",
                   "S3_N","S3_T","S7_N","S7_T","S12_N","S12_T","S14_N","S14_T",
                   "S5_N","S5_T","S11_N","S11_T","S15_N","S15_T","S16_N","S16_T")
tncol <- rep(brewer.pal(9,'Set1')[4], length(samp_plot_ord))
tncol[grep('_N',samp_plot_ord)] <- brewer.pal(9,'Set1')[2]
tncol[grep('_T',samp_plot_ord)] <- brewer.pal(9,'Set1')[1]

tcrcounts_plot <- tcrcounts[, samp_plot_ord]

#png(filename = '~/Desktop/clono_pie_bar.png', width=50, height=10, units = 'in', bg='white', res=300)
pdf(file = '~/Desktop/clono_pie_bar.pdf', width=45, height=10)
layout(matrix(1:56, nrow=2))
par(mar=c(0.1,3,3,0.1))
mypie(tcrcounts_plot[,1], main = colnames(tcrcounts_plot)[1], radius=1, cex.main=4, col.main = tncol[1])
par(mar=c(0.5,3.5,0.5,0.1))
mybar(tcrcounts_plot[,1], ylim = c(0,max(colSums(tcrcounts_plot))))
for(ii in 2:ncol(tcrcounts_plot)){
    par(mar=c(0.1,3,3,0.1))
    mypie(tcrcounts_plot[,ii], main = colnames(tcrcounts_plot)[ii], radius=1, cex.main=4, col.main = tncol[ii])
    par(mar=c(0.5,3.5,0.5,0.1))
    mybar(tcrcounts_plot[,ii], ylim = c(0,max(colSums(tcrcounts_plot))), axes=FALSE)
}
dev.off()
par(mar=c(5,4,4,2)+.1)
layout(1)

# png(filename = '~/Desktop/p1.png', width=50, height=2, units = 'in', bg='white', res=300)
# layout(matrix(1:28, nrow=1))
# for(ii in 1:ncol(tcrcounts)){
#     par(mar=c(0.1,3,3,0.1))
#     mypie(tcrcounts[,ii], main = colnames(tcrcounts)[ii], radius=1, cex.main=3.5, col.main = tncol[ii])
# }
# dev.off()
# png(filename = '~/Desktop/p2.png', width=50, height=4, units = 'in', bg='white', res=300)
# layout(matrix(1:28, nrow=1))
# for(ii in 1:ncol(tcrcounts)){
#     par(mar=c(0.5,3.5,0.5,0.5))
#     mybar(tcrcounts[,ii], ylim = c(0,max(colSums(tcrcounts))), axes=FALSE)
#     abline(h=(0:5)*2000, col='grey75',lty=2)
#     mybar(tcrcounts[,ii], add=TRUE, axes=FALSE)
# }
# dev.off()
# png(filename = '~/Desktop/p3.png', width=2, height=4, units = 'in', bg='white', res=300)
# layout(1)
# par(mar=c(0.5,9.5,0.5,0.5))
# plot(c(0,0), c(0,max(colSums(tcrcounts))), col='transparent', axes=FALSE, xlab='',ylab='')
# axis(2, cex.axis = 3, las=1)
# dev.off()









#######
# 4 D #
#######
####----------------------------------------------------------------------------
sce <- loadHDF5SummarizedExperiment('Data/CD8trajectory/')
sds <- readRDS('Data/CD8trajectory_sds.rds')


ind <- which(sce$slingCurveWeights_1 > .5)
pst <- sce$slingPseudotime_1[ind]
tcr <- sce$clonotype[ind]
samp <- sce$sample[ind]
le1 <- loess.entropy(pst,tcr, span = .5)
ind <- which(sce$slingCurveWeights_2 > .5)
pst <- sce$slingPseudotime_2[ind]
tcr <- sce$clonotype[ind]
samp <- sce$sample[ind]
le2 <- loess.entropy(pst,tcr, span = .5)

#png(filename = '~/Desktop/loessDiv_lin1.png', width=6, height=6, units = 'in', res = 200)
pdf(file = '~/Desktop/loessDiv_lower.pdf', width=6, height=6)
plot(umap[,c(1,2)], asp=1, col = 'grey80', xlab='UMAP-1', ylab='UMAP-2', pch=16, main='Lower Lineage')
points(umap[which(sce$slingCurveWeights_1 > .5),], col=colorby(le1$x), pch=16)
lines(sds, dims = c(1,2))
#
li <- as.raster(matrix(colorby(seq(0,1,length.out=50)), nrow=1))
rasterImage(li, -5.5, 2.7, -2.4, 3)
lines(c(-5.5,-5.5),c(3,5)); text(c(-5.5,-5.5),c(3,5),pos=2, labels=c(0,1), cex=.75)
text(-6.3,4,'entropy',srt=90,cex=.75)
lines(scaleAB(sort(le1$x), -5.5,-2.4), 
      scaleAB(c(0,1,le1$fitted[order(le1$x)]), 3,5)[-c(1:2)])
dev.off()

#png(filename = '~/Desktop/loessDiv_lin2.png', width=6, height=6, units = 'in', res = 200)
pdf(file = '~/Desktop/loessDiv_upper.pdf', width=6, height=6)
plot(umap[,c(1,2)], asp=1, col = 'grey80', xlab='UMAP-1', ylab='UMAP-2', pch=16, main='Upper Lineage')
points(umap[which(sce$slingCurveWeights_2 > .5),], col=colorby(le2$x), pch=16)
lines(sds, dims = c(1,2))
#
li <- as.raster(matrix(colorby(seq(0,1,length.out=50)), nrow=1))
rasterImage(li, -5.5, 2.7, -2.4, 3)
lines(c(-5.5,-5.5),c(3,5)); text(c(-5.5,-5.5),c(3,5),pos=2, labels=c(0,1), cex=.75)
text(-6.3,4,'entropy',srt=90,cex=.75)
lines(scaleAB(sort(le2$x), -5.5,-2.4), 
      scaleAB(c(0,1,le2$fitted[order(le2$x)]), 3,5)[-c(1:2)])
dev.off()


#######
# 4 F #
#######
####----------------------------------------------------------------------------
# Data matrix with each clonotype,
#      which sample it came from
#      number of cells in clonotype
#      entropy (wrt clusTCELL),
#      Normalize Shannon Index, and 
#      proportion in each cluster (for figure 4F)

clono <- data.frame(clonotype = unique(paired$clonotype))

clonoStats <- t(sapply(clono$clonotype, function(ct){
    ind <- which(paired$clonotype==ct)
    barcodes <- paired$barcode[ind]
    ncellrna <- sum(barcodes %in% colnames(sce))
    ncellTcell <- sum(barcodes %in% colnames(sceTCELL))
    if(length(ind)==1){
        majSamp <- paired$sample[ind]
        ncell <- 1
        PmajSamp <- 1
        clusEnt <- 0
        clusNSE <- 0
        clusPct <- rep(0, 19)
        if(ncellTcell){
            clusPct[sce[,barcodes]$clusTCELL+1] <- 1
        }
        names(clusPct) <- 0:18
    }else{
        tab <- table(paired$sample[ind])
        majSamp <- names(tab)[which.max(tab)]
        ncell <- sum(tab)
        PmajSamp <- max(tab)/ncell
        clus <- sceTCELL$clusTCELL[which(colnames(sceTCELL) %in% barcodes)]
        if(ncellTcell > 0){
            clusEnt <- entropy(table(clus)/length(clus))
            clusNSE <- entropy(table(clus)/length(clus), base = length(unique(clus)))
            clusPct <- table(factor(clus, levels=0:18)) / length(clus)
        }else{
            clusEnt <- 0
            clusNSE <- 0
            clusPct <- rep(0, 19)
            names(clusPct) <- 0:18
        }
    }
    return(c(majoritySample = majSamp, PCTmajoritySample = PmajSamp,
             nCells = ncell, nCellsRNA = ncellrna, nTcells = ncellTcell,
             entropy = clusEnt, normShannon = clusNSE, 
             clusPct))
}))
clono <- data.frame(clono, clonoStats)
clono$normShannon[is.nan(clono$normShannon)] <- 0


paired$sample <- factor(paired$sample)

clono <- data.frame(clonotype = unique(paired$clonotype))

clonoStatsSample <- t(sapply(clono$clonotype, function(ct){
    ind <- which(paired$clonotype==ct)
    barcodes <- paired$barcode[ind]
    ncellrna <- sum(barcodes %in% colnames(sce))
    ncellTcell <- sum(barcodes %in% colnames(sceTCELL))
    
    ind2 <- which((paired$clonotype == ct) & 
                      (paired$barcode %in% colnames(sceTCELL)))

    tab <- table(paired$sample[ind2])
    ncell <- sum(tab)
    sampPct <- tab/ncell
    
    return(c(nCells = length(ind), nCellsRNA = ncellrna, nTcells = ncellTcell,
             sampPct))
}))
clono <- data.frame(clono, clonoStatsSample)


#######
# 4 ? # Supplement
#######
####----------------------------------------------------------------------------
# TCR clonotype scatter plot of tumor vs. normal for each patient (e.g. for
# patient S1, proportion of clonotype in tumor on the y-axis, and proportion in
# normal on the x-axis)
paired <- readRDS('Data/pairedTCR_ALL.rds')

patient <- 'S1'

samps <- paste0(patient, c('_N','_T'))
# if(all(samps %in% unique(paired$sample)))

temp <- paired[paired$sample %in% samps, ]
tab <- table(temp$clonotype, temp$sample)
tab <- t(t(tab) / colSums(tab))
tab <- as.matrix(tab)

plot(log(tab[,1]), log(tab[,2]), col=rgb(0,0,0,.5))


eps <- .01
x <- log(tab[,1]+eps)
y <- log(tab[,2]+eps)
minval <- log(min(tab[tab!=0])+eps)
x[tab[,1]==0] <- runif(sum(tab[,1]==0), minval-.7, minval-.3)
y[tab[,2]==0] <- runif(sum(tab[,2]==0), minval-.7, minval-.3)
col <- sqrt(tab) / rowSums(sqrt(tab))


pdf(file=paste0('~/Desktop/clonoCor_',patient,'.pdf'), width=4,height = 4)

plot(x,y, asp=1, axes=FALSE, main=patient, xlab='Normal',ylab='Tumor',
     col = rgb(col[,2],0,col[,1], alpha = .5), pch=16)
abline(h = minval - .2, v = minval - .2)
lines(c(minval-.2,100), c(minval-.2,100))
points(x,y, col = rgb(col[,2],0,col[,1]), pch=1)
axis(1, at = c(minval-.5,log(c(.001,.01,.1,1)+eps)),
     labels = c(0,.001,.01,.1,1))
axis(2, at = c(minval-.5,log(c(.001,.01,.1,1)+eps)),
     labels = c(0,.001,.01,.1,1))

dev.off()




for(patient in paste0('S',1:16)){
    samps <- paste0(patient, c('_N','_T'))
    if(all(samps %in% unique(paired$sample))){
        
        
        temp <- paired[paired$sample %in% samps, ]
        tab <- table(temp$clonotype, temp$sample)
        tab <- t(t(tab) / colSums(tab))
        tab <- as.matrix(tab)
        
        plot(log(tab[,1]), log(tab[,2]), col=rgb(0,0,0,.5))
        
        
        eps <- .01
        x <- log(tab[,1]+eps)
        y <- log(tab[,2]+eps)
        minval <- log(min(tab[tab!=0])+eps)
        maxval <- log(max(tab[tab!=0])+eps)
        x[tab[,1]==0] <- runif(sum(tab[,1]==0), minval-.7, minval-.3)
        y[tab[,2]==0] <- runif(sum(tab[,2]==0), minval-.7, minval-.3)
        col <- sqrt(tab) / rowSums(sqrt(tab))
        
        
        pdf(file=paste0('~/Desktop/clonoCor_',patient,'.pdf'), width=5,height = 5)
        par(mar=c(4,4,4,4))
        plot(c(minval-.7,maxval),c(minval-.7,maxval), asp=1, axes=FALSE, main=patient, xlab='Normal',ylab='Tumor',
             col = 'white')
        abline(h = minval - .2, v = minval - .2)
        lines(c(minval-.2,100), c(minval-.2,100))
        points(x,y, asp=1, col = rgb(col[,2],0,col[,1], alpha = .5), pch=16)
        points(x,y, col = rgb(col[,2],0,col[,1]), pch=1)
        axis(1, at = c(minval-.5,log(c(.001,.01,.1,1)+eps)),
             labels = c(0,.001,.01,.1,1))
        axis(2, at = c(minval-.5,log(c(.001,.01,.1,1)+eps)),
             labels = c(0,.001,.01,.1,1))
        
        dev.off()
        par(mar=c(5,4,4,2)+.1)
    }
}



# tables of proportions in each part of the plot
for(patient in paste0('S',1:16)){
    samps <- paste0(patient, c('_N','_T'))
    if(all(samps %in% unique(paired$sample))){
        
        temp <- paired[paired$sample %in% samps, ]
        tab <- table(temp$clonotype, temp$sample)
        tab <- as.matrix(tab)
        
        
        
    }
}

