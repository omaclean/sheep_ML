########################################
dpi='1'
strain='SMI8'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SMI8=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi',
                        dpi,'.csv',sep=''))
smi8.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SMI8dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SMI8=tru.SMI8[tru.SMI8$X%in%smi8.3prim$X,]
smi8.3prim=smi8.3prim[match(tru.SMI8$X,smi8.3prim$X),]

smi8.3prim=smi8.3prim[order(smi8.3prim$adj.P.Val,decreasing = T),]
tru.SMI8=tru.SMI8[match(smi8.3prim$X,tru.SMI8$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(smi8.3prim))
cols[smi8.3prim$adj.P.Val<0.05]=2
cols[tru.SMI8$adj.P.Val<.05]=3
cols[smi8.3prim$adj.P.Val<0.05&tru.SMI8$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'dpi',dpi,'.png',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'dpi',dpi,'.pdf',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################

########################################
dpi='3'

library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SMI8=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi',
                        dpi,'.csv',sep=''))
smi8.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SMI8dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SMI8=tru.SMI8[tru.SMI8$X%in%smi8.3prim$X,]
smi8.3prim=smi8.3prim[match(tru.SMI8$X,smi8.3prim$X),]

smi8.3prim=smi8.3prim[order(smi8.3prim$adj.P.Val,decreasing = T),]
tru.SMI8=tru.SMI8[match(smi8.3prim$X,tru.SMI8$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(smi8.3prim))
cols[smi8.3prim$adj.P.Val<0.05]=2
cols[tru.SMI8$adj.P.Val<.05]=3
cols[smi8.3prim$adj.P.Val<0.05&tru.SMI8$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.png',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.pdf',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################

########################################
dpi='1'
strain='SLI6'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SLI6=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi',
                        dpi,'.csv',sep=''))
SLI6.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SLI6dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SLI6=tru.SLI6[tru.SLI6$X%in%SLI6.3prim$X,]
SLI6.3prim=SLI6.3prim[match(tru.SLI6$X,SLI6.3prim$X),]

SLI6.3prim=SLI6.3prim[order(SLI6.3prim$adj.P.Val,decreasing = T),]
tru.SLI6=tru.SLI6[match(SLI6.3prim$X,tru.SLI6$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(SLI6.3prim))
cols[SLI6.3prim$adj.P.Val<0.05]=2
cols[tru.SLI6$adj.P.Val<.05]=3
cols[SLI6.3prim$adj.P.Val<0.05&tru.SLI6$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.png',sep=''))

plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.pdf',sep=''))

plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################

########################################
dpi='3'
strain='SLI6'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SLI6=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi',
                        dpi,'.csv',sep=''))
SLI6.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SLI6dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SLI6=tru.SLI6[tru.SLI6$X%in%SLI6.3prim$X,]
SLI6.3prim=SLI6.3prim[match(tru.SLI6$X,SLI6.3prim$X),]

SLI6.3prim=SLI6.3prim[order(SLI6.3prim$adj.P.Val,decreasing = T),]
tru.SLI6=tru.SLI6[match(SLI6.3prim$X,tru.SLI6$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(SLI6.3prim))
cols[SLI6.3prim$adj.P.Val<0.05]=2
cols[tru.SLI6$adj.P.Val<.05]=3
cols[SLI6.3prim$adj.P.Val<0.05&tru.SLI6$adj.P.Val<.05]=4
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.pdf',sep=''))
plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.png',sep=''))
plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim with controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################




########################################
dpi='1'
strain='SMI8'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SMI8=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi',
                        dpi,'.csv',sep=''))
smi8.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/no_controls_infect_dpi0_other/SMI8dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SMI8=tru.SMI8[tru.SMI8$X%in%smi8.3prim$X,]
smi8.3prim=smi8.3prim[match(tru.SMI8$X,smi8.3prim$X),]

smi8.3prim=smi8.3prim[order(smi8.3prim$adj.P.Val,decreasing = T),]
tru.SMI8=tru.SMI8[match(smi8.3prim$X,tru.SMI8$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(smi8.3prim))
cols[smi8.3prim$adj.P.Val<0.05]=2
cols[tru.SMI8$adj.P.Val<.05]=3
cols[smi8.3prim$adj.P.Val<0.05&tru.SMI8$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.no.controls.dpi',dpi,'.png',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.no.controls.dpi',dpi,'.pdf',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################

########################################
dpi='3'
strain='SMI8'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SMI8=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi',
                        dpi,'.csv',sep=''))
smi8.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/no_controls_infect_dpi0_other/SMI8dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SMI8=tru.SMI8[tru.SMI8$X%in%smi8.3prim$X,]
smi8.3prim=smi8.3prim[match(tru.SMI8$X,smi8.3prim$X),]

smi8.3prim=smi8.3prim[order(smi8.3prim$adj.P.Val,decreasing = T),]
tru.SMI8=tru.SMI8[match(smi8.3prim$X,tru.SMI8$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(smi8.3prim))
cols[smi8.3prim$adj.P.Val<0.05]=2
cols[tru.SMI8$adj.P.Val<.05]=3
cols[smi8.3prim$adj.P.Val<0.05&tru.SMI8$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.no.controls.dpi',dpi,'.png',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.no.controls.dpi',dpi,'.pdf',sep=''))

plot(smi8.3prim$logFC,tru.SMI8$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SMI8 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SMI8 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################

########################################
dpi='1'
strain='SLI6'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SLI6=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi',
                        dpi,'.csv',sep=''))
SLI6.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/no_controls_infect_dpi0_other/SLI6dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SLI6=tru.SLI6[tru.SLI6$X%in%SLI6.3prim$X,]
SLI6.3prim=SLI6.3prim[match(tru.SLI6$X,SLI6.3prim$X),]

SLI6.3prim=SLI6.3prim[order(SLI6.3prim$adj.P.Val,decreasing = T),]
tru.SLI6=tru.SLI6[match(SLI6.3prim$X,tru.SLI6$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(SLI6.3prim))
cols[SLI6.3prim$adj.P.Val<0.05]=2
cols[tru.SLI6$adj.P.Val<.05]=3
cols[SLI6.3prim$adj.P.Val<0.05&tru.SLI6$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.no.controls.dpi',dpi,'.png',sep=''))

plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.no.controls.dpi',dpi,'.pdf',sep=''))

plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################

########################################
dpi='3'
strain='SLI6'
library(RColorBrewer)
list.files('~/Documents/sheep_megadata/RNA_Seq_5.7.21/')
tru.SLI6=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi',
                        dpi,'.csv',sep=''))
SLI6.3prim=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/no_controls_infect_dpi0_other/SLI6dpi',
                          dpi,'.csv',sep=''))

colpal=c('#444444',brewer.pal(3,'Dark2'))

tru.SLI6=tru.SLI6[tru.SLI6$X%in%SLI6.3prim$X,]
SLI6.3prim=SLI6.3prim[match(tru.SLI6$X,SLI6.3prim$X),]

SLI6.3prim=SLI6.3prim[order(SLI6.3prim$adj.P.Val,decreasing = T),]
tru.SLI6=tru.SLI6[match(SLI6.3prim$X,tru.SLI6$X),]

par(mfrow=c(1,1))

cols=rep(1,nrow(SLI6.3prim))
cols[SLI6.3prim$adj.P.Val<0.05]=2
cols[tru.SLI6$adj.P.Val<.05]=3
cols[SLI6.3prim$adj.P.Val<0.05&tru.SLI6$adj.P.Val<.05]=4
png(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.no.controls.png',sep=''))

plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()

pdf(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/3prim_vs_tru/',strain,'.dpi',dpi,'.no.controls.pdf',sep=''))

plot(SLI6.3prim$logFC,tru.SLI6$logFC,pch=19,
     col=scales::alpha(colpal,.5)[cols],xlim=c(-2,4),ylim=c(-2,4),
     xlab=paste('SLI6 dpi',dpi,'_3prime sequencing logFC',sep=''),
     ylab=paste('SLI6 dpi',dpi,' tru SEQ logFC',sep=''),
     main='all genes, 3prim WITHOUT controls')

legend(legend=paste(c('FDR>0.05 in both; N=','FDR<0.05 in 3prime; N=',
                      'FDR<0.05 in tru-seq; N=','FDR<0.05 in both; N='),
                    sapply(1:4,function(x) length(which(cols==x))),sep=''),
       bty='n',pch=19,col=colpal,
       x='topleft')
abline(0,1)
dev.off()
########################################
########################################






