####################################################################
##########################################################################

logFC_threshold=3
FDR_threshold=0.001

dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE'
setwd(dir)



files=list.files(dir,pattern='.csv')
files=c(files[grepl('Sassari.control',files)],
        files[grepl('Teramo.control',files)],
        files[grepl('BTV8',files)],
        files[grepl('Sassari.2013',files)],
        files[grepl('Teramo.2013',files)],
        files[grepl('Sassari.2006',files)],
        files[grepl('Teramo.2006',files)]
) 
gene_list=c()
for(i in files[!grepl('control',files)]){
  dat=read.csv(i,stringsAsFactors = F)
  gene_list=c(gene_list,dat$ensembl_ID[abs(as.numeric(dat[,4]))>logFC_threshold&
                                         as.numeric(dat[,3])<FDR_threshold])
}
gene_list=unique(gene_list)
length(gene_list)
gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

location1='Sassari'
location2='Teramo'
dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all'
setwd(dir)
###################################
dpis_files=sapply(files,function(x) strsplit(strsplit(x,'.vs.')[[1]][2],'.csv')[[1]][1])
#files=files[c(which(dpis_files=='1'),which(dpis_files=='3'),which(dpis_files=='7'))]
logFCs=FDRs=matrix(ncol=length(files),nrow=length(unique(gene_list)))
for(i in 1:length(files)) {
  dat=read.csv(files[i],stringsAsFactors = F)
  logFCs[,i]=dat[match(gene_list,dat$ensembl_ID),4]
  FDRs[,i]=dat[match(gene_list,dat$ensembl_ID),3]
}


colnames(logFCs)=sapply(files,function(x) gsub('.dpi','\ndpi',strsplit(x,'.csv')[[1]][1]))
# 
colnames(logFCs)=sapply(files,function(x) gsub('.dpi','',strsplit(x,'.csv')[[1]][1]))
colnames(logFCs)=colnames(FDRs)=sapply(colnames(logFCs),function(x) paste(substr((strsplit(x,'\\.')[[1]]),1,4),collapse='.'))



 
orths=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
rownames(logFCs)=rownames(FDRs)=orths$gene_names_new[match(gene_list,orths$ens_IDs)]
library(gplots)
library(pheatmap)


write.csv(FDRs[(hclust(dist(logFCs))$order),],file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/Heatmap/FDR_table.Massimo_FDR',FDR_threshold,'_logFC',logFC_threshold,'.csv',sep=''))
write.csv(logFCs[(hclust(dist(logFCs))$order),],file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/Heatmap/logFC_table.Vanessa_FDR',FDR_threshold,'_logFC',logFC_threshold,'.csv',sep=''))


###################################

rownames(logFCs)[which(gene_list=='ENSOARG00000006727')]
which(gene_list=='ENSOARG00000006727')
# pheatmap(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],
#          cluster_rows=T,cluster_cols=F,treeheight_row = 0, treeheight_col = 0)
png(filename = paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/Heatmap/logFC.heatmap.Vanessa_FDR',FDR_threshold,'_logFC',logFC_threshold,'.png',sep=''),width=360,height=round(nrow(logFCs)*8))


pheatmap(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],
         cluster_rows=T,cluster_cols=F,treeheight_row = 0, treeheight_col = 0,
         color= colorRampPalette(rev(RColorBrewer::brewer.pal(7,'RdYlBu'))[c(1:3,5:7)])(60))

length(which(rownames(logFCs)==''|is.na(rownames(logFCs))))
dev.off()
##
png(filename = paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/Heatmap/logFC.heatmap.w.tree.Vanessa_FDR',FDR_threshold,'_logFC',logFC_threshold,'.png',sep=''),width=500,height=round(nrow(logFCs)*8))
pheatmap(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],
         cluster_rows=T,cluster_cols=F,treeheight_row = 100, treeheight_col = 0,
         color= colorRampPalette(rev(RColorBrewer::brewer.pal(7,'RdYlBu'))[c(1:3,5:7)])(60))
length(which(rownames(logFCs)==''|is.na(rownames(logFCs))))
dev.off()



