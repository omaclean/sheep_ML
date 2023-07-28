dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs.fromArtur.csv')

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
counts=counts[rowSums(counts[,2:ncol(counts)>0])>14,]

mat=matrix(ncol=2,nrow=500000)
total=1

##fix up BTMS- remove duplicates and melt 
for(i in 1:ncol(dat)){
  genes_i=unique(dat[2:nrow(dat),i])
  genes_i=as.character(genes_i[genes_i!=' '])
  len=length(genes_i)
  
  mat[total:(total+len-1),]=cbind(as.character(rep(names(dat)[i],len)),genes_i)
  if(any(duplicated(genes_i))){
    print(c(i,length(which(duplicated(genes_i)))))}
  total=total+len-1
}
mat=mat[!is.na(mat[,1]),]

for(i in unique(mat[,1])){
  if(any(duplicated((mat[mat[,1]==i,2])))){
    print(i)
  }
}
###########################


##########

library(clusterProfiler)

counts=counts[rowSums(counts[2:ncol(counts)]>0)>14,]

order=unique(counts$Geneid)
order=sample(order,length(order),replace=F)
logFC=rnorm(length(order))
order2=sort(logFC,decreasing = T)
names(order2)=order
test=GSEA(order2,TERM2GENE = mat)



dat=read.csv('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/Teramo.2006.dpi.0.vs.3.csv')
order=sort(dat$Teramo_SLI6_logFC,decreasing=T)
names(order)=dat$ensembl_ID
test=GSEA(order,TERM2GENE = mat)
dotplot(test)
dat=read.csv('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/Teramo.2006.dpi.0.vs.1.csv')
order=sort(dat$Teramo_SLI6_logFC,decreasing=T)
names(order)=dat$ensembl_ID
test=GSEA(order,TERM2GENE = mat)
dotplot(test)

condition='Sassari.BTV8.dpi.0.vs.1'
dat=read.csv(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',condition,'.csv',sep=''))
order=sort(dat[,4],decreasing=T)
names(order)=dat$ensembl_ID
test=GSEA(order,TERM2GENE = mat)
dotplot(test,title=condition)

#####################
mat2=mat
mat2[,1]=gsub('type.I','type-I',mat[,1])
mat2[,1]=gsub('\\.\\.','\n',mat2[,1])
mat2[,1]=gsub('\\.','\n',mat2[,1])
mat2[,1]=gsub('\n\n','',mat2[,1])
mat2[,1]=gsub('\n \n','',mat2[,1])
mat2[,1]=sapply(mat2[,1],function(x) paste(strsplit(x,'\n')[[1]][1:4],collapse='\n'))

#

conditions=c('Teramo.control.dpi.0.vs.1','Teramo.control.dpi.0.vs.3','Teramo.control.dpi.0.vs.7',
             'Sassari.control.dpi.0.vs.1','Sassari.control.dpi.0.vs.3','Sassari.control.dpi.0.vs.7')

plots=list()
for(i in 1:length(conditions)){
  condition=conditions[i]
dat=read.csv(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',condition,'.csv',sep=''))
order=sort(dat[,4],decreasing=T)
names(order)=dat$ensembl_ID
test=GSEA(order,TERM2GENE = mat2)
if(nrow(test)==0){
  plots[[i]]=ggplot(dat)+ggtitle(condition)
  }else{
  plots[[i]]=dotplot(test,title=condition)
  }
}

png(filename = paste('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/GSEA_BTM/',condition,'.png',sep=''),width=900,height=900)
plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]])
dev.off()

##############################################################
conditions=c('Teramo.2006.dpi.0.vs.1','Teramo.2013.dpi.0.vs.1',
             'Sassari.2006.dpi.0.vs.1','Sassari.2013.dpi.0.vs.1')
plots=list()
for(i in 1:length(conditions)){
  condition=conditions[i]
  dat=read.csv(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',condition,'.csv',sep=''))
  order=sort(dat[,4],decreasing=T)
  names(order)=dat$ensembl_ID
  test=GSEA(order,TERM2GENE = mat2)
  if(nrow(test)==0){
    plots[[i]]=ggplot(dat)+ggtitle(condition)
  }else{
    plots[[i]]=dotplot(test,title=condition)
  }
}
png(filename = paste('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/GSEA_BTM/',condition,'.png',sep=''),width=900,height=900)
plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]])
dev.off()

############################
conditions=c('Teramo.2006.dpi.0.vs.3','Teramo.2013.dpi.0.vs.3',
             'Sassari.2006.dpi.0.vs.3','Sassari.2013.dpi.0.vs.3')
plots=list()
for(i in 1:length(conditions)){
  condition=conditions[i]
  dat=read.csv(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',condition,'.csv',sep=''))
  order=sort(dat[,4],decreasing=T)
  names(order)=dat$ensembl_ID
  test=GSEA(order,TERM2GENE = mat2)
  if(nrow(test)==0){
    plots[[i]]=ggplot(dat)+ggtitle(condition)
  }else{
    plots[[i]]=dotplot(test,title=condition)
  }
}
png(filename = paste('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/GSEA_BTM/',condition,'.png',sep=''),width=900,height=900)
plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]])
dev.off()

############################
conditions=c('Teramo.2006.dpi.0.vs.7','Teramo.2013.dpi.0.vs.7',
             'Sassari.2006.dpi.0.vs.7','Sassari.2013.dpi.0.vs.7')
plots=list()
for(i in 1:length(conditions)){
  condition=conditions[i]
  dat=read.csv(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',condition,'.csv',sep=''))
  order=sort(dat[,4],decreasing=T)
  names(order)=dat$ensembl_ID
  test=GSEA(order,TERM2GENE = mat2)
  if(nrow(test)==0){
    plots[[i]]=ggplot(dat)+ggtitle(condition)
  }else{
    plots[[i]]=dotplot(test,title=condition)
  }
}
png(filename = paste('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/GSEA_BTM/',condition,'.png',sep=''),width=900,height=900)
plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]])
dev.off()
############################
conditions=c('Sassari.BTV8.dpi.0.vs.1','Sassari.BTV8.dpi.0.vs.3','Sassari.BTV8.dpi.0.vs.7')
plots=list()
for(i in 1:length(conditions)){
  condition=conditions[i]
  dat=read.csv(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',condition,'.csv',sep=''))
  order=sort(dat[,4],decreasing=T)
  names(order)=dat$ensembl_ID
  test=GSEA(order,TERM2GENE = mat2)
  if(nrow(test)==0){
    plots[[i]]=ggplot(dat)+ggtitle(condition)
  }else{
    plots[[i]]=dotplot(test,title=condition)
  }
}
png(filename = paste('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/GSEA_BTM/',condition,'.png',sep=''),width=900,height=900)
plot_grid(plots[[1]],plots[[2]],plots[[3]])
dev.off()



