library(psych);library(ggplot2);library(caret);library(randomForest)
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]

for(i in 2:ncol(counts)){
  counts[,i]=counts[,i]/sum(counts[,i])
}
 
for(i in 1:nrow(counts)){
  counts[i,2:ncol(counts)]=counts[i,2:ncol(counts)]/sum(counts[i,2:ncol(counts)])
}

any(apply(counts[,2:ncol(counts)],2,function(x) sum(x))==0)

## melt BTMs
mat=matrix(ncol=2,nrow=500000)
total=1
for(i in 1:ncol(dat)){
  genes_i=unique(dat[2:nrow(dat),i])
  genes_i=as.character(genes_i[genes_i!=' '])
  len=length(genes_i)
  mat[total:(total+len-1),]=cbind(as.character(rep(names(dat)[i],len)),genes_i)
  total=total+len-1
}
mat=mat[!is.na(mat[,1]),]
mat=mat[mat[,2]%in%counts[,1],]
## get geometric means out for counts

#take only BTMs with at least 2 genes with counts
BTMs_work=unique(mat[duplicated(mat[,1]),1])
mat2=matrix(ncol=ncol(counts[2:ncol(counts)]),nrow=length(BTMs_work))
colnames(mat2)=colnames(counts[,2:ncol(counts)])
rownames(mat2)=BTMs_work
for(i in 1:length(BTMs_work)){
  for(j in 2:ncol(counts)){
    mat2[i,j-1]=geometric.mean(1+counts[counts$Geneid%in%
                                          mat[mat[,1]==BTMs_work[i],2]
                                        ,j])-1
  }
}
dim(dat)
which(is.na(mat2),arr.ind=T)
any(apply(mat2,1,function(x) sum(x))==0)

mat2[is.na(mat2[,1]),]
###########################
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
par(mfrow=c(3,2),mar=c(4,4,1,1))
mat3=mat2[rowSums(mat2>0)>1,]
pca_btm=prcomp(t(mat3),retx=T,scale.=T,center=T)

for(dpi_plot in c('1','3','7')){
  
  types_convert=list()
  types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
  types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8_dpi',dpi_plot,sep='')
  types_convert[['SMI13']]=types_convert[['SLI13']]=paste('2013_dpi',dpi_plot,sep='')
  types_convert[['SMI6']]=types_convert[['SLI6']]=paste('2006_dpi',dpi_plot,sep='')
  types2=as.character(types_convert[types])
  
  dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
  
  types2[dpis==0]='uninfected'
  
  

  to_plot=(dpis==dpi_plot)|(types2=='uninfected')
  types2_fact=factor(types2,levels=unique(unlist(types_convert)))
  
  plot(pca_btm$x[to_plot,1],pca_btm$x[to_plot,2],pch=19,
       col=(brewer.pal(length(unique(types2)),'Dark2')[as.numeric(types2_fact)])[to_plot],
       xlab=paste('PCA1 prop dev =',100*round(pca_btm$sdev[1]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA2 prop dev =',100*round(pca_btm$sdev[2]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-13,44.5),ylim=c(-13,14))
  legend(x='bottomright',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
  
  plot(pca_btm$x[to_plot,3],pca_btm$x[to_plot,4],pch=19,
       col=(brewer.pal(length(unique(types2)),'Dark2')[as.numeric(types2_fact)])[to_plot],
       xlab=paste('PCA3 prop dev =',100*round(pca_btm$sdev[3]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA4 prop dev =',100*round(pca_btm$sdev[4]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-13,14),ylim=c(-13,14))
  legend(x='bottomright',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
  print(range(pca_btm$x[,1:4]))
}

mat2_cor=cor(t(mat2))
cors=rep(NA,nrow(mat2_cor)^2)
counter=1
for(i in 1:(ncol(mat2_cor)-1)){
  for(j in (i+1):ncol(mat2_cor)){
    cors[counter]=mat2_cor[i,j]
    counter=counter+1
  }
}
cors=cors[!is.na(cors)]
par(mfrow=c(1,1),mar=c(6,4,3,2))
hist(cors,xlab='correlation between differnt BTMs across conditions',xlim=c(-1,1),
     main=paste('histogram:',round(100*length(which(abs(cors)>0.95))/length(cors),2),
                '% comparisons >95% correlated & ',round(100*length(which(abs(cors)>0.90))/length(cors),2),
                '% comparisons >90% correlated\n'
                ,round(100*sum(apply(mat2_cor,2,function(x) length(which(abs(x)>0.95))>1))/nrow(mat2_cor),2),
                '% of BTMs 95+% correlated to at least one other BTM &',
                round(100*sum(apply(mat2_cor,2,function(x) length(which(abs(x)>0.9))>1))/nrow(mat2_cor),2),
                '90%'),breaks=20)

round(100*sum(apply(mat2_cor,2,function(x) length(which(abs(x)>0.95))>1))/nrow(mat2_cor),5)
round(100*sum(apply(mat2_cor,2,function(x) length(which(abs(x)>0.9))>1))/nrow(mat2_cor),5)

rownames(mat2)[apply(mat2_cor,2,function(x) length(which(abs(x)>0.95))>1)]


mat2_cor[colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))],
         colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))]]

pheatmap(mat2_cor[colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))],
                  colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))]]
         ,cluster_rows=F,cluster_cols=F)

###########################################################
BTMs_to_plot=(names(pca_btm$rotation[,3])[head(order(abs(pca_btm$rotation[,3]),decreasing = T),20)])

mat_plot=mat2[match(BTMs_to_plot,rownames(mat2)),]

dpis=sapply(colnames(mat_plot),function(x) strsplit(x,'_')[[1]][3])
treatments=sapply(colnames(mat_plot),function(x) strsplit(x,'\\.')[[1]][1])
colnames(mat_plot)=paste(treatments,dpis,sep='_')

cond_list_order=c('SMC_0','SMC_1','SMC_3','SMC_7',
                  'SLC_0','SLC_1','SLC_3','SLC_7',
                  'SMI8_0','SMI13_0','SMI6_0',
                  'SLI13_0','SLI6_0',
                  'SMI8_1','SMI13_1','SMI6_1',
                  'SLI13_1','SLI6_1',
                  'SMI8_3','SMI13_3','SMI6_3',
                  'SLI13_3','SLI6_3',
                  'SMI8_7','SMI13_7','SMI6_7',
                  'SLI13_7','SLI6_7')
types_convert_plot=list()
types_convert_plot[['SMC']]=types_convert_plot[['SLC']]='control'
types_convert_plot[['SLI8']]=types_convert_plot[['SMI8']]=paste('BTV8',sep='')
types_convert_plot[['SMI13']]=types_convert_plot[['SLI13']]=paste('2013',sep='')
types_convert_plot[['SMI6']]=types_convert_plot[['SLI6']]=paste('2006',sep='')


plots=list()
for(i in 1:nrow(mat_plot)){
  BTM_to_plot=i
  to_plot=data.frame(vals=mat_plot[BTM_to_plot,],conds=colnames(mat_plot),treatment=treatments)
  to_plot$conds=factor(to_plot$conds,levels = cond_list_order)
  to_plot$treatment_col=factor(types_convert_plot[to_plot$treatment],levels=c('control','BTV8','2013','2006'))
  plots[[i]]=ggplot(to_plot,aes(x=conds,y=vals,fill=treatment_col))+
    geom_dotplot(binaxis='y',stackdir = 'center',alpha=0.3,dotsize = 0.9)+theme_bw()+
    ggtitle(rownames(mat_plot)[i])+
    geom_vline(xintercept=(13.5))+
    geom_vline(xintercept=(18.5))+
    geom_vline(xintercept=(23.5))
}
library(gridExtra)
png(filename=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn/PCA3_top_20_variables.png',sep=''),width=1100,height=200*length(BTMs_to_plot))

grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]]
             ,plots[[5]], plots[[6]],plots[[7]],plots[[8]],plots[[9]]
             ,plots[[10]], plots[[11]],plots[[12]],plots[[13]],plots[[14]]
             ,plots[[16]], plots[[17]],plots[[18]],plots[[19]],plots[[20]]
             ,ncol=1)
dev.off()


