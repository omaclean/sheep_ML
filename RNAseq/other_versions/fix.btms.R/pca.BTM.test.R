library(psych);library(ggplot2);library(caret);library(randomForest);library(RColorBrewer)
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
redundant_BTMS=c("Rho.GTPase.cycle..M4.14.","mitotic.cell.division..M6.",
                 "enriched.in.NK.cells..KIR.cluster...M61.1.","enriched.in.B.cells..VI...M69.")

dat=dat[,!colnames(dat)%in%redundant_BTMS]


for(i in 2:ncol(counts)){
  counts[,i]=counts[,i]/sum(counts[,i])
}
 
# for(i in 1:nrow(counts)){
#   counts[i,2:ncol(counts)]=counts[i,2:ncol(counts)]/sum(counts[i,2:ncol(counts)])
# }

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
  types_convert[['SMI13']]=paste('Sassari_2013_dpi',dpi_plot,sep='')
    types_convert[['SLI13']]=paste('Teramo_2013_dpi',dpi_plot,sep='')
  types_convert[['SMI6']]=paste('Sassari_2006_dpi',dpi_plot,sep='')
    types_convert[['SLI6']]=paste('Teramo_2006_dpi',dpi_plot,sep='')
  types2=as.character(types_convert[types])
  
  dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
  types2[dpis==0]='uninfected'
  
  to_plot=(dpis==dpi_plot)|(types2=='uninfected')
  types2_fact=factor(types2,levels=unique(unlist(types_convert)))
  
  plot(pca_btm$x[to_plot,1],pca_btm$x[to_plot,2],pch=19,cex=1.3,
       col=alpha((brewer.pal(length(unique(types2)),'Dark2')[c(1,2,3,6,4,5)])[as.numeric(types2_fact)],.5)[to_plot],
       xlab=paste('PCA1 prop dev =',100*round(pca_btm$sdev[1]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA2 prop dev =',100*round(pca_btm$sdev[2]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-20,20),ylim=c(-30,13))
  # legend(x='bottomright',bty='n',legend=levels(types2_fact),pch=19,
  #        col=brewer.pal(length(unique(types2)),'Dark2'))
  
  plot(pca_btm$x[to_plot,3],pca_btm$x[to_plot,4],pch=19,cex=1.3,
       col=alpha((brewer.pal(length(unique(types2)),'Dark2')[c(1,2,3,6,4,5)])[as.numeric(types2_fact)],.5)[to_plot],
       xlab=paste('PCA3 prop dev =',100*round(pca_btm$sdev[3]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA4 prop dev =',100*round(pca_btm$sdev[4]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-30,19),ylim=c(-14,20))
  legend(x='topleft',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
  
  print(paste(round(range(pca_btm$x[,1]),1),round(range(pca_btm$x[,2]),1),round(range(pca_btm$x[,3]),1),
              round(range(pca_btm$x[,4]),1),sep='; ',collapse='{' ))
}


###########################
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
par(mfrow=c(3,2),mar=c(4,4,1,1))
mat3=mat2[rowSums(mat2>0)>1,]
pca_btm=prcomp(t(mat3),retx=T,scale.=T,center=T)
for(dpi_plot in c('1','3','7')){
  types_convert=list()
  types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
  types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8_dpi',dpi_plot,sep='')
  types_convert[['SMI13']]=paste('Sassari_2013_dpi',dpi_plot,sep='')
  types_convert[['SLI13']]=paste('Teramo_2013_dpi',dpi_plot,sep='')
  types_convert[['SMI6']]=paste('Sassari_2006_dpi',dpi_plot,sep='')
  types_convert[['SLI6']]=paste('Teramo_2006_dpi',dpi_plot,sep='')
  types2=as.character(types_convert[types])
  
  dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
  types2[dpis==0]='uninfected'
  
  to_plot=(dpis==dpi_plot)|(types2=='uninfected')
  types2_fact=factor(types2,levels=unique(unlist(types_convert)))
  
  plot(pca_btm$x[to_plot,5],pca_btm$x[to_plot,6],pch=19,cex=1.3,
       col=alpha((brewer.pal(length(unique(types2)),'Dark2')[c(1,2,3,6,4,5)])[as.numeric(types2_fact)],.5)[to_plot],
       xlab=paste('PCA5 prop dev =',100*round(pca_btm$sdev[5]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA6 prop dev =',100*round(pca_btm$sdev[6]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-10,10),ylim=c(-10,26))
  # legend(x='bottomright',bty='n',legend=levels(types2_fact),pch=19,
  #        col=brewer.pal(length(unique(types2)),'Dark2'))
  
  plot(pca_btm$x[to_plot,7],pca_btm$x[to_plot,8],pch=19,cex=1.3,
       col=alpha((brewer.pal(length(unique(types2)),'Dark2')[c(1,2,3,6,4,5)])[as.numeric(types2_fact)],.5)[to_plot],
       xlab=paste('PCA7 prop dev =',100*round(pca_btm$sdev[7]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA8 prop dev =',100*round(pca_btm$sdev[8]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-10,10),ylim=c(-14,10))
  legend(x='bottomleft',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
  
  print(paste(round(range(pca_btm$x[,5]),1),round(range(pca_btm$x[,6]),1),round(range(pca_btm$x[,7]),1),
              round(range(pca_btm$x[,8]),1),collapse='{' ))
}






to_plot=(dpis=='3')|(types2=='uninfected')
plot_ly(x=pca_btm$x[to_plot,3],y=pca_btm$x[to_plot,4], 
        color=types2_fact[to_plot],colors=brewer.pal(4,'Dark2'),alpha=0.85, 
        text=~paste(sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis)[to_plot])%>%
   add_markers(marker=list(size=12))%>%
    layout(xaxis=list(title=paste('PCA3 prop dev =',100*round(pca_btm$sdev[3]^2/sum(pca_btm$sdev^2),4),'%')),
          yaxis=list(title=paste('PCA4 prop dev =',100*round(pca_btm$sdev[4]^2/sum(pca_btm$sdev^2),4),'%')),
          legend=list(font=list(size=20)))
plot_ly(x=pca_btm$x[to_plot,3],y=pca_btm$x[to_plot,4], 
        color=types2_fact[to_plot],colors=brewer.pal(4,'Dark2'),alpha=0.85, 
        text=~paste(sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis)[to_plot])%>%
  add_markers(marker=list(size=12))
 

mat2_cor=cor(t(mat2))
cors=rep(NA,nrow(mat2_cor)^2)
counter=1
for(i in 1:(ncol(mat2_cor)-1)){
  for(j in (i+1):ncol(mat2_cor)){
    cors[counter]=mat2_cor[i,j]
    counter=counter+1
  }
}











#############################
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

BTMs_to_plot=c( "innate.antiviral.response..M150."                       
                , "antiviral.IFN.signature..M75."                          
                , "type.I.interferon.response..M127."                      
                , "enriched.in.activated.dendritic.cells..II...M165."      
                ,"Activated..LPS..dendritic.cell.surface.signature..S11." 
                , "myeloid..dendritic.cell.activation.via.NFkB..I...M43.0."
                ,"Ran.mediated.mitosis..M15."   )

mat2_cor[colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))],
         colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))]]

pheatmap(mat2_cor[colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))],
                  colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))]]
         ,cluster_rows=F,cluster_cols=F)




