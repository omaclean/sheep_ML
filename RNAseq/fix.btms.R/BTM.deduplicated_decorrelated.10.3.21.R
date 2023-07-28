library(psych);library(ggplot2);library(caret);library(randomForest)
dpi_plot='7'

dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
#filter counts & BTMs so that gene is expressed in at least 2 conditions (on average)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
####################
#remove BTMs which are contained in another (i.e. full subsets)
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

## melt BTMs generate 2xN matrix of BTM name & genes insode
mat=matrix(ncol=2,nrow=500000)
total=1
for(i in 1:ncol(dat)){
  genes_i=unique(dat[2:nrow(dat),i])
  genes_i=as.character(genes_i[genes_i!=' '])
  len=length(genes_i)
  mat[total:(total+len-1),]=cbind(as.character(rep(names(dat)[i],len)),genes_i)
  total=total+len-1
}
#remove NAs from initiliased matrix (unused rows)
mat=mat[!is.na(mat[,1]),]
#filter down BTM describer to only includes genes which met the >14 count threshold above
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




#################

mat2_cor=cor(t(mat2))
pairs=matrix(ncol=2,nrow=0)
for(i in 1:(ncol(mat2_cor)-1)){
  count=0
  for(j in (i+1):ncol(mat2_cor)){
    if(mat2_cor[i,j]>0.9){
      pairs=rbind(pairs,c(i,j))
    }
  }
}

to_combine=list()
for(i in 1:nrow(pairs)){
    to_combine[[i]]=pairs[i,]
}
# in case ordering of above makes things not work
to_remove=c()

hits=T
while(hits==T){
  hits=F
  if(length(to_combine)>1){
    for(i in 1:(length(to_combine)-1)){
      for(j in (i+1):length(to_combine)){
        if(any(to_combine[[i]]%in%to_combine[[j]])){
          hits=T
          to_combine[[i]]=unique(c(to_combine[[i]],to_combine[[j]]))
          to_combine[[j]]=NULL
          to_remove=c(to_remove,j)
          break
        }
      }
      if(hits==T){break}
    }
  }
}
######################################################
#mat_copy=mat

for(i in 1:length(to_combine)){
  mat[mat[,1]%in%BTMs_work[to_combine[[i]]] ,1]   =paste(c('combine(',BTMs_work[to_combine[[i]]]  ,')'),collapse=' ')
}
to_del=c()
for(i in unique(mat[,1])){
  hits=which(mat[,1]==i)
  to_del=c(to_del,hits[duplicated(mat[hits,2])])
}
if(length(to_del)>0){
  mat=mat[-to_del,]
}

#remove NAs from initiliased matrix (unused rows)
mat=mat[!is.na(mat[,1]),]

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




####################



mat3=mat2

####
# convert colnames from counts into 'types' &extract DPI
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
types_convert=list()
types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8'
types_convert[['SMI13']]=types_convert[['SLI13']]='2013'
types_convert[['SMI6']]=types_convert[['SLI6']]='2006'
types2=as.character(types_convert[types])
dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
types3=paste(types2,'_dpi',dpis)

rfedat=t(mat3[,to_plot])
to_plot=(dpis==dpi_plot)|(types2=='uninfected'|dpis=='0')
BTV8DPI1=types2
BTV8DPI1[dpis=='0']='uninfected'
BTV8DPI1=as.factor(BTV8DPI1[to_plot])

############################################################
################################################################
################################################################################################




RF=randomForest(y= BTV8DPI1,x= rfedat,
                importance=T,do.trace = F)

RF2=randomForest(y= BTV8DPI1,x= rfedat[,rownames(RF$importance)[order(RF$importance[,4],decreasing = T)[1:20]]],
                 importance=T,do.trace = F)
RF2$confusion


to_check=BTV8DPI1!='BTV8'

RF3=randomForest(y= as.factor(as.character(BTV8DPI1[to_check])),x= rfedat[to_check,],
                importance=T,do.trace = F)

RF4=randomForest(y= as.factor(as.character(BTV8DPI1[to_check])),x= rfedat[to_check,rownames(RF$importance)[order(RF$importance[,4],decreasing = T)[1:20]]],
                 importance=T,do.trace = F)
# 
# ##############################
# infected=!(types2=='uninfected'|dpis=='0')[to_plot]
# excess_var=rep(0,ncol(rfedat))
# for(i in 1:ncol(rfedat)){
#   if(var(rfedat[!infected,i])>var(rfedat[,i])){
#     excess_var[i]=1
#   }
# }
# sum(excess_var)/length(excess_var)
# rfedat=rfedat[,!excess_var]
# 
# 
# RF5=randomForest(y= BTV8DPI1,x= rfedat,
#                  importance=T,do.trace = F)
# 
# RF6=randomForest(y= BTV8DPI1,x= rfedat[,rownames(RF5$importance)[order(RF5$importance[,4],decreasing = T)[1:20]]],
#                  importance=T,do.trace = F)


RF$confusion
RF2$confusion
RF3$confusion
RF4$confusion
RF5$confusion
RF6$confusion

#######################
par(mfrow=c(3,2))
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
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-16,18),ylim=c(-15,32))
  # legend(x='bottomright',bty='n',legend=levels(types2_fact),pch=19,
  #        col=brewer.pal(length(unique(types2)),'Dark2'))
  
  plot(pca_btm$x[to_plot,3],pca_btm$x[to_plot,4],pch=19,cex=1.3,
       col=alpha((brewer.pal(length(unique(types2)),'Dark2')[c(1,2,3,6,4,5)])[as.numeric(types2_fact)],.5)[to_plot],
       xlab=paste('PCA3 prop dev =',100*round(pca_btm$sdev[3]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA4 prop dev =',100*round(pca_btm$sdev[4]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-16,12),ylim=c(-19,15))
  legend(x='bottomleft',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
  
  print(paste(round(range(pca_btm$x[,1]),1),round(range(pca_btm$x[,2]),1),round(range(pca_btm$x[,3]),1),
              round(range(pca_btm$x[,4]),1),sep='; ',collapse='{' ))
}

types_convert=list()
types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8_dpi',sep='')
types_convert[['SMI13']]=paste('Sassari_2013_dpi',sep='')
types_convert[['SLI13']]=paste('Teramo_2013_dpi',sep='')
types_convert[['SMI6']]=paste('Sassari_2006_dpi',sep='')
types_convert[['SLI6']]=paste('Teramo_2006_dpi',sep='')
types2=as.character(types_convert[types])
types2[dpis=='0']='uninfected'

pal=(brewer.pal(length(unique(types2)),'Dark2')[c(1,2,3,6,4,5)])
plot_ly(x=pca_btm$x[,1],y=pca_btm$x[,2],z=pca_btm$x[,3],
        color=factor(types2,levels=unique(unlist(types_convert))),colors=pal,alpha=0.85,
        text=~paste(colnames(mat2),sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))





###############
types2=paste('dpi',dpis)
types2[grepl('SMC|SLC',colnames(mat2))|types2=='dpi 0']='uninfected'

pal=brewer.pal(length(unique(types2)),'Dark2')
plot_ly(x=pca_btm$x[,1],y=pca_btm$x[,2],z=pca_btm$x[,3],
        color=factor(types2,levels=unique(c('uninfected',types2))),
        colors=pal,alpha=0.85,
        text=~paste(colnames(mat2),sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))


plot_ly(x=pca_btm$x[,4],y=pca_btm$x[,5],z=pca_btm$x[,6],
        color=factor(types2,levels=unique(c('uninfected',types2))),
        colors=pal,alpha=0.85,
        text=~paste(colnames(mat2),sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA6 (var prop=",round(100*PCA$sdev[6]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))
#####################
rfedat=t(mat3)


SMCdpi7=rep('no',nrow(rfedat))
#SMCdpi7[grepl('SMC',rownames(rfedat))&dpis=='7'|grepl('SMI3',rownames(rfedat))&dpis=='7']='yes'
SMCdpi7[grepl('SMI8',rownames(rfedat))&dpis=='1']='yes'



RF5=randomForest(y= as.factor(SMCdpi7[grepl('SMC',colnames(mat2))]),
                 x= rfedat[grepl('SMC',rownames(rfedat)),],importance=T,do.trace = F)


RF5=randomForest(y= as.factor(SMCdpi7),
                 x= rfedat,importance=T,do.trace = F)

dim(rfedat)

RF6=randomForest(y= as.factor(SMCdpi7),x= rfedat[,rownames(RF5$importance)[order(RF5$importance[,4],decreasing = T)[1:20]]],
                 importance=T,do.trace = F)
BTMs_to_plot=rownames(RF5$importance)[order(RF5$importance[,4],decreasing = T)[1:10]]
                         
###################
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
png(filename=paste('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn.8.3.21/SMI8_1.png',sep=''),width=1100,height=200*length(plots))

grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]]
             ,plots[[5]], plots[[6]],plots[[7]],plots[[8]],plots[[9]]
             ,plots[[10]]#,plots[[11]]#,plots[[12]]
             ,ncol=1)
dev.off()
print(length(BTMs_to_plot))
                      