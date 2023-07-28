setwd("/home/oscar/Documents/sheep_megadata/11.11.20/")
library(ggplot2)
library(plotly)
library(gridExtra)
library(RColorBrewer)
library(Amelia)
library(reshape2)

data=read.csv('/home/oscar/Documents/sheep_megadata/17.2.21/combined/combined_sheet.csv')

to_remove=c("tongue_pos_cells" ,"lung_pos_cells","infected_skin_superficial_dermis_pos_cells",
            "infected_skin_deep_dermis_pos_cells",            "non.infected_skin_superficial_dermis_pos_cells",
            "spleen_pos_cells",            "lymph_nodes_cortex_pos_cells" ,
            "rectal..temperature",            "number.of.follicles.mean.per.prescapular.ln",
            "percentpositivenucleiFOXP3cortexnofollicles",            "percentpositivenucleiFOXP3wholefollicle",
            "percentpositivenucleiFOXP3wholelymphnode",            "percentpositivenucleiFOXP3medulla",
            "percentpositivenucleiFOXP3folliclecentres")
data=data[,!colnames(data)%in%to_remove]
data=data[data$dpi=='7',]
length(which(is.na(data[,3:ncol(data)])))/(ncol((data[,3:ncol(data)]))*nrow((data[,3:ncol(data)])))

cols=c("#20B2AA",brewer.pal(6,"Reds")[1+c(1,3,5,4,2)])


pca_dat=data[,3:ncol(data)]
length(which(is.na(data[,3:ncol(data)])))/(ncol((data[,3:ncol(data)]))*nrow((data[,3:ncol(data)])))
for(i in 1:ncol(pca_dat)){
  nas=is.na(pca_dat[,i])
  nas[pca_dat[,i]==""]=T
  if(any(nas)){
    pca_dat[nas,i]=mean(as.numeric((pca_dat[,i])[!nas]))
  }
  pca_dat[,i]=as.numeric(pca_dat[,i])
}
rownames(pca_dat)=paste(data$ID,data$dpi,sep='-')


################################BTMS
  if(!exists('BTM_check')){
  library(psych);library(ggplot2);library(caret);library(randomForest)
    dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
    
  for(animal in unique(grep('SLC',dat$ID,value=T))){
    for(j in 3:ncol(dat)){
      if(is.na(dat[dat$ID==animal&dat$dpi=='7',j])&!is.na(dat[dat$ID==animal&dat$dpi=='21',j])){
        dat[dat$ID==animal&dat$dpi=='7',j]=dat[dat$ID==animal&dat$dpi=='21',j]
      }
    }
  }
  
  counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
  counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
  #redundant_BTMS=c("Rho.GTPase.cycle..M4.14.","mitotic.cell.division..M6.",
  #                 "enriched.in.NK.cells..KIR.cluster...M61.1.","enriched.in.B.cells..VI...M69.")
  
  dat=dat[,!colnames(dat)%in%redundant_BTMS]
  
  
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
  
  dim(mat2[is.na(mat2[,1]),])
  }
BTM_check=T
###########################
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
par(mfrow=c(3,2),mar=c(4,4,1,1))
mat3=mat2[rowSums(mat2>0)>1,]

#####################################

pca_btm=prcomp(t(mat3),retx=T,scale.=T,center=T)
pca_met=prcomp((pca_dat),retx=T,scale.=T,center=T)
#pca_met$x=pca_met$x[grepl('7',rownames(pca_met$x)),]
pca_met$x=pca_met$x[data$dpi=='7',]
pca_met$x=pca_met$x[c(grep('SLC',rownames(pca_met$x)),grep('SMC',rownames(pca_met$x)),
                      grep('SMI13',rownames(pca_met$x)),grep('SLI13',rownames(pca_met$x)),
                      grep('SMI6',rownames(pca_met$x)),grep('SLI6',rownames(pca_met$x)),
                      grep('SMI8',rownames(pca_met$x))),]
dim(pca_btm$x)


rownames(pca_met$x)#=paste(rownames(pca_met$x),data$dpi,sep='-')
#rownames(pca_met$x)[grepl('SLC',rownames(pca_met$x))]=gsub('-21','-7',rownames(pca_met$x)[grepl('SLC',rownames(pca_met$x))])
rownames(pca_btm$x)=gsub('\\.','-',sapply(rownames(pca_btm$x),
                                        function(x) paste(strsplit(x,'_')[[1]][c(1,3)],
                                                          collapse='-')))
dim(pca_btm$x[rownames(pca_btm$x)%in%rownames(pca_met$x),])

pca_btm2=pca_btm$x[rownames(pca_btm$x)%in%rownames(pca_met$x),]
pca_met2=pca_met$x[rownames(pca_met$x)%in%rownames(pca_btm2),]
pca_btm2=pca_btm2[match(rownames(pca_met2),rownames(pca_btm2)),]

col_levs=unique(sapply(rownames(pca_met2),function(x) strsplit(x,'-')[[1]][1]))
cols_i=rep(0,nrow(pca_met2))
for(i in 1:nrow(pca_met2)){
  cols_i[i]=which(col_levs== strsplit(rownames(pca_met2)[i],'-')[[1]][1])
}
cols=brewer.pal(length(col_levs),'Dark2')[c(1,5,6,7,3,4,2)]

library(scales)
par(mfrow=c(1,2),mar=c(5,4,2,1))


plot(pca_btm2[,1],pca_met2[,1],pch=19,col=scales::alpha(cols[cols_i],.7),
     ylab='metabolic PCA axis 1',xlab='RNA seq BTM axis 1')
#legend(x='topright',bty='n',legend=col_levs,pch=19,col=cols)
plot(pca_btm2[,3],pca_met2[,1],pch=19,col=scales::alpha(cols[cols_i],.7),
     ylab='metabolic PCA axis 1',xlab='RNA seq BTM axis 3')
legend(x='bottomleft',bty='n',legend=col_levs,pch=19,col=cols)

par(mfrow=c(2,2),mar=c(5,4,4,1))

plot(pca_btm2[,1],pca_met2[,2],pch=19,col=scales::alpha(cols[cols_i],.8))
#legend(x='topright',bty='n',legend=col_levs,pch=19,col=cols)
plot(pca_btm2[,3],pca_met2[,2],pch=19,col=scales::alpha(cols[cols_i],.8))
legend(x='bottomleft',bty='n',legend=col_levs,pch=19,col=cols)

plot(pca_btm2[,2],pca_met2[,2],pch=19,col=scales::alpha(cols[cols_i],.8))
#legend(x='topright',bty='n',legend=col_levs,pch=19,col=cols)
plot(pca_btm2[,4],pca_met2[,2],pch=19,col=scales::alpha(cols[cols_i],.8))
legend(x='bottomleft',bty='n',legend=col_levs,pch=19,col=cols)

par(mfrow=c(1,2))

plot(pca_met2[grepl('SL',rownames(pca_met$x)),1],
     pca_dat$Proteine.Totali[grepl('SL',rownames(pca_met$x))],
     col=scales::alpha(cols[cols_i],.8)[grepl('SL',rownames(pca_met$x))],pch=19)
plot(pca_met2[grepl('SM',rownames(pca_met$x)),1],
     pca_dat$Proteine.Totali[grepl('SM',rownames(pca_met$x))],
     col=scales::alpha(cols[cols_i],.8)[grepl('SM',rownames(pca_met$x))],pch=19)
legend(x='topright',bty='n',legend=col_levs,pch=19,col=cols)

for(i in 1:10){
  for(j in 1:10){
    if(abs(cor(pca_met2[,i],pca_btm2[,j]))>.5){
      print(c(i,j,cor(pca_met2[,i],pca_btm2[,j])))
    }
  }
}
cor(pca_btm2[,1],pca_btm2[,3])
plot(pca_btm2[,1],pca_btm2[,3])
