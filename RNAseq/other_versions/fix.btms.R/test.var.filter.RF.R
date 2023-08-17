library(psych);library(ggplot2);library(caret);library(randomForest)

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
mat2_copy=mat2


mat2_cor=cor(t(mat2))
cors=rep(0,nrow(mat2_cor))

for(i in 1:(ncol(mat2_cor)-1)){
  count=0
  for(j in (i+1):ncol(mat2_cor)){
    if(mat2_cor[i,j]>0.9){
      cors[i]=cors[i]+1
      rownames(mat2)[j]=paste('mean(',rownames(mat2)[i],rownames(mat2)[j],')')
      mat2[j,]=sapply(1:ncol(mat2),function(x) mean(mat2[j,x],mat2[i,x]))
    }
  }
}
table(cors)
mat2=mat2[cors==0,]

mat3=mat2


####
# convert colnames from counts into 'types' &extract DPI
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
types_convert=list()
types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8_dpi'
types_convert[['SMI13']]=types_convert[['SLI13']]='2013_dpi'
types_convert[['SMI6']]=types_convert[['SLI6']]='2006_dpi'
types2=as.character(types_convert[types])
dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
types3=paste(types2,'_dpi',dpis)


############################################################
################################################################
################################################################################################
dpi_plot='7'
conds=c('SMI8','SLI6','SMI6','SLI13','SMI13')

to_plot=(dpis==dpi_plot)|(types2=='uninfected'|dpis=='0')
animal_states=rep(1,length(which(to_plot)))
rfedat=t(mat3[,to_plot])
count=0
for(cond in conds){
  count=count+1
  animal_states[grepl(cond,rownames(rfedat))&dpis[to_plot]==dpi_plot]=count+1 
}


BTV8DPI1=factor(c('uninf',conds)[animal_states],levels=c('uninf',conds))






RF=randomForest(y= BTV8DPI1,x= rfedat,
                importance=T,do.trace = F)

RF2=randomForest(y= BTV8DPI1,x= rfedat[,rownames(RF$importance)[order(RF$importance[,4],decreasing = T)[1:20]]],
                 importance=T,do.trace = F)


##############################
library(psych);library(ggplot2);library(caret);library(randomForest)

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

for(i in 1:nrow(counts)){
  counts[i,2:ncol(counts)]=counts[i,2:ncol(counts)]/sum(counts[i,2:ncol(counts)])
}

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
mat2_copy=mat2

mat2_cor=cor(t(mat2))
cors=rep(0,nrow(mat2_cor))

for(i in 1:(ncol(mat2_cor)-1)){
  count=0
  for(j in (i+1):ncol(mat2_cor)){
    if(mat2_cor[i,j]>0.9){
      cors[i]=cors[i]+1
      rownames(mat2)[j]=paste('mean(',rownames(mat2)[i],rownames(mat2)[j],')')
      mat2[j,]=sapply(1:ncol(mat2),function(x) mean(mat2[j,x],mat2[i,x]))
    }
  }
}
table(cors)
mat2=mat2[cors==0,]

mat3=mat2


####
# convert colnames from counts into 'types' &extract DPI


to_plot=(dpis==dpi_plot)|(types2=='uninfected'|dpis=='0')
BTV8DPI1=factor(c('uninf',conds)[animal_states],levels=c('uninf',conds))
################################################################################################
dpi_plot='7'
conds=c('SMI8','SLI6','SMI6','SLI13','SMI13')

to_plot=(dpis==dpi_plot)|(types2=='uninfected'|dpis=='0')

animal_states=rep(1,length(which(to_plot)))
rfedat=t(mat3[,to_plot])
count=0
for(cond in conds){
  count=count+1
  animal_states[grepl(cond,rownames(rfedat))&dpis[to_plot]==dpi_plot]=count+1 
}
BTV8DPI1=factor(c('uninf',conds)[animal_states],levels=c('uninf',conds))






RF3=randomForest(y= BTV8DPI1,x= rfedat,
                 importance=T,do.trace = F)

RF4=randomForest(y= BTV8DPI1,x= rfedat[,rownames(RF3$importance)[order(RF3$importance[,4],decreasing = T)[1:20]]],
                 importance=T,do.trace = F)


##############################
infected=!(types2=='uninfected'|dpis=='0')[to_plot]
excess_var=rep(0,ncol(rfedat))
for(i in 1:ncol(rfedat)){
  if(var(rfedat[!infected,i])>var(rfedat[,i])){
    excess_var[i]=1
  }
}
sum(excess_var)/length(excess_var)
rfedat=rfedat[,!excess_var]
##


RF5=randomForest(y= as.factor(BTV8DPI1),x= rfedat,
                 importance=T,do.trace = F)

RF6=randomForest(y= BTV8DPI1,x= rfedat[,rownames(RF5$importance)[order(RF5$importance[,4],decreasing = T)[1:20]]],
                 importance=T,do.trace = F)


#RF$confusion
RF2$confusion
#RF3$confusion
RF4$confusion
#RF5$confusion
RF6$confusion
