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




mat2_cor=cor(t(mat2))
hit_len=length(mat2_cor[mat2_cor>0.9& mat2_cor<1])

#################
while(hit_len>0){
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
  mat2_save=mat2
  ########################
  
  mat2_cor=cor((mat2))
  pairs=matrix(ncol=2,nrow=0)
  for(i in 1:(ncol(mat2_cor)-1)){
    count=0
    for(j in (i+1):ncol(mat2_cor)){
      if(mat2_cor[i,j]>0.9){
        pairs=rbind(pairs,c(i,j))
      }
    }
  }
  
  hit_len=length(mat2_cor[mat2_cor>0.9& mat2_cor<1])
}

write.csv(mat2,'/home/oscar/Documents/sheep_megadata/13.4.21/dat_plots_filt/BTMs.final.csv')
#########################################
