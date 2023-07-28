
test=read.delim('/home/oscar/Downloads/BTM_sheep.txt',sep  ='\t',fill = TRUE)
write.csv(test[,2:ncol(test)],'~/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs.fromArtur.csv',
               row.names=F)

test$ENS_ID=gsub(' ','',test$ENS_ID)
test$ENS_ID=gsub('\n','',test$ENS_ID)


grep(' ',test$ENS_ID,value=T)[1]

hist(as.numeric(sapply(test$ENS_ID,function(x) length(unique(strsplit(x,',')[[1]]))))/test$NR_F.OUND,breaks=280)
     

length(which(as.numeric(sapply(test$ENS_ID,function(x) length(unique(strsplit(x,',')[[1]]))))/test$NR_F.OUND==1))
dim(test)

as.numeric(sapply(test$NR_F.OUND,function(x) length(strsplit(x,',')[[1]])))


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

dpi=7
counts=counts[,c(1,grep(paste(dpi,'_dpi',sep=''),colnames(counts)))]

types_list=c("control","2013","2006","BTV8")
types_headers=c('SLC|SMC','SLI13|SMI13','SLI6|SMI6','SMI8')
#types_array=
animal_names=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) 
    paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))
  
inf_status=factor(inf_status,levels=types_list)
counts2=counts[,2:ncol(counts)]
counts2=counts2[rowSums(counts2>0)>7,]
d3=DGEList(counts=counts2,group= inf_status)
design <- model.matrix(~inf_status)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)

library(edgeR)
for(i in 1:(length(types_list)-1)){
  for(j in (i+1):length(types_list)){
    
 
    de2 = exactTest(d3, pair = c(types_list[i],types_list[j]))
    test= topTags(de2, n=nrow(counts))
    print(c(length(which(test$table$FDR<0.05)),types_list[i],types_list[j]))
  }
}


#################################################################################
################

dpi=7
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

counts=counts[,c(1,grep(paste(dpi,'_dpi',sep=''),colnames(counts)))]

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]
###########
types_list=c("control","2013","2006","BTV8")
types_headers=c('SLC|SMC','SLI13|SMI13','SLI6|SMI6','SMI8')
#types_array=
animal_names=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) 
  types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))

inf_status=factor(inf_status,levels=types_list)
locations=c("short","long")



location=factor(locations[sapply(colnames(counts2),function(x) 
  grep(substr(strsplit(x,'\\.')[[1]][1],1,2) ,c('SM','SL') ))],levels=locations)

d3=DGEList(counts=counts2,group= inf_status)
design <- model.matrix(~inf_status+location)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))
for(i in 1:(length(types_list)-1)){
  for(j in (i+1):length(types_list)){
   contrast=as.numeric(types_list%in%types_list[j])
   if(types_list[i]!='control'){
     contrast[which(types_list%in%types_list[i])]=-1
   }

    de2=glmQLFTest(d3,contrast=c(contrast,0))
    test= topTags(de2, n=nrow(counts))
    print(c(length(which(test$table$FDR<0.05)),types_list[i],types_list[j]))
  }
}
dim(counts2)


#### no location


d3=DGEList(counts=counts2,group= inf_status)
design <- model.matrix(~inf_status)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))
for(i in 1:(length(types_list)-1)){
  for(j in (i+1):length(types_list)){
    contrast=as.numeric(types_list%in%types_list[j])
    if(types_list[i]!='control'){
      contrast[which(types_list%in%types_list[i])]=-1
    }
    
    de2=glmQLFTest(d3,contrast=c(contrast))
    test= topTags(de2, n=nrow(counts))
    print(c(length(which(test$table$FDR<0.05)),types_list[i],types_list[j]))
  }
}
dim(counts2)


#################################################################################
################

dpi=0
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

counts=counts[,c(1,grep(paste(dpi,'_dpi|SMC|SLC',sep=''),colnames(counts)))]
types_list=c("short","long")
types_headers=c('SL','SM')
#####################################################
counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]
#types_array=
animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts2),function(x) types_list[grep(substr(strsplit(x,'\\.')[[1]][1],1,2),types_headers)]))

inf_status=factor(inf_status,levels=types_list)


d3=DGEList(counts=counts2,group= inf_status)
design <- model.matrix(~inf_status)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)

library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))
for(i in 1:(length(types_list)-1)){
  for(j in (i+1):length(types_list)){
    de2 = exactTest(d3, pair = c(types_list[i],types_list[j]))
    test= topTags(de2, n=nrow(counts))
    print(c(length(which(test$table$FDR<0.05)),types_list[i],types_list[j]))
  }
}











###############

dpi=0
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

counts=counts[,c(1,grep(paste(dpi,'_dpi|0_dpi',sep=''),colnames(counts)))]

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]
###########
types_list=c("control","2013","2006","BTV8")
types_headers=c('SLC|SMC','SLI13|SMI13','SLI6|SMI6','SMI8')
#types_array=
animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) 
  types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))
inf_status[grep('0_dpi',)]

inf_status=factor(inf_status,levels=types_list)
locations=c("short","long")



location=factor(locations[sapply(colnames(counts2),function(x) 
  grep(substr(strsplit(x,'\\.')[[1]][1],1,2) ,c('SM','SL') ))],levels=locations)

d3=DGEList(counts=counts2,group= inf_status)
design <- model.matrix(~location)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))
for(i in 1:(length(types_list)-1)){
  for(j in (i+1):length(types_list)){
    contrast=as.numeric(types_list%in%types_list[j])
    if(types_list[i]!='control'){
      contrast[which(types_list%in%types_list[i])]=-1
    }
    
    de2=glmQLFTest(d3,contrast=c(0,1))
    test= topTags(de2, n=nrow(counts))
    print(c(length(which(test$table$FDR<0.05)),types_list[i],types_list[j]))
  }
}
dim(counts2)


