
####four examples 

library(edgeR)


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
###########
types_list=c("2006")
types_headers=c('SLI6')


counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]


animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts2),function(x) 
  types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))
inf_status[grep('0_dpi',colnames(counts2))]="control"

dpi_levels=c('0','1','3','7')
dpis=factor(sapply(colnames(counts2),function(x) tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1)),levels=dpi_levels)

inf_status=factor(inf_status,levels=types_list)
locations=c("short","long")

location=factor(locations[sapply(colnames(counts2),function(x) 
  grep(substr(strsplit(x,'\\.')[[1]][1],1,2) ,c('SM','SL') ))],levels=locations)

d3=DGEList(counts=counts2,group= inf_status)
animal_names2=as.factor(animal_names)
################################################
design <- model.matrix(~dpis)
##################
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results=matrix(ncol=3,nrow=choose(length(unique(dpis)),2)+2)
results[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )

cont=rep(0,length(design[1,]))
cont[length(cont)]=1
de2=glmQLFTest(d3,contrast=cont)
test= topTags(de2, n=nrow(counts))

results[2,]=c(length(which(test$table$FDR<0.05)),'baseline',tail(colnames(design),1))
count=count_init=3
genes_hits=list()
for(i in 1:(length(unique(dpis))-1)){
  for(j in (i+1):length(unique(dpis))){
    contrast=as.numeric(dpi_levels%in%dpis[j])
    if(colnames(design)[i]!='(Intercept)'){
      contrast[i]=-1
    }
    print(c(colnames(design)[i],contrast,rep(0,length(design[1,])-length(contrast))))
    de2=glmQLFTest(d3,contrast=c(contrast,rep(0,length(design[1,])-length(contrast))))
    test= topTags(de2, n=nrow(counts))
    genes_hits[[count-count_init+1]]=rownames(test$table)[test$table$FDR<0.05]
    results[count,]=c(length(which(test$table$FDR<0.05)),colnames(design)[i],colnames(design)[j])
    count=count+1
  }
}
print(results)



###############################################

d3=DGEList(counts=counts2)
animal_names2=as.factor(animal_names)
################################################
design <- model.matrix(~dpis+animal_names2)
##################
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results2=matrix(ncol=3,nrow=choose(length(unique(dpis)),2)+2)
results2[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )

cont=rep(0,length(design[1,]))
cont[length(cont)]=1
de2=glmQLFTest(d3,contrast=cont)
test= topTags(de2, n=nrow(counts))

results2[2,]=c(length(which(test$table$FDR<0.05)),'baseline',tail(colnames(design),1))
count=count_init=3
genes_hits_indiv=list()
for(i in 1:(length(unique(dpis))-1)){
  for(j in (i+1):length(unique(dpis))){
    contrast=rep(0,length(design[1,]))
      contrast[which(dpi_levels%in%dpis[j])]=1
    if(colnames(design)[i]!='(Intercept)'){
      contrast[i]=-1
    }
    print(c(colnames(design)[i],contrast,rep(0,length(design[1,])-length(contrast))))
    de2=glmQLFTest(d3,contrast=c(contrast,rep(0,length(design[1,])-length(contrast))))
    test= topTags(de2, n=nrow(counts))
    genes_hits_indiv[[count-count_init+1]]=rownames(test$table)[test$table$FDR<0.05]
    results2[count,]=c(length(which(test$table$FDR<0.05)),colnames(design)[i],colnames(design)[j])
    count=count+1
  }
}
print(results)
print(results2)




###################################################################################################
library(edgeR)


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
###########
types_list=c("2006")
types_headers=c('SMI6')


counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]


animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts2),function(x) 
  types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))
inf_status[grep('0_dpi',colnames(counts2))]="control"

dpi_levels=c('0','1','3','7')
dpis=factor(sapply(colnames(counts2),function(x) tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1)),levels=dpi_levels)

inf_status=factor(inf_status,levels=types_list)
locations=c("short","long")

location=factor(locations[sapply(colnames(counts2),function(x) 
  grep(substr(strsplit(x,'\\.')[[1]][1],1,2) ,c('SM','SL') ))],levels=locations)

d3=DGEList(counts=counts2,group= inf_status)
animal_names2=as.factor(animal_names)
################################################
design <- model.matrix(~dpis)
##################
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results_1=matrix(ncol=3,nrow=choose(length(unique(dpis)),2)+2)
results_1[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )

cont=rep(0,length(design[1,]))
cont[length(cont)]=1
de2=glmQLFTest(d3,contrast=cont)
test= topTags(de2, n=nrow(counts))

results_1[2,]=c(length(which(test$table$FDR<0.05)),'baseline',tail(colnames(design),1))
count=count_init=3
genes_hits2=list()
for(i in 1:(length(unique(dpis))-1)){
  for(j in (i+1):length(unique(dpis))){
    contrast=as.numeric(dpi_levels%in%dpis[j])
    if(colnames(design)[i]!='(Intercept)'){
      contrast[i]=-1
    }
    print(c(colnames(design)[i],contrast,rep(0,length(design[1,])-length(contrast))))
    de2=glmQLFTest(d3,contrast=c(contrast,rep(0,length(design[1,])-length(contrast))))
    test= topTags(de2, n=nrow(counts))
    genes_hits2[[count-count_init+1]]=rownames(test$table)[test$table$FDR<0.05]
    results_1[count,]=c(length(which(test$table$FDR<0.05)),colnames(design)[i],colnames(design)[j])
    count=count+1
  }
}
print(results_1)



###############################################

d3=DGEList(counts=counts2)
animal_names2=as.factor(animal_names)
################################################
design <- model.matrix(~dpis+animal_names2)
##################
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results_12=matrix(ncol=3,nrow=choose(length(unique(dpis)),2)+2)
results_12[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )

cont=rep(0,length(design[1,]))
cont[length(cont)]=1
de2=glmQLFTest(d3,contrast=cont)
test= topTags(de2, n=nrow(counts))

results_12[2,]=c(length(which(test$table$FDR<0.05)),'baseline',tail(colnames(design),1))
count=count_init=3
genes_hits2_indiv=list()
for(i in 1:(length(unique(dpis))-1)){
  for(j in (i+1):length(unique(dpis))){
    contrast=rep(0,length(design[1,]))
    contrast[which(dpi_levels%in%dpis[j])]=1
    if(colnames(design)[i]!='(Intercept)'){
      contrast[i]=-1
    }
    print(c(colnames(design)[i],contrast,rep(0,length(design[1,])-length(contrast))))
    de2=glmQLFTest(d3,contrast=c(contrast,rep(0,length(design[1,])-length(contrast))))
    test= topTags(de2, n=nrow(counts))
    genes_hits2_indiv[[count-count_init+1]]=rownames(test$table)[test$table$FDR<0.05]
    results_12[count,]=c(length(which(test$table$FDR<0.05)),colnames(design)[i],colnames(design)[j])
    count=count+1
  }
}
print(results_1)
print(results_12)


library(VennDiagram)

for(i in 3:nrow(results)){
  comp=paste(results[i,2:3],collapse='_v_')
  if(comp==paste(results2[i,2:3],collapse='_v_')&
     comp==paste(results_12[i,2:3],collapse='_v_')&
     comp==paste(results_1[i,2:3],collapse='_v_')){
    
    venns=list(genes_hits[[i-2]],genes_hits_indiv[[i-2]],
               genes_hits2[[i-2]],genes_hits2_indiv[[i-2]]
    )
    venn.diagram(x=venns,
                 category.names=paste(c('Teramo simple\n','Teramo individual\n',
                                           'Sassari simple\n','Sassari individual\n'),rep(' n=',4),
                                        unlist(lapply(venns,function(x) length(x))),sep=''),filename = 
                   paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/',comp,'.png'),
                 )
    
  }
}

