
library(edgeR)
## two comparisons 1) control animal A vs B. 2) control animal dpi 0 vs 7

#1)
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
###########
types_list=c("control")
types_headers=c('SLC')


counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]


animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))

counts2=counts2[,animal_names%in%unique(animal_names)[1:2]]
animal_names=as.factor(as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-'))))


d3=DGEList(counts=counts2,group= animal_names)

################################################
design <- model.matrix(~animal_names)
##################
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results=matrix(ncol=3,nrow=2)
results[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )



genes_hits=list()
contrast=rep(0,ncol(design))
contrast[2]=1
de2=glmQLFTest(d3,contrast=contrast)
test= topTags(de2, n=nrow(counts2))
genes_hits[[1]]=rownames(test$table)[test$table$FDR<0.05]
results[2,]=c(length(which(test$table$FDR<0.05)),colnames(design)[1],colnames(design)[2])


print(results)


##########################################################################################




#2)
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
###########
types_list=c("control")
types_headers=c('SLC')
genes_hits=list()
counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]

dpi_levels=c('0'#,'1','3'
             ,'3')
dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))


counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]

dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

dpis=factor(dpis,levels=dpi_levels)
animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))

animal_names2=as.factor(animal_names)

################################################
design <- model.matrix(~dpis+animal_names2)
##################
d3=DGEList(counts=counts2)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=paste(animal_names,'dpi',as.character(dpis))
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results=matrix(ncol=3,nrow=2)
results[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )


contrast=rep(0,ncol(design))
contrast[2]=1
de2=glmQLFTest(d3,contrast=contrast)
test= topTags(de2, n=nrow(counts2))
FDR_threshold=0.05
genes_hits[[4]]=rownames(test$table)
genes_hits[[5]]=test$table$FDR
genes_hits[[6]]=test$table$logFC
results[2,]=c(length(which(test$table$FDR<FDR_threshold)),colnames(design)[1],colnames(design)[2])


print(results)
print(0.01*nrow(counts2))



plot(genes_hits[[2]],genes_hits[[5]][match(genes_hits[[1]],genes_hits[[4]])])
dim(counts2)
plot(genes_hits[[3]],genes_hits[[6]][match(genes_hits[[1]],genes_hits[[4]])])

head(cbind(genes_hits[[1]],genes_hits[[4]][match(genes_hits[[1]],genes_hits[[4]])]))






