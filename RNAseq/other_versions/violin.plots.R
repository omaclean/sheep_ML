library(edgeR);library(scales)
genes_hits_logFC=genes_hits_genes=genes_hits_FDR=list()
dpi_levels=c('1','3','7')
cond_levels=c('SMC','SLC','SMI6','SLI6','SMI13','SLI13','SMI8')
location1='Sassari'
types_list1=c("control")
types_headers1=c('SMC')

location2='Teramo'
types_list2=c("control")
types_headers2=c('SLC')

for(dpi in dpi_levels){
  ###########################
  genes_hits_genes[[dpi]]=genes_hits_FDR[[dpi]]=genes_hits_logFC[[dpi]]=list()
  dpi_levels=c('0',dpi)
  for(condition in cond_levels){
    
  
  
  gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)
  
  counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
  rownames(counts)=counts$Geneid
  ###########
  types_list=condition
  types_headers=condition
  
  counts2=counts[,2:ncol(counts)]
  counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
  counts2=counts2[rowSums(counts2>0)>7,]
  colnames(counts2)
  
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
  
  
  contrast=rep(0,ncol(design))
  contrast[2]=1
  de2=glmQLFTest(d3,contrast=contrast)
  test= topTags(de2, n=nrow(counts2))
  FDR_threshold=0.05
  genes_hits_logFC[[dpi]][[condition]]=test$table$logFC[test$table$FDR<0.05]
  genes_hits_FDR[[dpi]][[condition]]=test$table$FDR[test$table$FDR<0.05]
  genes_hits_genes[[dpi]][[condition]]=rownames(test$table)[test$table$FDR<0.05]
  }
}
dpi_levels=c('1','3','7')
library(ggplot2);library(gridExtra)
for(dpi in c('1','3','7')){
  data=data.frame(logFC=NULL,gene_name=NULL,condition=NULL)
  for(cond in cond_levels){
    data=rbind(data,cbind(genes_hits_logFC[[dpi]][[cond]],genes_hits_genes[[dpi]][[cond]],
                          rep(paste(cond,'\nn=',length(genes_hits_logFC[[dpi]][[cond]])),
                              length(genes_hits_logFC[[dpi]][[cond]]))))
  }
  data=rbind(data,cbind(rnorm(1650),rep('crap\nn=1650',1650),rep('crap\nn=1650',1650)))
  dim(data)
  colnames(data)=c('logFC','gene_name','condition')
  data$logFC=as.numeric(as.character(data$logFC))
  str(data)
  if(dpi==dpi_levels[1]){
    g1=ggplot(data,aes(y=logFC,x=condition))+geom_violin(scale='count')+
      theme_bw()+ggtitle('dpi1 vs dpi0')+ylim(-5,10)
  }
  if(dpi==dpi_levels[2]){
    g2=ggplot(data,aes(y=logFC,x=condition))+geom_violin(scale='count')+
      theme_bw()+ggtitle('dpi3 vs dpi0')+ylim(-5,10)
  }
  if(dpi==dpi_levels[3]){
    g3=ggplot(data,aes(y=logFC,x=condition))+geom_violin(scale='count')+
      theme_bw()+ggtitle('dpi7 vs dpi0')+ylim(-5,10)
  }
}
grid.arrange(g1,g2,g3)
dim(data)
