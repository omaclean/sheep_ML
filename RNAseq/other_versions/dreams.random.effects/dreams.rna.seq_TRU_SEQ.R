library('variancePartition')
library('BiocParallel')
library('edgeR')

counts=read.table('/home/oscar/Documents/sheep_megadata/tru_seq_29.6.21/All_Count.txt',header=T)

dpi=7

types_list=c("control","2013","2006","BTV8")
pairs=rbind(c('SLC','SLI6'),
            c('SMC','SMI8'))

for(dpi in c('1','3')){
  for(i in 1:nrow(pairs)){
    counts2=counts[,c(grep(paste('0','_dpi',sep=''),colnames(counts)),
                      grep(paste(dpi,'_dpi',sep=''),colnames(counts)))]
    rownames(counts2)=counts[,1]
    
    types_headers=pairs[i,]
    counts2=counts2[,grepl(types_headers[2],colnames(counts2))]
    
    counts2=counts2[rowSums(counts2>0)>7,]
    animal_names=as.character(sapply(colnames(counts2),function(x) 
      paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
    dpis=as.character(sapply(colnames(counts2),function(x) 
      strsplit(x,'\\.|_')[[1]][4]))
    inf_status=rep(0,ncol(counts2))
    
    inf_status[dpis>0&grepl(types_headers[2],colnames(counts2))]=1
    
    
    
    geneExpr = DGEList( counts2 )  
    
    geneExpr = calcNormFactors( geneExpr )
    
    meta_data=data.frame(inf_status,animal_names,dpis)
    rownames(meta_data)=colnames(counts2)
    form <- ~ inf_status + (1|animal_names) 
    
    meta_data
    
    # estimate weights using linear mixed model of dream
    vobjDream = voomWithDreamWeights( geneExpr, form, meta_data )
    fitmm = dream( vobjDream, form, meta_data )
    hits=topTable(fitmm,coef='inf_status',number=nrow(counts2))
    print(hits[1,])
    gene_names=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
    
    hits2=data.frame(human_gene_names=  gene_names$gene_names_new[match(rownames(hits),gene_names$ens_IDs)],
                     ensembl_ID=rownames(hits),FDR=hits$adj.P.Val ,logFC=hits$logFC)
    write.csv(hits,file=paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQ',
                              pairs[i,2],'dpi',dpi,'.csv',sep=''))
  }
}


for(dpi in c('1','3')){
  for(i in 1:nrow(pairs)){
    counts2=counts[,c(grep(paste('0','_dpi',sep=''),colnames(counts)),
                      grep(paste(dpi,'_dpi',sep=''),colnames(counts)))]
    rownames(counts2)=counts[,1]
    
    types_headers=pairs[i,]
    counts2=counts2[,grepl(types_headers[2],colnames(counts2))]
    
    counts2=counts2[rowSums(counts2>0)>7,]
    animal_names=as.character(sapply(colnames(counts2),function(x) 
      paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
    dpis=as.character(sapply(colnames(counts2),function(x) 
      strsplit(x,'\\.|_')[[1]][4]))
    inf_status=rep(0,ncol(counts2))
    
    inf_status[dpis>0&grepl(types_headers[2],colnames(counts2))]=1
    
    
    
    geneExpr = DGEList( counts2 )  
    
    geneExpr = calcNormFactors( geneExpr )
    
    meta_data=data.frame(inf_status,animal_names,dpis)
    rownames(meta_data)=colnames(counts2)
    form <- ~ inf_status #+ (1|animal_names) 
    
    d3=DGEList(counts=counts2,group= inf_status)
    design <- model.matrix(~inf_status)
    d3=calcNormFactors(d3)
    d3 <- estimateDisp(d3,design)
    d3 <- glmQLFit(d3,design)
    de2=glmQLFTest(d3,contrast=c(0,1))
   hits= topTags(de2, n=nrow(counts))
    
    # estimate weights using linear mixed model of dream
    # vobjDream = voomWithDreamWeights( geneExpr, form, meta_data )
    # # fitmm = dream( vobjDream, form, meta_data )
    # 
    #  hits=topTable(fitmm,coef='inf_status',number=nrow(counts2))
   # print(hits[1:100,])
    print(c(dpi,pairs[i,],length(which(hits$table$FDR<0.05))))
    gene_names=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
    
    hits2=data.frame(human_gene_names=  gene_names$gene_names_new[match(rownames(hits$table),gene_names$ens_IDs)],
                     ensembl_ID=rownames(hits$table),FDR=hits$table$FDR ,logFC=hits$table$logFC)
    write.csv(hits,file=paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQ_norandom',
                              pairs[i,2],'dpi',dpi,'.csv',sep=''))
  }
}





for(dpi in c('1','3')){
  for(i in 1:nrow(pairs)){
    counts2=counts[,c(grep(paste('0','_dpi',sep=''),colnames(counts)),
                      grep(paste(dpi,'_dpi',sep=''),colnames(counts)))]
    rownames(counts2)=counts[,1]
    
    types_headers=pairs[i,]
    counts2=counts2[,grepl(types_headers[2],colnames(counts2))]
    
    counts2=counts2[rowSums(counts2>0)>7,]
    animal_names=as.character(sapply(colnames(counts2),function(x) 
      paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
    dpis=as.character(sapply(colnames(counts2),function(x) 
      strsplit(x,'\\.|_')[[1]][4]))
    inf_status=rep(0,ncol(counts2))
    
    inf_status[dpis>0&grepl(types_headers[2],colnames(counts2))]=1
    
    
    
    geneExpr = DGEList( counts2 )  
    
    geneExpr = calcNormFactors( geneExpr )
    
    meta_data=data.frame(inf_status,animal_names,dpis)
    rownames(meta_data)=colnames(counts2)
    form <- ~ inf_status #+ (1|animal_names) 
    
    d3=DGEList(counts=counts2,group= inf_status)
    design <- model.matrix(~inf_status+animal_names)
    d3=calcNormFactors(d3)
    d3 <- estimateDisp(d3,design)
    d3 <- glmQLFit(d3,design)
    de2=glmQLFTest(d3,contrast=c(0,1,rep(0,length(unique(animal_names))-1)))
    hits= topTags(de2, n=nrow(counts))
    
    # estimate weights using linear mixed model of dream
    # vobjDream = voomWithDreamWeights( geneExpr, form, meta_data )
    # # fitmm = dream( vobjDream, form, meta_data )
    # 
    #  hits=topTable(fitmm,coef='inf_status',number=nrow(counts2))
    print(c(dpi,pairs[i,],length(which(hits$table$FDR<0.05))))
    gene_names=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
    
    hits2=data.frame(human_gene_names=  gene_names$gene_names_new[match(rownames(hits$table),gene_names$ens_IDs)],
                     ensembl_ID=rownames(hits$table),FDR=hits$table$FDR ,logFC=hits$table$logFC)
    write.csv(hits,file=paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQ_norandom_fixed_indiv',
                              pairs[i,2],'dpi',dpi,'.csv',sep=''))
  }
}



isexpr = rowSums(counts2>0) >= 7


geneExpr = calcNormFactors( geneExpr )


form <- ~ Disease + (1|Individual) 

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, metadata )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata )

topTable(fitmm,coef='Disease1',number=nrow(geneExpr))
