dpi_levels=c('0','3')
gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
library(DESeq2)
types_headers=c('SMC')

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]

dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]
dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='_')))

dat_info=data.frame(dpi=as.factor(dpis),animal_name=as.factor(animal_names))

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = dat_info,
                              design= ~ animal_name+dpi)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="dpi_3_vs_0")
length(which(res$padj<0.05))
hist(res$log2FoldChange)
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="dpi_3_vs_0", type="apeglm")
hist(res$log2FoldChange)
length(which(res$padj<0.05))


##########################################


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
library(DESeq2)
types_headers=c('SLC')

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]

dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]
dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='_')))

dat_info=data.frame(dpi=as.factor(dpis),animal_name=as.factor(animal_names))

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = dat_info,
                              design= ~ animal_name+dpi)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="dpi_3_vs_0")
length(which(res$padj<0.05))
hist(res$log2FoldChange)
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="dpi_3_vs_0", type="apeglm")
hist(res$log2FoldChange)
length(which(res$padj<0.05))

#################################################################################################
#################################################################################################
dpi_levels=c('0','7')
gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
library(DESeq2)
types_headers=c('SMC')

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]

dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]
dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='_')))

dat_info=data.frame(dpi=as.factor(dpis),animal_name=as.factor(animal_names))

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = dat_info,
                              design= ~ animal_name+dpi)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="dpi_7_vs_0")
length(which(res$padj<0.05))
hist(res$log2FoldChange)
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="dpi_7_vs_0", type="apeglm")
hist(res$log2FoldChange)
length(which(res$padj<0.05))


##########################################


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
library(DESeq2)
types_headers=c('SLC')

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]

dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]
dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='_')))

dat_info=data.frame(dpi=as.factor(dpis),animal_name=as.factor(animal_names))

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = dat_info,
                              design= ~ animal_name+dpi)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="dpi_7_vs_0")
length(which(res$padj<0.05))
hist(res$log2FoldChange)
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="dpi_7_vs_0", type="apeglm")
hist(res$log2FoldChange)
length(which(res$padj<0.05))




