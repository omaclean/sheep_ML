library(edgeR);library(scales)

ISGs=read.csv('/home/oscar/Documents/Alex_RNA_seq/Ov_ISGlist.csv')
orths=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')

genes_hits=list()
dpi_levels=c('0',NA)
location1='Sassari'
location2='Teramo'

types_loop1=c('C','I13','I6','I8')
types_loop2=c('C','I13','I6','I13')
types_list1_loop=c("control",'2013','2006','BTV8')
types_list2_loop=c("control",'2013','2006','2013')
for(lev2 in c('1','3',7)){
  dpi_levels[2]=lev2  
  for(i in 1:length(types_loop1)){
    
    types_headers1=paste('SM',types_loop1[i],sep='')
    types_headers2=paste('SL',types_loop2[i],sep='')
    types_list1=types_list1_loop[i]
    types_list2=types_list2_loop[i]
    
    ###########################
    
    gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)
    
    counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
    rownames(counts)=counts$Geneid
    ###########
    types_list=types_list1
    types_headers=types_headers1
    
    counts2=counts[,2:ncol(counts)]
    counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
    counts2=counts2[rowSums(counts2>0)>7,]

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
    
    
    genes_hits[[1]]=rownames(test$table)
    genes_hits[[2]]=test$table$FDR
    genes_hits[[3]]=test$table$logFC
    
    #2)
    ###########
    types_list=types_list2
    
    types_headers=types_headers2
    
    counts2=counts[,2:ncol(counts)]
    counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
    counts2=counts2[rowSums(counts2>0)>7,]
    
    
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
    
    
    genes_hits[[4]]=rownames(test$table)
    genes_hits[[5]]=test$table$FDR
    genes_hits[[6]]=test$table$logFC
    
    
  
    names_plot=gene_names$gene_name[match(genes_hits[[1]],gene_names$sheep_ID)]
    
    
    #######################################################
    dat_write=cbind(names_plot,genes_hits[[1]],genes_hits[[2]],genes_hits[[3]])
    names_plot=gene_names$gene_name[match(genes_hits[[1]],gene_names$sheep_ID)]
    colnames(dat_write)=c('human_gene_names','ensembl_ID',
                          paste(location1,types_headers1,c('FDR','logFC'),sep='_'))
    for(i in 1:nrow(dat_write)){
      if(is.na(dat_write[i,1])){
        if(dat_write[i,2]%in%orths$ens_IDs){
          dat_write[i,1]=orths$gene_names_new[orths$ens_IDs==dat_write[i,2]] 
        }
      }
    }
    dat_write=cbind(dat_write,dat_write[,2]%in%ISGs$Gene.ID)
    colnames(dat_write)[ncol(dat_write)]='is.ISG'
    mutual_DE=which(genes_hits[[2]]<0.05)
    write.csv(dat_write[mutual_DE,],
              file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE/',
                         location1,'.',types_list1,'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),row.names = F)
    write.csv(dat_write,
              file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',
                         location1,'.',types_list1,'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),row.names = F)
    
    #############
    names_plot=gene_names$gene_name[match(genes_hits[[4]],gene_names$sheep_ID)]
    
    dat_write=cbind(names_plot,genes_hits[[4]],genes_hits[[5]],genes_hits[[6]])
    
    colnames(dat_write)=c('human_gene_names','ensembl_ID',
                          paste(location2,types_headers2,c('FDR','logFC'),sep='_'))
    mutual_DE=which(genes_hits[[5]]<0.05)
    write.csv(dat_write[mutual_DE,],
              file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE/',
                         location2,'.',types_list2,'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),
              row.names = F)
    write.csv(dat_write,
              file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',
                         location2,'.',types_list2,'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),
              row.names = F)
    print(c(types_list1,types_list2,dpi_levels))
  }
}
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE'
setwd(dir)

logFC_threshold=3
gene_list=c()
files=list.files(dir,pattern='.csv')
files=c(files[!grepl('control',files)],files[grepl('control',files)])
for(i in files[!grepl('control',files)]){
  dat=read.csv(i,stringsAsFactors = F)
  gene_list=c(gene_list,dat$ensembl_ID[abs(as.numeric(dat[,4]))>logFC_threshold])
}
gene_list=unique(gene_list)
length(gene_list)

gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

location1='Sassari'
location2='Teramo'
dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all'
setwd(dir)
logFCs=matrix(ncol=length(files),nrow=length(unique(gene_list)))
for(i in 1:length(files)) {
  dat=read.csv(files[i],stringsAsFactors = F)
  logFCs[,i]=dat[match(gene_list,dat$ensembl_ID),4]
}

colnames(logFCs)=sapply(files,function(x) gsub('.dpi','.\ndpi',strsplit(x,'.csv')[[1]][1]))
rownames(logFCs)=gene_names$gene_name[match(gene_list,gene_names$sheep_ID)]
heatmap(logFCs, Colv = NA, Rowv = NA, scale="column",
       col= colorRampPalette(c('#354bed','#f2e292','#e42825'))(25))

##########################################################################
##########################################################################
##########################################################################
##########################################################################

logFC_threshold=2
FDR_threshold=0.05

dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE'
setwd(dir)



files=list.files(dir,pattern='.csv')
files=c(files[grepl('Sassari.control',files)],
         files[grepl('Teramo.control',files)],
         files[grepl('BTV8',files)],
            files[grepl('Sassari.2013',files)],
        files[grepl('Teramo.2013',files)],
        files[grepl('Sassari.2006',files)],
        files[grepl('Teramo.2006',files)]
       ) 
gene_list=c()
for(i in files[!grepl('control',files)]){
  dat=read.csv(i,stringsAsFactors = F)
  gene_list=c(gene_list,dat$ensembl_ID[abs(as.numeric(dat[,4]))>logFC_threshold&
                as.numeric(dat[,3])<FDR_threshold])
}
gene_list=unique(gene_list)
length(gene_list)
gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

location1='Sassari'
location2='Teramo'
dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all'
setwd(dir)

logFCs=FDRs=matrix(ncol=length(files),nrow=length(unique(gene_list)))
for(i in 1:length(files)) {
  dat=read.csv(files[i],stringsAsFactors = F)
  logFCs[,i]=dat[match(gene_list,dat$ensembl_ID),4]
  FDRs[,i]=dat[match(gene_list,dat$ensembl_ID),3]
}

colnames(logFCs)=sapply(files,function(x) gsub('.dpi','\ndpi',strsplit(x,'.csv')[[1]][1]))
# 
 colnames(logFCs)=sapply(files,function(x) gsub('.dpi','',strsplit(x,'.csv')[[1]][1]))
 colnames(logFCs)=colnames(FDRs)=sapply(colnames(logFCs),function(x) paste(substr((strsplit(x,'\\.')[[1]]),1,4),collapse='.'))

rownames(logFCs)=rownames(logFCs)=gene_names$gene_name[match(gene_list,gene_names$sheep_ID)]

library(gplots)
library(pheatmap)
library(biomaRt)


sheepgenes=useMart("ensembl",dataset="oaries_gene_ensembl")
test_orthologues=getBM(attributes=c("ensembl_gene_id",
                                     "hsapiens_homolog_associated_gene_name","hsapiens_homolog_orthology_confidence","hsapiens_homolog_ensembl_gene"),
                        filters = 'ensembl_gene_id', values = gene_list,mart=sheepgenes)

for(i in 1:length(gene_list)){
  if(is.na(rownames(logFCs)[i])){
    if(gene_list[i]%in%test_orthologues$ensembl_gene_id){
      trials=test_orthologues$hsapiens_homolog_associated_gene_name[test_orthologues$ensembl_gene_id==gene_list[i]]
      if(any(trials!='')){
        rownames(logFCs)[i]=trials[which(trials!='')[1]]
      }
    }
  }
}
extra_genes=read.csv('/home/oscar/Documents/orthologue_sets/Alex/Ov_fullmerge_final.csv')
for(i in 1:length(gene_list)){
  if(is.na(rownames(logFCs)[i])){
    if(gene_list[i]%in%extra_genes$Gene.ID_ov){
      trials=c(extra_genes$Gene_name[extra_genes$Gene.ID_ov==gene_list[i]],
               extra_genes$Gene_name_human[extra_genes$Gene.ID_ov==gene_list[i]])
      if(any(trials!='')){
        rownames(logFCs)[i]=trials[which(trials!='')[1]]
      }
    }
  }
}

write.csv(FDRs,file='~/FDR_table.Massimo_FDR0.05_logFC2.csv')
write.csv(logFCs,file='~/logFC_table.Massimo_FDR0.05_logFC2.csv')

rownames(logFCs)[which(gene_list=='ENSOARG00000006727')]
which(gene_list=='ENSOARG00000006727')
pheatmap(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],
         cluster_rows=T,cluster_cols=F,treeheight_row = 0, treeheight_col = 0)
pheatmap(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],
         cluster_rows=T,cluster_cols=F,treeheight_row = 0, treeheight_col = 0,
          color= colorRampPalette(rev(RColorBrewer::brewer.pal(7,'RdYlBu'))[c(1:3,5:7)])(60))
(gene_list)[is.na(rownames(logFCs))|rownames(logFCs)=='']
length(which(rownames(logFCs)==''|is.na(rownames(logFCs))))

#####
#####
#####
#####
#####
#####
#####
bigloves=c("ENSOARG00000010070","CD163", "CD274","CD28", "CD4","CD69","CD8A","CTLA4","FOXP3","IL2RA","NCAM1","NCR1","S1PR1","S1PR2","S1PR3","S1PR4")


setwd(dir)
files=list.files(dir,pattern='.csv')
logFCs=matrix(ncol=length(files),nrow=length(unique(gene_list)))
gene_list=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt')[,1]


gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)
all_genes_names=gene_names$gene_name[match(gene_list,gene_names$sheep_ID)]


sheepgenes=useMart("ensembl",dataset="oaries_gene_ensembl")
test_orthologues=getBM(attributes=c("ensembl_gene_id",
                                    "hsapiens_homolog_associated_gene_name","hsapiens_homolog_orthology_confidence","hsapiens_homolog_ensembl_gene"),
                       filters = 'ensembl_gene_id', values = gene_list,mart=sheepgenes)
test_orthologues$hsapiens_homolog_associated_gene_name
for(i in 1:length(gene_list)){
  if(is.na(all_genes_names[i])){
    if(gene_list[i]%in%test_orthologues$ensembl_gene_id){
      trials=test_orthologues$hsapiens_homolog_associated_gene_name[test_orthologues$ensembl_gene_id==gene_list[i]]
      if(any(trials!='')){
        all_genes_names[i]=trials[which(trials!='')[1]]
      }
    }
  }
}
extra_genes=read.csv('/home/oscar/Documents/orthologue_sets/Alex/Ov_fullmerge_final.csv')
for(i in 1:length(gene_list)){
  if(is.na(all_genes_names[i])){
    if(gene_list[i]%in%extra_genes$Gene.ID_ov){
      trials=c(extra_genes$Gene_name[extra_genes$Gene.ID_ov==gene_list[i]],
               extra_genes$Gene_name_human[extra_genes$Gene.ID_ov==gene_list[i]])
      if(any(trials!='')){
        all_genes_names[i]=trials[which(trials!='')[1]]
      }
    }
  }
}

bigloves_ens=gene_list[match(bigloves,all_genes_names)]


bigloves_ens[1]=bigloves[1]

write.csv(rbind(bigloves_ens,bigloves),file='~/Downloads/for_Vanessa.csv')
bigloves_ensna=bigloves_ens[!is.na(bigloves_ens)]

dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all'
setwd(dir)
logFCs=matrix(ncol=length(files),nrow=length(unique(bigloves_ensna)))
files=list.files(dir,pattern='.csv')
files=c(files[grepl('Sassari.control',files)],
        files[grepl('Teramo.control',files)],
        files[grepl('BTV8',files)],
        files[grepl('Sassari.2013',files)],
        files[grepl('Teramo.2013',files)],
        files[grepl('Sassari.2006',files)],
        files[grepl('Teramo.2006',files)]
) 
for(i in 1:length(files)) {
  dat=read.csv(files[i],stringsAsFactors = F)
  logFCs[,i]=dat[match(bigloves_ensna,dat$ensembl_ID),4]
}
rownames(logFCs)=bigloves[!is.na(bigloves_ens)]
colnames(logFCs)=sapply(files,function(x) gsub('.dpi','',strsplit(x,'.csv')[[1]][1]))
colnames(logFCs)=sapply(colnames(logFCs),function(x) paste(substr((strsplit(x,'\\.')[[1]]),1,4),collapse='.'))
pheatmap(logFCs,
         cluster_rows=F,cluster_cols=F,treeheight_row = 0, treeheight_col = 0)



###########
###########
# ###########
# ###########
# ###########
# heatmap(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),], Colv = NA, Rowv = NA, scale="column",
#         col= colorRampPalette(c('#354bed','#AAAAAA','#e42825'))(25))
# logFCs[3,]
# heatmap.2(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],density.info = "none",  scale="column",
#           col= colorRampPalette(c('#354bed','#AAAAAA','#e42825'))(25),trace="none")
# ?heatmap.2
# 
# library(RColorBrewer)
# 
# heatmap.2(logFCs[apply(logFCs,1,function(x) length(which(!is.na(x)))>15),],density.info = "none",  scale="column",
#           col= colorRampPalette(brewer.pal(6,'RdYlBu'))(25),trace="none")


