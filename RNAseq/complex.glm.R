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
for(lev2 in c('1','3','7')){
  dpi_levels[2]=lev2
  for(i in 2:(length(types_loop1)-1)){
    for(j in (i+1):length(types_loop1)){
      types_headers1=paste('SM',c('C',types_loop1[i],types_loop1[j]),sep='')
      types_list1=c(types_list1_loop[1],types_list1_loop[i],types_list1_loop[j])
      
      ###########################
      
      gene_names=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
      
      counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
      rownames(counts)=counts$Geneid
      ###########
      types_list=types_list1
      types_headers=c(types_headers1)
      
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
      infected=! (grepl(types_headers1[1],animal_names)| dpis==0)
      strain2=grepl(types_headers1[3],animal_names)&infected
      ################################################
      design <- model.matrix(~infected+strain2+animal_names2)
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
      test= topTags(de2, n=nrow(counts2),sort.by = 'n')
      FDR_threshold=0.05
      genes_hits[[1]]=rownames(test$table)
      genes_hits[[2]]=test$table$FDR
      genes_hits[[3]]=test$table$logFC
      ##################
      
      contrast=rep(0,ncol(design))
      contrast[2]=1
      contrast[3]=1
      de2=glmQLFTest(d3,contrast=contrast)
      test= topTags(de2, n=nrow(counts2),sort.by = 'n')
      FDR_threshold=0.05
      genes_hits[[4]]=rownames(test$table)
      genes_hits[[5]]=test$table$FDR
      genes_hits[[6]]=test$table$logFC
      ####################
      contrast=rep(0,ncol(design))
      
      contrast[3]=1
      de2=glmQLFTest(d3,contrast=contrast)
      test= topTags(de2, n=nrow(counts2),sort.by = 'n')
      FDR_threshold=0.05
      genes_hits[[7]]=rownames(test$table)
      genes_hits[[8]]=test$table$FDR
      genes_hits[[9]]=test$table$logFC
      
      
      #############
      
      #######################################################
      names_plot=gene_names$gene_names_new[match(genes_hits[[1]],gene_names$ens_IDs)]
      
      for(test_i in 1:length(genes_hits[[1]])){
        if(genes_hits[[1]][i]!=genes_hits[[4]][i]){
          print('FAILURE')
        }
      }
      
      dat_write=cbind(names_plot,rownames(genes_hits),genes_hits[[1]],genes_hits[[2]],genes_hits[[3]]
                      ,genes_hits[[5]],genes_hits[[6]],genes_hits[[8]],genes_hits[[9]])
      
      colnames(dat_write)=c('human_gene_names','ensembl_ID',
                            paste(location1,types_headers1[2],c('FDR','logFC'),sep='_'),
                            paste(location1,types_headers1[3],c('FDR','logFC'),sep='_'),
                            paste(location1,types_headers1[2],'vs',types_headers1[3],c('FDR','logFC'),sep='_'))
      write.csv(dat_write,
                file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/3_tier_GLM/sep.all/',
                           location1,'.',types_list1[2],'_',types_list1[3],'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),row.names = F)
    }
  }
  
}

rmultinom(1,prob=rep(3*10^-6,3))

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

ISGs=read.csv('/home/oscar/Documents/Alex_RNA_seq/Ov_ISGlist.csv')
orths=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')

genes_hits=list()
dpi_levels=c('0',NA)
location1='Teramo'

types_loop1=c('C','I13','I6')
types_list1_loop=c("control",'2013','2006')
for(lev2 in c('1','3','7')){
  dpi_levels[2]=lev2
  for(i in 2:(length(types_loop1)-1)){
    for(j in (i+1):length(types_loop1)){
      types_headers1=paste('SL',c('C',types_loop1[i],types_loop1[j]),sep='')
      types_list1=c(types_list1_loop[1],types_list1_loop[i],types_list1_loop[j])
      
      ###########################
      
      gene_names=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
      
      counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
      rownames(counts)=counts$Geneid
      ###########
      types_list=types_list1
      types_headers=c(types_headers1)
      
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
      infected=! (grepl(types_headers1[1],animal_names)| dpis==0)
      strain2=grepl(types_headers1[3],animal_names)&infected
      ################################################
      design <- model.matrix(~infected+strain2+animal_names2)
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
      test= topTags(de2, n=nrow(counts2),sort.by = 'n')
      FDR_threshold=0.05
      genes_hits[[1]]=rownames(test$table)
      genes_hits[[2]]=test$table$FDR
      genes_hits[[3]]=test$table$logFC
      ##################
      
      contrast=rep(0,ncol(design))
      contrast[2]=1
      contrast[3]=1
      de2=glmQLFTest(d3,contrast=contrast)
      test= topTags(de2, n=nrow(counts2),sort.by = 'n')
      FDR_threshold=0.05
      genes_hits[[4]]=rownames(test$table)
      genes_hits[[5]]=test$table$FDR
      genes_hits[[6]]=test$table$logFC
      ####################
      contrast=rep(0,ncol(design))
      
      contrast[3]=1
      de2=glmQLFTest(d3,contrast=contrast)
      test= topTags(de2, n=nrow(counts2),sort.by = 'n')
      FDR_threshold=0.05
      genes_hits[[7]]=rownames(test$table)
      genes_hits[[8]]=test$table$FDR
      genes_hits[[9]]=test$table$logFC
      
      
      #############
      
      #######################################################
      names_plot=gene_names$gene_names_new[match(genes_hits[[1]],gene_names$ens_IDs)]
      
      for(test_i in 1:length(genes_hits[[1]])){
        if(genes_hits[[1]][i]!=genes_hits[[4]][i]){
          print('FAILURE')
        }
      }
      
      dat_write=cbind(names_plot,rownames(genes_hits),genes_hits[[1]],genes_hits[[2]],genes_hits[[3]]
                      ,genes_hits[[5]],genes_hits[[6]],genes_hits[[8]],genes_hits[[9]])
      
      colnames(dat_write)=c('human_gene_names','ensembl_ID',
                            paste(location1,types_headers1[2],c('FDR','logFC'),sep='_'),
                            paste(location1,types_headers1[3],c('FDR','logFC'),sep='_'),
                            paste(location1,types_headers1[2],'vs',types_headers1[3],c('FDR','logFC'),sep='_'))
      write.csv(dat_write,
                file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/3_tier_GLM/sep.all/',
                           location1,'.',types_list1[2],'_',types_list1[3],'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),row.names = F)
    }
  }
  
}

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


