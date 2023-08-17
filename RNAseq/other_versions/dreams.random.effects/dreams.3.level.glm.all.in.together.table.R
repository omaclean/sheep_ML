library('variancePartition')
library('BiocParallel')
library('edgeR')

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

males=read.csv('/home/oscar/Documents/sheep_megadata/temperature_humidity_gender.csv')
males=males$animal.number[grepl('SMI6',males$animal.number)&
                            males$gender=='m']
counts=counts[,!grepl(paste(males,collapse='_'),colnames(counts))]

genes_hits=list()
dpi_levels=c('0',NA)
location1='Teramo'
colnames(counts)

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
      rownames(counts2)=counts[,1]
      counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
      counts2=counts2[rowSums(counts2>0)>7,]
      
      dpis=sapply(colnames(counts2),function(x) 
        tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
      
      
      counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]
      
      dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]
      
      dpis=factor(dpis,levels=dpi_levels)
      animal_names=as.character(sapply(colnames(counts2),function(x) 
        paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
      
      animal_names=as.factor(animal_names)
      inf_status=! (grepl(types_headers1[1],animal_names)| dpis==0)
      strain=grepl(types_headers1[3],animal_names)&inf_status
      ################################################
      
      
      geneExpr = DGEList( counts2 )  
      
      geneExpr = calcNormFactors( geneExpr )
      
      meta_data=data.frame(inf_status,strain,animal_names)
      rownames(meta_data)=colnames(counts2)
      form <- ~ inf_status + strain+ (1|animal_names) 
      
      
      # estimate weights using linear mixed model of dream
      vobjDream = voomWithDreamWeights( geneExpr, form, meta_data )
      fitmm = dream( vobjDream, form, meta_data )
      
      hits=topTable(fitmm,coef='inf_statusTRUE',number=nrow(counts2),sort.by = 'n')
      hits2=topTable(fitmm,coef=c(2,3),number=nrow(counts2),sort.by = 'n')
      hits2=hits2[match(rownames(hits),rownames(hits2)),]
      hits3=topTable(fitmm,coef='strainTRUE',number=nrow(counts2),sort.by = 'n')
      hits3=hits3[match(rownames(hits),rownames(hits3)),]
      genes_hits[[1]]=rownames(hits)
      genes_hits[[2]]=hits$adj.P.Val
      genes_hits[[3]]=hits$logFC
      ##################
      
      
      genes_hits[[4]]=rownames(hits2)
      genes_hits[[5]]=hits2$adj.P.Val
      genes_hits[[6]]=hits2[,1]+hits2[,2]
      
      
      genes_hits[[7]]=rownames(hits3)
      genes_hits[[8]]=hits3$adj.P.Val
      genes_hits[[9]]=hits3$logFC
      for(test_I in 1:nrow(hits3)){
        if(rownames(hits3)[test_I]!=rownames(hits)[test_I]|rownames(hits2)[test_I]!=rownames(hits)[test_I]){
          print('errro')
        }
      }
      
      gene_names=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
      names_plot=  gene_names$gene_names_new[match(rownames(hits),gene_names$ens_IDs)]
      dat_write=cbind(names_plot,rownames(genes_hits),genes_hits[[1]],genes_hits[[2]],genes_hits[[3]]
                      ,genes_hits[[5]],genes_hits[[6]],genes_hits[[8]],genes_hits[[9]])
      
      
      
      
      colnames(dat_write)=c('human_gene_names','ensembl_ID',
                            paste(location1,types_headers1[2],c('FDR','logFC'),sep='_'),
                            paste(location1,types_headers1[3],c('FDR','logFC'),sep='_'),
                            paste(location1,types_headers1[2],'vs',types_headers1[3],c('FDR','logFC'),sep='_'))
      write.csv(dat_write,
                file=paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/3.tier.GLM/',
                           location1,'.',types_list1[2],'_',types_list1[3],'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''),row.names = F)
    }
  }
}