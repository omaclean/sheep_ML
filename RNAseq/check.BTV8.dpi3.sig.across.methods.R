
ISG_list_ov=read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv",stringsAsFactors = F)

setwd('~/Documents/sheep_megadata/RNA_Seq_1.12.20/')

for(dpi in c('1','3','7')){
  if(length(readLines(paste('Quan.exact.test/SMI8dpi0-VS-SMI8dpi',dpi,'_siggenes_edgeR_final.csv',sep='')))>0){
    
  
  exact=read.csv(paste('Quan.exact.test/SMI8dpi0-VS-SMI8dpi',dpi,'_siggenes_edgeR_final.csv',sep=''),header=F)
  }else{exact=NULL}
  dat2=read.csv(paste('mixed_effects/controls_infect_dpi0_other/SMI8dpi',dpi,'.csv',sep=''),header=T)
  dream_no_cont=read.csv(paste('mixed_effects/no_controls_infect_dpi0_other/SMI8dpi',dpi,'.csv',sep=''))
  fixed=read.csv(paste('fixed.effects/sep.all/Sassari.BTV8.dpi.0.vs.',dpi,'.csv',sep=''))
  print(c(length(which(exact[,1]%in%ISG_list_ov$Gene.ID&exact[,3]<0.05)),
          length(which(dat2$X%in%ISG_list_ov$Gene.ID&dat2$adj.P.Val<0.05)),
          length(which(dream_no_cont$X%in%ISG_list_ov$Gene.ID&dream_no_cont$adj.P.Val<0.05)),
          length(which(fixed$ensembl_ID%in%ISG_list_ov$Gene.ID&fixed$Sassari_SMI8_FDR<0.05))
            ))
  print(c(length(which(exact[,3]<0.05)),
          length(which(dat2$adj.P.Val<0.05)),
          length(which(dream_no_cont$adj.P.Val<0.05)),
          length(which(fixed$Sassari_SMI8_FDR<0.05))
  ))
  print(c('###############dpi',dpi))
}
str(dat)
