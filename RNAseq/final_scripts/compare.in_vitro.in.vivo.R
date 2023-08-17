
conds_invit=c('BTV.8-2017','BTV.1-2013','BTV.1-2006')
conds_inviv=c('SMI8','SMI13','SMI6')
dpi_comp=1
png(paste('/home/oscar/Pictures/plots/Sheep_megadata/11.5.22/RNA_Seq/invivo_vs_invitro_dpi'
          ,dpi_comp,'.png',sep=''),height=900,width=600)
par(mfrow=c(3,2),mar=c(4,4,1,1))

for(i in 1:length(conds_invit)){
  for(time in c('6h','12h')){
  strain=conds_invit[i]
  in_vitro=read.csv(paste('~/Documents/sheep_megadata/RNA.seq.lab.nov.21/GLM.results.merge.tab.all.batch.combatseq/',
                          conds_invit[i],'_',time,'.csv',sep=''))
  
  in_vivo=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/',
                         conds_inviv[i],'dpi',dpi_comp,'.csv',sep=''))
  
  
  in_vivo=in_vivo[in_vivo$X%in%in_vitro$X,]
  table(in_vivo$adj.P.Val<0.05)
  
  in_vitro=in_vitro[in_vitro$X%in%in_vivo$X,]
  in_vitro=in_vitro[match(in_vivo$X,in_vitro$X),]
  cols=rep(1,nrow(in_vitro))
  cols[in_vitro$FDR<0.05]=2
  cols[in_vivo$adj.P.Val<0.05]=3
  cols[in_vivo$adj.P.Val<0.05&in_vitro$FDR<0.05]=4
  
  in_vitro=in_vitro[order(cols),]
  in_vivo=in_vivo[order(cols),]
  cols=cols[order(cols)]
  cols=cols[order(cols)]
  col_nums=as.numeric(table(cols))
  colpal=RColorBrewer::brewer.pal(4,'Dark2')
  cols=scales::alpha(colpal[cols],cols/4)
  plot(in_vivo$logFC,in_vitro$logFC,col=cols,pch=19,
       main=paste(strain,' in vivo 3dpi vs in vitro',time))
  legend(x='topleft',legend=paste(
    c('notDE','DE in vitro only','DE in vivo only ','DE both'),'| n=',
    col_nums)
    ,col=colpal,pch=19,bty='n')

  }
}
dev.off()


pdf(paste('/home/oscar/Pictures/plots/Sheep_megadata/11.5.22/RNA_Seq/invivo_vs_invitro_dpi'
          ,dpi_comp,'.pdf',sep=''),height=9.00,width=6.00)
par(mfrow=c(3,2),mar=c(4,4,1,1))

for(i in 1:length(conds_invit)){
  for(time in c('6h','12h')){
    strain=conds_invit[i]
    in_vitro=read.csv(paste('~/Documents/sheep_megadata/RNA.seq.lab.nov.21/GLM.results.merge.tab.all.batch.combatseq/',
                            conds_invit[i],'_',time,'.csv',sep=''))
    
    in_vivo=read.csv(paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/',
                           conds_inviv[i],'dpi',dpi_comp,'.csv',sep=''))
    
    
    in_vivo=in_vivo[in_vivo$X%in%in_vitro$X,]
    table(in_vivo$adj.P.Val<0.05)
    
    in_vitro=in_vitro[in_vitro$X%in%in_vivo$X,]
    in_vitro=in_vitro[match(in_vivo$X,in_vitro$X),]
    cols=rep(1,nrow(in_vitro))
    cols[in_vitro$FDR<0.05]=2
    cols[in_vivo$adj.P.Val<0.05]=3
    cols[in_vivo$adj.P.Val<0.05&in_vitro$FDR<0.05]=4
    
    in_vitro=in_vitro[order(cols),]
    in_vivo=in_vivo[order(cols),]
    cols=cols[order(cols)]
    cols=cols[order(cols)]
    col_nums=as.numeric(table(cols))
    colpal=RColorBrewer::brewer.pal(4,'Dark2')
    cols=scales::alpha(colpal[cols],cols/4)
    plot(in_vivo$logFC,in_vitro$logFC,col=cols,pch=19,
         main=paste(strain,' in vivo 3dpi vs in vitro',time))
    legend(x='topleft',legend=paste(
      c('notDE','DE in vitro only','DE in vivo only ','DE both'),'| n=',
      col_nums)
      ,col=colpal,pch=19,bty='n')
    
  }
}
dev.off()