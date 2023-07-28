
library(dplyr);library(ggplot2);library(reshape2)
mat2=read.csv('/home/oscar/Documents/sheep_megadata/30.6.21/dat_plots_filt/BTMs.final.csv',row.names = 1)

conds=c('SLC','SMC','SMI8','SLI6','SMI6','SLI13','SMI13')
dpis=c('0','1','3','7')


BTMs_to_plot=c('RIG')
ggplotme=function(data,condition_name,ylim1,class1){
  ggplot(data,aes(y=log2(read_prop+1),fill=condition,x=ORF_condition))+#geom_linerange( aes(ymin=min,ymax=max,col=species,fill=species,alpha=0.8))+
 #   geom_abline(slope=0,intercept=log2(12.98+1) ,colour="#000000")+
    theme_bw()+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
    geom_crossbar(aes(y=log2(mean+1),ymin=log2(mean+1),ymax=log2(mean+1),col=condition,fill=condition,alpha=1))+
    geom_dotplot(binaxis = "y",stackdir = "center",dotsize=2,stackratio=0.8,binwidth=ylim1/20,alpha=.7)+
    ylab("log2(read_prop+1)")+
    ggtitle(paste(class1,condition_name))+xlab(paste(class1,"gene and infection state"))+ylim(0,log2(ylim+1))
  
  #ggplot(data,aes(y=log(cpm+1),fill=species,x=condition_casual))+geom_dotplot(aes(alpha=0.8),binaxis = "y", 
  #                                                                           stackdir = "center",dotsize=0.5,stackratio=0.4)+
  #  ggtitle(paste("Interferon alpha",condition_name))+xlab("Interferon copy and infection state")+theme_bw()+ylim(0,4)
}


  for(i in 1:4){
    counts2=mat2[grep(paste(BTMs_to_plot,collapse='|'),rownames(mat2)),]
    data=melt(counts2,factorsAsStrings=F)
    dpi_i=dpis[i]
  
    #data,condition_name,ylim1,class1
    data$dpi=sapply(as.character(data[,1]),function(x) tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
    data$condition=sapply(as.character(data[,1]),function(x) strsplit(x,'\\.')[[1]][1])
    data[,1]=rep('comb_BTM(RIG1)',nrow(data))
    colnames(data)=c('gene','read_prop','dpi','condition')
    data=data[unlist(sapply(conds,function(x)which(data$condition==x))),]
    data$ORF_condition=factor(paste(data$gene,"\n",data$condition),
                              levels=unique(paste(data$gene,"\n",data$condition)))
    ylim=max(data$read_prop)
    data=data[data$dpi==dpi_i,]
    library(psych)
    data$mean=sapply(1:nrow(data),function(x) 
      (-1)+geometric.mean(1+data$read_prop[data$ORF_condition==data$ORF_condition[x]]))
    data$max=sapply(1:nrow(data),function(x) 
      max(data$read_prop[data$ORF_condition==data$ORF_condition[x]]))
    data$min=sapply(1:nrow(data),function(x) 
      min(data$read_prop[data$ORF_condition==data$ORF_condition[x]]))
 
  class1='test'
  condition_name=paste('dpi',dpis[i])
  if(i==1){
    p1=ggplotme(data,condition_name,ylim,class1)
  }
  if(i==2){
    p2=ggplotme(data,condition_name,ylim,class1)
  }
  if(i==3){
    p3=ggplotme(data,condition_name,ylim,class1)
  }
  if(i==4){
    p4=ggplotme(data,condition_name,ylim,class1)
  }

}

grid.arrange(p1,p2,p3,p4,nrow=2)
