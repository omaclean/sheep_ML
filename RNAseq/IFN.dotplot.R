###iFN alpha
library(dplyr);library(ggplot2);library(reshape2)
orthologues=read.csv(file="/home/oscar/Documents/orthologue_sets/sheep_and_cow_ID_only.csv",stringsAsFactors = F)

sheep_IFNA_names=c("ENSOARG00000008723", "ENSOARG00000008754", "ENSOARG00000014373","ENSOARG00000018611",
                   "ENSOARG00000010517")

#### ifn beta
sheep_IFNB_names=c("ENSOARG00000008675","ENSOARG00000008791","ENSOARG00000008697")



 IFNs=list()
# IFNs[['IFNa']]=sheep_IFNA_names
# IFNs[['IFNb']]=sheep_IFNB_names

IFNs[['Van_genes']]=c("ENSOARG00000006292" ,"ENSOARG00000010278")
########################
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T,stringsAsFactors = F)

rownames(counts)=counts$Geneid

counts2=apply(counts[,2:ncol(counts)],2,function(x) (10^6)*x/sum(x))

length(as.numeric(as.matrix(counts[,2:ncol(counts)])))

###########


counts2=counts2[rownames(counts2)%in%c(unlist(IFNs)),]
length(which(rownames(counts2)%in%c(unlist(IFNs))))
which(counts2>0,arr.ind = T)
dim(counts2)

#test
library(tidyverse)


conditions=c('SMC','SLC','SMI13','SLI13','SMI6','SLI6','SMI8')

# A) # southern lights# aurora #



ggplotme=function(data,condition_name,ylim1,class1){
  ggplot(data,aes(y=log2(read_prop+1),fill=condition,x=ORF_condition))+#geom_linerange( aes(ymin=min,ymax=max,col=species,fill=species,alpha=0.8))+
    geom_abline(slope=0,intercept=log2(12.98+1) ,colour="#000000")+
    theme_bw()+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
    geom_crossbar(aes(y=log2(mean+1),ymin=log2(mean+1),ymax=log2(mean+1),col=condition,fill=condition,alpha=1))+
    geom_dotplot(binaxis = "y",stackdir = "center",dotsize=1.4,stackratio=0.6,binwidth=ylim1/50)+ylab("log2(read_prop+1)")+
    ggtitle(paste(class1,condition_name))+xlab(paste(class1,"gene and infection state"))+ylim(0,ylim)
  
  #ggplot(data,aes(y=log(cpm+1),fill=species,x=condition_casual))+geom_dotplot(aes(alpha=0.8),binaxis = "y", 
  #                                                                           stackdir = "center",dotsize=0.5,stackratio=0.4)+
  #  ggtitle(paste("Interferon alpha",condition_name))+xlab("Interferon copy and infection state")+theme_bw()+ylim(0,4)
}

#gene_names$casual=sapply(1:nrow(gene_names),function(x) paste(gene_names$GeneName[x],'_O',1+
#                                                                x-which(gene_names$GeneName==gene_names$GeneName[x])[1] ,sep=''))



range(as.numeric(unlist(c(counts[counts$Geneid%in%c(unlist(IFNs)),2:ncol(counts)]))))



plot(1,1)
dpi_i=c(0,1,3,7)
class1='test'
for(condition_name in names(IFNs)){

  for(i in 1:4){
    
    data=melt(counts2,factorsAsStrings=F)
    data[,1]=as.character(data$Var1)
    data$dpi=sapply(as.character(data[,2]),function(x) tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))
    data[,2]=sapply(as.character(data[,2]),function(x) strsplit(x,'\\.')[[1]][1])
    colnames(data)=c('gene','condition','read_prop','dpi')
    data=data[data$gene%in%IFNs[[condition_name]],]
   
  
  #rownames(counts)
  
    
    # order=rep(0,nrow(data))
    # ordercount=1
    # for(condit_i in conditions){
    #   
    #   hits=which(data$condition==condit_i)
    #   
    #   hits=hits[unlist(sapply(gene_names$casual,function(x) which(data$gene_casual[hits]==x)))]
    #   order[ordercount:(ordercount+length(hits)-1)]=hits
    #   ordercount=ordercount+length(hits)
    # }
    #data=data[order,]
    data$ORF_condition=factor(paste(data$gene,"\n",data$condition),
                              levels=unique(paste(data$gene,"\n",data$condition)))
    
    
    library(psych)
    data$mean=sapply(1:nrow(data),function(x) 
      (-1)+geometric.mean(1+data$read_prop[data$ORF_condition==data$ORF_condition[x]]))
    data$max=sapply(1:nrow(data),function(x) 
      max(data$read_prop[data$ORF_condition==data$ORF_condition[x]]))
    data$min=sapply(1:nrow(data),function(x) 
      min(data$read_prop[data$ORF_condition==data$ORF_condition[x]]))
    ylim=max(data$read_prop)
    
    
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
}
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2)
counts2
