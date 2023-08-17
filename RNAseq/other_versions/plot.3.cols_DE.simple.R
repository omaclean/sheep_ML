####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#################################       new
##
###
#####
#################################
virus='2006'
location='Sassari'
library(scales)
par(mfrow=c(1,1))
T2006d1=read.csv(paste("~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/",location,".",
  virus,".dpi.0.vs.1.csv",sep=''),stringsAsFactors = F)
T2006d3=read.csv(paste("~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/",location,".",
                              virus,".dpi.0.vs.3.csv",sep=''),stringsAsFactors = F)
T2006d7=read.csv(paste("~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/",location,".",
                               virus,".dpi.0.vs.7.csv",sep=''),stringsAsFactors = F)


ISG_list_ov=read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv",stringsAsFactors = F)

#out_dir="/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/all_genes.11.11.19"#no/ doesn't work
common_ISGS=unique(c(ISG_list_ov$Gene.ID))
data=list(T2006d1,T2006d3,T2006d7)
names(data)=paste(location,virus,'dpi',c(1,3,7),sep='_')

sig_gene_names=matrix(nrow=10000,ncol=4,dat=NA)
to_name=c("TMEM140","CDADC1","RICTOR","CD274","SLC25A30", "ELMO2","FCGR3A","SP140","SASS6","CCDC120" ,
          "MED17","ISG20" ,"TIFA","ZNF135","MED25", "IDO1"  )
axes_lim=c(-4.5,7)
par(mfrow=c(2,2),mar=c(4,4,1,1))
for(i in 1:2){
  for(j in (i+1):3){
  dat1=data[[i]]
  dat2=data[[j]]
  dat2=dat2[dat2[,2]%in%common_ISGS,]
  dat1=dat1[dat1[,2]%in%dat2[,2],]
  dat2=dat2[dat2[,2]%in%dat1[,2],]
  

  col=rep(1,nrow(dat1))
  #col[(dat1[,3]<0.05|dat2[,3]<0.05)]=1
  
  
  col[(dat2[,4]<dat1[,4]&dat1[,3]<0.05)]=2
  col[(dat2[,4]>dat1[,4]&dat2[,3]<0.05)]=3
  
  
  
  
  dat1=dat1[!is.na(col),]
  dat2=dat2[!is.na(col),]
  col=col[!is.na(col)]
  
  #sig_gene_names[1:length(which(col>1)),i]=dat$gene_name[col>1]
  #1            #2                        #3                        #4             #5
  legend_set=c("not sig either",paste(names(data)[i],"gt",names(data)[j], ' N=',length(which(col==2)))
               ,paste(names(data)[i],"gt",names(data)[j],' N=',length(which(col==3))))
  
  
  print(range(c(dat1[,4],dat2[,4])))
  colpal=scales::alpha(c("#FFCB0A","#FF0606","#0B4FF9"),0.6)#adad85#e1cccc
  colpal=scales::alpha(c("#AAAAAA","#FF0606","#0B4FF9"),0.6)#adad85#e1cccc
  ###d69d00
  #  colpal=alpha(c("#000000","#e63900","#0000b3","#66ccff","#cda7a7","#8A8A8A"),1)#8A8A8A
  #alpha(col,0.7)
  #    png(paste("/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/7.11.19/",condition_name,".png",sep=""),width = 700, height = 700)
  
  labs=rep('',nrow(dat))
  for(names in to_name){
    labs[which(dat1[,1]==names)]=names
  }
  
  plot(dat1[col==1,4],dat2[col==1,4],
       col=colpal[1],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
       main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
  par(new=T)
  plot(dat1[col==2,4],dat2[col==2,4],
       col=colpal[2],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
       main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
  par(new=T)
  plot(dat1[col==3,4],dat2[col==3,4],
       col=colpal[3],pch=19,xlab=paste(names(data)[i],"log2FC  vs Uninfected"),
       ylab=paste(names(data)[j],"log2FC  vs Uninfected"),ylim=axes_lim,xlim=axes_lim,
       cex=1.3,frame=F,xaxt ='n',yaxt ='n')
  text(dat1[,4],dat2[,4],labels=labs)
  axis(side=1,labels=c(-2:5)*2.5,at=c(-2:5)*2.5)
  axis(side=2,labels=c(-2:5)*2.5,at=c(-2:5)*2.5)
  abline(0,1,lty=2)
  legend(col=scales::alpha(colpal,0.9),legend=legend_set,pch=19,x="bottomleft",bty="n")
  #  legend(col=alpha(colpal,0.9),legend=paste(legend_set,rep("( N=",length(legend_set)),sapply(1:length(legend_set),
  #      function(x) length(which(col==x))),rep(")",length(legend_set)),sep=""),
  #       lwd=15,x="topleft",bty="n")
  #  dev.off()
  
  }
}
plot(1)
# 
# 