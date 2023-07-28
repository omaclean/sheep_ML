####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#################################       new
##
###
#####
#################################

library(scales)
par(mfrow=c(1,1))


ISG_list_ov=read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv",stringsAsFactors = F)

#out_dir="/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/all_genes.11.11.19"#no/ doesn't work
common_ISGS=unique(c(ISG_list_ov$Gene.ID))

axes_lim=c(-4.5,7)

location='Teramo'
types_loop1=c('C','I13','I6')
types_list1=c("control",'2013','2006')
dpi_levels=c('0','2')
par(mfrow=c(2,2),mar=c(3,4,1,0))
for(dat_i in 2:(length(types_loop1)-1)){
   for(dat_j in (dat_i+1):length(types_loop1)){
    types_headers1=paste('SL',c('C',types_loop1[dat_i],types_loop1[dat_j]),sep='')
    types_list1=c(types_list1[1],types_list1[dat_i],types_list1[dat_j])
      for(lev2 in c('1','3','7')){
        dpi_levels[2]=lev2
        dat1=read.csv(file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/3_tier_GLM/sep.all/',
                           location,'.',types_list1[2],'_',types_list1[3],'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''))
    print(lev2)
      # for(i in 1:2){
      #   for(j in (i+1):3){
          dat1=dat1[dat1[,2]%in%common_ISGS,]
          
          col=rep(1,nrow(dat1))
          #col[(dat1[,3]<0.05|dat2[,3]<0.05)]=1
          
          
          col[(dat1[,3]<0.05|dat1[,5]<0.05)]=2
          col[dat1[,7]<0.05&dat1[,8]<0&dat1[,3]<0.05]=3
          col[dat1[,7]<0.05&dat1[,8]>0&dat1[,5]<0.05]=4
          
          colnames(dat1)
          
          #sig_gene_names[1:length(which(col>1)),i]=dat$gene_name[col>1]
          #1            #2                        #3                        #4             #5
          legend_set=c(paste("not sig either N=",length(which(col==1))),
                        paste('sig in at least one condition N=',length(which(col==2))),
                             paste(types_list1[dat_i],"signif gt",types_list1[dat_j], ' N=',length(which(col==3))),
                       paste(types_list1[dat_j],"signif gt",types_list1[dat_i],' N=',length(which(col==4))))
          
          
          colpal=scales::alpha(c("#FFCB0A","#FF0606","#0B4FF9"),0.6)#adad85#e1cccc
          colpal=c(scales::alpha("#AAAAAA",0.2),scales::alpha(c("#FFCB0A","#FF0606","#0B4FF9"),0.6))#adad85#e1cccc
          ###d69d00
          #  colpal=alpha(c("#000000","#e63900","#0000b3","#66ccff","#cda7a7","#8A8A8A"),1)#8A8A8A
          #alpha(col,0.7)
          #    png(paste("/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/7.11.19/",condition_name,".png",sep=""),width = 700, height = 700)
          
          labs=rep('',nrow(dat1))
          if(any(col==3|col==4)){
            threshold=sort(abs(dat1[col==3|col==4,8]),decreasing=T)[
              min(c(length(which(col==3|col==4)),10))]
            labs[(col==3|col==4)&abs(dat1[,8])>=threshold]=dat1[(col==3|col==4)&
                  abs(dat1[,8])>=threshold,1]
          }
          
          plot(dat1[col==1,4],dat1[col==1,6],
               col=colpal[1],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
               main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
          par(new=T)
          plot(dat1[col==2,4],dat1[col==2,6],
               col=colpal[2],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
               main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
          par(new=T)
          plot(dat1[col==3,4],dat1[col==3,6],
               col=colpal[3],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
               main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
          par(new=T)
          plot(dat1[col==4,4],dat1[col==4,6],
               col=colpal[4],pch=19,
               xlab=paste(colnames(dat1)[4] ,"log2FC  vs Uninfected\n"),
               ylab=paste(colnames(dat1)[6] ,"log2FC  vs Uninfected"),
               ylim=axes_lim,xlim=axes_lim,
               main=paste(location, 'dpi', lev2),
               cex=1.3,frame=F,xaxt ='n',yaxt ='n')
          
          text(dat1[,4]+rnorm(nrow(dat1),sd=.3),dat1[,6]+rnorm(nrow(dat1),sd=.3),labels=labs)
          axis(side=1,labels=c(-2:5)*2.5,at=c(-2:5)*2.5)
          axis(side=2,labels=c(-2:5)*2.5,at=c(-2:5)*2.5)
          abline(0,1,lty=2)
          legend(col=colpal,legend=legend_set,pch=19,x="bottomleft",bty="n")
          #  legend(col=alpha(colpal,0.9),legend=paste(legend_set,rep("( N=",length(legend_set)),sapply(1:length(legend_set),
          #      function(x) length(which(col==x))),rep(")",length(legend_set)),sep=""),
          #       lwd=15,x="topleft",bty="n")
          #  dev.off()
        }
    plot(1)
    #   }
    # }
  }
}


location='Sassari'
types_loop1=c('C','I13','I6','I8')
types_list1=c("control",'2013','2006','BTV8')
dpi_levels=c('0','2') #initiator - wrong deliberary
par(mfrow=c(3,3),mar=c(3,4,1,0))
for(dat_i in 2:(length(types_loop1)-1)){
  for(dat_j in (dat_i+1):length(types_loop1)){
    types_headers1=paste('SM',c('C',types_loop1[dat_i],types_loop1[dat_j]),sep='')
    
    for(lev2 in c('1','3','7')){
      dpi_levels[2]=lev2
      dat1=read.csv(file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/3_tier_GLM/sep.all/',
                               location,'.',types_list1[dat_i],'_',types_list1[dat_j],'.dpi.',dpi_levels[1],'.vs.',dpi_levels[2],'.csv',sep=''))
      print(lev2)
      # for(i in 1:2){
      #   for(j in (i+1):3){
      dat1=dat1[dat1[,2]%in%common_ISGS,]
      
      col=rep(1,nrow(dat1))
      #col[(dat1[,3]<0.05|dat2[,3]<0.05)]=1
      
      col[(dat1[,3]<0.05|dat1[,5]<0.05)]=2
      col[dat1[,7]<0.05&dat1[,8]<0&dat1[,3]<0.05]=3
      col[dat1[,7]<0.05&dat1[,8]>0&dat1[,5]<0.05]=4
      
      #sig_gene_names[1:length(which(col>1)),i]=dat$gene_name[col>1]
      #1            #2                        #3                        #4             #5
      legend_set=c(paste("not sig either N=",length(which(col==1))),
                   paste('sig in at least one condition N=',length(which(col==2))),
                   paste(types_list1[dat_i],"signif gt",types_list1[dat_j], ' N=',length(which(col==3))),
                   paste(types_list1[dat_j],"signif gt",types_list1[dat_i],' N=',length(which(col==4))))
      
      
      colpal=scales::alpha(c("#FFCB0A","#FF0606","#0B4FF9"),0.6)#adad85#e1cccc
      colpal=c(scales::alpha("#AAAAAA",0.2),scales::alpha(c("#FFCB0A","#FF0606","#0B4FF9"),0.6))#adad85#e1cccc
      ###d69d00
      #  colpal=alpha(c("#000000","#e63900","#0000b3","#66ccff","#cda7a7","#8A8A8A"),1)#8A8A8A
      #alpha(col,0.7)
      #    png(paste("/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/7.11.19/",condition_name,".png",sep=""),width = 700, height = 700)
      labs=rep('',nrow(dat1))
      if(any(col==3|col==4)){
        threshold=sort(abs(dat1[col==3|col==4,8]),decreasing=T)[
          min(c(length(which(col==3|col==4)),10))]
        labs[(col==3|col==4)&abs(dat1[,8])>=threshold]=dat1[(col==3|col==4)&abs(dat1[,8])>=threshold,1]
      }
      
      plot(dat1[col==1,4],dat1[col==1,6],
           col=colpal[1],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
           main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
      par(new=T)
      plot(dat1[col==2,4],dat1[col==2,6],
           col=colpal[2],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
           main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
      par(new=T)
      plot(dat1[col==3,4],dat1[col==3,6],
           col=colpal[3],pch=19,xlab="",ylab="",ylim=axes_lim,xlim=axes_lim,
           main="",cex=1.3,frame=F,xaxt ='n',yaxt ='n')
      par(new=T)
      plot(dat1[col==4,4],dat1[col==4,6],
           col=colpal[4],pch=19,
           xlab=paste(colnames(dat1)[4] ,"log2FC  vs Uninfected\n"),
           ylab=paste('\n\n',colnames(dat1)[6] ,"log2FC  vs Uninfected"),
           ylim=axes_lim,xlim=axes_lim,
           main=paste(location,'dpi', lev2),
           cex=1.3,frame=F,xaxt ='n',yaxt ='n')
      
      text(dat1[,4]+rnorm(nrow(dat1),sd=.3),dat1[,6]+rnorm(nrow(dat1),sd=.3),labels=labs)
      axis(side=1,labels=c(-2:5)*2.5,at=c(-2:5)*2.5)
      axis(side=2,labels=c(-2:5)*2.5,at=c(-2:5)*2.5)
      abline(0,1,lty=2)
      legend(col=colpal,legend=legend_set,pch=19,x="bottomleft",bty="n")
      #  legend(col=alpha(colpal,0.9),legend=paste(legend_set,rep("( N=",length(legend_set)),sapply(1:length(legend_set),
      #      function(x) length(which(col==x))),rep(")",length(legend_set)),sep=""),
      #       lwd=15,x="topleft",bty="n")
      #  dev.off()
    }
    #plot(1)
    #   }
    # }
  }
}



# 
# 