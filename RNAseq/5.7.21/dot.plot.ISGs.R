##################################################################################################################################
#################################       new
##
###
#####
#################################

library(scales);library(RColorBrewer);library(ggplot2);library(gridExtra)
par(mfrow=c(1,1))


ISG_list_ov=read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv",stringsAsFactors = F)
dpi_levels=c('0',NA)
#out_dir="/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/all_genes.11.11.19"#no/ doesn't work
common_ISGS=unique(c(ISG_list_ov$Gene.ID))
# 
# axes_lim=c(-4.5,7)
# plot_nonsig=F
# locations=c('Teramo','Sassari')
# types_loop1=c('C','I13','I6')
# viruses=c('2013','2006',"BTV8")
# 
# plots=list()
# plot_i=0
# for(lev2 in c('1','3','7')){
#   plot_i=plot_i+1
#   plot_dat=matrix(nrow=0,ncol=4)
#   for(virus in viruses){
#     for(location in locations){
#       if(file.exists(paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',
#                            location,'.',virus,'.dpi.0.vs.',lev2,'.csv',sep=''))){
#         print('yay')
#         dpi_levels[2]=lev2
#         dat1=read.csv(file=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/',
#                                  location,'.',virus,'.dpi.0.vs.',lev2,'.csv',sep=''))
#         
#         
#         if(!plot_nonsig){
#           dat1=dat1[dat1[,3]<0.05,]
#         }
#         # for(i in 1:2){
#         #   for(j in (i+1):3){
#         dat1=dat1[dat1[,2]%in%common_ISGS,]
#         plot_dat=rbind(plot_dat,as.matrix(cbind(dat1[,4],dat1[,3]<0.05,
#                                                 rep(paste(location,'.',virus,sep=''),nrow(dat1)),
#                                                 rep(paste(location,'.',virus,'\n(N=',
#                                                           length(which(dat1[,3]<0.05)),')',sep=''),nrow(dat1))
#         )))
#         print(paste(location,'.',virus,lev2,nrow(dat1)))
#       }
#     }
#   }
#   
#   plot_dat=rbind(plot_dat,
#                  as.matrix(cbind(rnorm(400),rep(T,400),
#                                  rep('randomcrap',400), rep('randomcrap',400))))
#   plot_dat_save=plot_dat
#   plot_dat=as.data.frame(plot_dat)
#   colnames(plot_dat)=c('logFC','significant','condition','condition_lab')
#   plot_dat$colour=plot_dat$condition
#   if(plot_nonsig){
#     plot_dat$colour[!as.logical(plot_dat$significant)]='nonsig'
#     plot_dat$colour=factor(plot_dat$colour,levels=c(unique(plot_dat$condition),'nonsig'))
#     plot_dat$condition=factor(plot_dat$condition,levels=unique(plot_dat$condition))
#     colpal=c(scales::alpha(RColorBrewer::brewer.pal(length(unique(plot_dat$condition)),'Dark2'),.5),
#              scales::alpha('#AAAAAA',.3))
#   }else{
#     
#     plot_dat$colour=factor(plot_dat$colour,levels=unique(plot_dat$condition))
#     plot_dat$condition=factor(plot_dat$condition,levels=unique(plot_dat$condition))
#     colpal=scales::alpha(RColorBrewer::brewer.pal(length(unique(plot_dat$condition)),'Dark2'),.5)
#   }
#   plot_dat$logFC=as.numeric(plot_dat$logFC)
#   # plots[[plot_i]]= ggplot(plot_dat,aes(x=condition_lab,y=logFC,color=colour))+
#   #   geom_jitter()+scale_colour_manual(values=colpal)+theme_bw()+ylim(-5,10)
#   plots[[plot_i]]= ggplot(plot_dat,aes(x=condition_lab,y=logFC,fill=colour))+
#     geom_violin(scale='count')+scale_fill_manual(values=colpal)+
#     theme_bw()+ylim(-5,10)+ggtitle(paste('dpi',lev2))
#   
# }
# 
# grid.arrange(plots[[1]],plots[[2]],plots[[3]],ncol=1)

##############################
plot_nonsig=T

plots=list()
plot_i=0
translator_loc=list()
translator_loc[['Sassari']]='SMI'
translator_loc[['Teramo']]='SLI'
viruses=c('8','6')
virus_names=c('BTV8','2006')
for(lev2 in c('1','3')){
  plot_i=plot_i+1
  plot_dat=matrix(nrow=0,ncol=4)
  for(combo in 1:2){
    if(combo==1){virus1='8';location='Sassari'}
    if(combo==2){virus1='6';location='Teramo'}
    
      if(file.exists(paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQ',
                           translator_loc[[location]],virus1,'dpi',lev2,'.csv',sep=''))){
        print(c('yay',location,virus1))
        dpi_levels[2]=lev2
        dat1=read.csv(file=paste('~/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQ',
                                 translator_loc[[location]],virus1,'dpi',lev2,'.csv',sep=''))
        virus=grep(virus1,virus_names,value=T)
        dat1=dat1[,c(1,1,6,2)]
        colnames(dat1)
        if(!plot_nonsig){
          dat1=dat1[dat1[,3]<0.05,]
        }
        # for(i in 1:2){
        #   for(j in (i+1):3){
        dat1=dat1[dat1[,2]%in%common_ISGS,]
        plot_dat=rbind(plot_dat,as.matrix(cbind(dat1[,4],dat1[,3]<0.05,
                                                rep(paste(location,'.',virus,sep=''),nrow(dat1)),
                                                rep(paste(location,'.',virus,'\n(N=',
                                                          length(which(dat1[,3]<0.05)),')',sep=''),nrow(dat1))
        )))
      }
    
  }
  #   
  # plot_dat=rbind(plot_dat,
  #                 as.matrix(cbind(rep(0,nrow(ISG_list_ov)),rep(T,nrow(ISG_list_ov)),
  #                                      rep('crap',nrow(ISG_list_ov)), rep('crap',nrow(ISG_list_ov)))))
  plot_dat_save=plot_dat
  
  plot_dat=plot_dat_save
  plot_dat=as.data.frame(plot_dat)
  colnames(plot_dat)=c('logFC','significant','condition','condition_lab')
  plot_dat$colour=plot_dat$condition
  plot_dat$colour[!as.logical(plot_dat$significant)]='nonsig'
  plot_dat$colour=factor(plot_dat$colour,levels=c('nonsig',unique(plot_dat$condition)))
  colpal=c(scales::alpha('#AAAAAA',.2),scales::alpha(RColorBrewer::brewer.pal(length(unique(plot_dat$condition)),'Dark2'),.5))
  plot_dat$condition=factor(plot_dat$condition,levels=unique(plot_dat$condition))
  
  plot_dat$condition_lab=factor(plot_dat$condition_lab,levels=unique(plot_dat$condition_lab))
  
  plot_dat$logFC=as.numeric(plot_dat$logFC)
  plot_dat=plot_dat[order(plot_dat$colour),]
  plots[[plot_i]]= ggplot(plot_dat,aes(x=condition_lab,y=logFC,color=colour))+
    geom_jitter()+scale_colour_manual(values=colpal[sort(unique(as.numeric(plot_dat$colour)))])+theme_bw()+ylim(-3,8)+
    xlab('virus:location')+
    ggtitle(paste('ISG expression dpi',lev2))
  # plots[[plot_i]]= ggplot(plot_dat,aes(x=condition_lab,y=logFC,fill=colour))+
  #   geom_violin(scale='count')+scale_fill_manual(values=colpal)+theme_bw()+ylim(-5,10)
  
  
}
library(gridExtra)
grid.arrange(plots[[1]],plots[[2]],ncol=1)

dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi1.csv')
hist(p.adjust(dat$P.Value,method='hochberg'))
hist(dat$adj.P.Val)



dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi1.csv')
hist((dat$P.Value))

dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi1.csv')
dat2=dat[dat[,1]%in%common_ISGS,]
var(dat2$logFC)
dat3=dat[!dat[,1]%in%common_ISGS,]
var(dat3$logFC)
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi1.csv')
dat3=dat[dat[,1]%in%common_ISGS,]
var(dat3$logFC)
dat3=dat[!dat[,1]%in%common_ISGS,]
var(dat3$logFC)
################################################################################################################################################
################################################################################################################################################
############################################################################################################
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi1.csv')
dat2=dat[dat[,1]%in%common_ISGS,]

dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi1.csv')
dat3=dat[dat[,1]%in%common_ISGS,]

par(mfrow=c(1,2))

dat2=dat2[dat2$X%in%dat3$X,]
dat3=dat3[match(dat2$X,dat3$X),]

plot(dat2$logFC,dat3$logFC,pch=19,col=scales::alpha('#AAAAAA',.4),xlim=c(-2,8),ylim=c(-2,3),
     ylab='BTV 8 dpi1 logFC',xlab='BTV-1 2006 dpi1 logFC',main='ISGs')

####################################
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi3.csv')
dat2=dat[dat[,1]%in%common_ISGS,]

dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi1.csv')
dat3=dat[dat[,1]%in%common_ISGS,]

# par(mfrow=c(1,2))

dat2=dat2[dat2$X%in%dat3$X,]
dat3=dat3[match(dat2$X,dat3$X),]

plot(dat2$logFC,dat3$logFC,pch=19,col=scales::alpha(c('#AAAAAA','#D95F02')[1+as.numeric(dat2$adj.P.Val<0.05)],.4),xlim=c(-2,8),ylim=c(-2,3),
     ylab='BTV 8 dpi1 logFC',xlab='BTV-1 2006 dpi3 logFC',main='ISGs')
################################################################################################################################################
# ################################################################################################################################################
# ############################################################################################################
# dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi1.csv')
# dat2=dat[!dat[,1]%in%common_ISGS,]
# 
# dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi1.csv')
# dat3=dat[!dat[,1]%in%common_ISGS,]
# 
# par(mfrow=c(1,2))
# 
# dat2=dat2[dat2$X%in%dat3$X,]
# dat3=dat3[match(dat2$X,dat3$X),]
# 
# plot(dat2$logFC,dat3$logFC,pch=19,col=scales::alpha('#AAAAAA',.4),xlim=c(-2,8),ylim=c(-3,4),
#      ylab='BTV 8 dpi1 logFC',xlab='BTV-1 2006 dpi1 logFC',main='non-ISGs')
# 
# ####################################
# dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi3.csv')
# dat2=dat[!dat[,1]%in%common_ISGS,]
# 
# dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSMI8dpi1.csv')
# dat3=dat[!dat[,1]%in%common_ISGS,]
# 
# # par(mfrow=c(1,2))
# 
# dat2=dat2[dat2$X%in%dat3$X,]
# dat3=dat3[match(dat2$X,dat3$X),]
# 
# plot(dat2$logFC,dat3$logFC,pch=19,col=scales::alpha('#AAAAAA',.4),xlim=c(-2,8),ylim=c(-3,4),
#      ylab='BTV 8 dpi1 logFC',xlab='BTV-1 2006 dpi3 logFC',main='non-ISGs')
