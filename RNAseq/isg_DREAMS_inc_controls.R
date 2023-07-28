####~/scripts/R/sheep_mega_data/rna_seq/violin_and_dotplot_dreams.R
#################################

library(scales);library(RColorBrewer);library(ggplot2);library(gridExtra)
par(mfrow=c(1,1))


ISG_list_ov=read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv",stringsAsFactors = F)
dpi_levels=c('0',NA)
#out_dir="/home/oscar/Pictures/plots/RNA-seq/ov_vs_bov/all_genes.11.11.19"#no/ doesn't work
common_ISGS=unique(c(ISG_list_ov$Gene.ID))


############################## PLUS controls
plot_nonsig=T

plots=list()
plot_i=0
translator_loc=list()
locations=c('Teramo','Sassari')
translator_loc[['Sassari']]='SM'
translator_loc[['Teramo']]='SL'
viruses=c('C','I8','I13','I6')
virus_names=c('control','BTV8','2013','2006')
for(lev2 in c('1','3','7')){
  plot_i=plot_i+1
  plot_dat=matrix(nrow=0,ncol=4)
  for(virus_1 in 1:length(viruses)){
    virus1=viruses[virus_1]
    for(location in locations){
      if(virus1!='C'){
        file_in=paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/',
                      translator_loc[[location]],virus1,'dpi',lev2,'.csv',sep='')
      }else{file_in=paste('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/no_controls_infect_dpi0_other/',
                          translator_loc[[location]],virus1,'dpi',lev2,'.csv',sep='')}
      if(file.exists(file_in)){
        print(c('yay',location,virus1))
        dpi_levels[2]=lev2
        dat1=read.csv(file=file_in)
        virus=virus_names[virus_1]
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
grid.arrange(plots[[1]],plots[[2]],plots[[3]],ncol=1)