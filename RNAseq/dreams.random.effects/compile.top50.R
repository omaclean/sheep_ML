


dir="~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/"
gene_names=read.csv('~/Documents/orthologue_sets/Sheep_mega_data.list.csv')
setwd(dir)

files=list.files('.',pattern='.csv')
i=files[1]
for(i in files){
  dat=read.csv(i)
  dat50=dat[dat$adj.P.Val<0.05,]
  dat50=dat[order(dat$logFC ,decreasing=T)[1:50],]
  dat50$gene_names=gene_names$gene_names_new[match(dat50$X,gene_names$ens_IDs)]
  dat50=dat50[,c(1,ncol(dat50),2:(ncol(dat50)-1))]
  write.csv(dat50,file=paste('top50/',i,sep=''))
}


dir="~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/top50"

setwd(dir)

files=list.files('.',pattern='.csv')
i=files[1]

hits=matrix(ncol=2,nrow=0)
for(i in files){
  dat=read.csv(i,row.names = 1)
  hits=rbind(hits,dat[,1:2])
}

hit_counts=cbind(names(table(hits[,1])),gene_names$gene_names_new[match(names(table(hits[,1])),gene_names$ens_IDs)],
                 table(hits[,1]))
hit_counts=hit_counts[order(hit_counts[,3],decreasing = T),]
head(hit_counts)
write.csv(hit_counts,'tophits.csv')


#############





cats=c('SMI8dpi1','SMI13dpi1','SMI6dpi1',
       'SLI13dpi1','SLI6dpi1',
       'SMI8dpi3','SMI13dpi3','SMI6dpi3',
       'SLI13dpi3','SLI6dpi3',
       'SMI8dpi7','SMI13dpi7','SMI6dpi7',
       'SLI13dpi7','SLI6dpi7')


clin=read.csv('/home/oscar/Documents/sheep_megadata/14.9.21/RF_input_data/clinical_score.data.csv')
four=read.csv('/home/oscar/Documents/sheep_megadata/14.9.21/RF_input_data/final_four.data.csv')
six=read.csv('/home/oscar/Documents/sheep_megadata/14.9.21/RF_input_data/final_seven.data.csv')

dir="~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/"
setwd(dir)
files=list.files('.',pattern='.csv')
files
comb=intersect(colnames(clin),intersect(colnames(four),colnames(six)))
top.btms=grep('BTM',comb,value=T)
mat=read.csv('/home/oscar/Documents/sheep_megadata/14.9.21/dat_plots_filt/BTMs.gene.sets.csv',row.names = 1)
mat[,1]=paste('BTM.',mat[,1],sep='')

mat[,1]=gsub('\\)|\\ |\\(','.',mat[,1])
head(mat)
for(i in 1:length(top.btms)){
  genes=mat[mat[,1]==top.btms[i],2]
  grep('RIG.1.like.receptor',mat[,1],value=T)[1]
  gene_num=length(genes)
  vals=data.frame(cats=rep(cats,gene_num),FDR=rep(NA,length(cats)*gene_num),
                  logFC=rep(NA,length(cats)*gene_num),
                  genes=rep(genes,each=length(cats)))
  plots=list()
  for(j in 1:length(cats)){
    dat_j=read.csv(paste(cats[j],'.csv',sep=''))
    dat_j=dat_j[dat_j$X%in%genes,]
    if(nrow(dat_j)==0){print(paste(cats[j],top.btms[i]));next}
    for(k in 1:gene_num){
      if(genes[k]%in%dat_j$X){
        vals[vals$cats==cats[j]&vals$genes==genes[k],2:3]=
          c(dat_j$adj.P.Val[dat_j$X==genes[k]],dat_j$logFC[dat_j$X==genes[k]])  
      }
    }
  }
  head(vals)
  for(j in 1:length(unique(vals$genes))){
    if(nrow(vals)==0){print(top.btms[i]);next}
    gene_j=unique(vals$genes)[j]
    vals_gene=vals[vals$genes==gene_j,]
    vals_gene=vals_gene[!is.na(vals_gene$FDR),]
    vals_gene$cats2=factor(vals_gene$cats,levels=cats)
    vals_gene$signif=vals_gene$FDR<0.05
    plots[[j]]=ggplot(vals_gene,aes(y=logFC,x=cats2,fill=signif))+
      geom_bar(stat='identity')+theme_bw()+
      ggtitle(paste(gene_j,gene_names$gene_names_new[gene_names$ens_IDs==gene_j]))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ylim(min(vals$logFC[!is.na(vals$logFC)]),max(vals$logFC[!is.na(vals$logFC)]))
  }
  png(paste('~/Pictures/plots/Sheep_megadata/14.9.21/RF_plots/BTM.top.hits.barplots/',
            substr(top.btms[i],1,30),'.png',sep=''),width=1000,height=100*length(plots))
  do.call("grid.arrange", c(plots,ncol=3))
  dev.off()
  
}
