counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

counts2=apply(counts[,2:ncol(counts)],2,function(x) x/sum(x))
library(plotly);library(RColorBrewer);library(scales)

cond_levels=c('SMC','SLC','SMI6','SLI6','SMI13','SLI13','SMI8')
conds=sapply(colnames(counts2),function(x) strsplit(x,'\\.')[[1]][1])
dpis=as.numeric(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][3]))
cond_i=rep(NA,ncol(counts2))
conditions=c('uninfected','dpi1','dpi3','dpi7')
clinical_score=rep(0,ncol(counts2))
animal_codes=sapply(colnames(counts2),function(x) gsub('\\.','-',
                                                       paste(strsplit(x,'_')[[1]][1])))
gender=read.csv('/home/oscar/Documents/sheep_megadata/temperature_humidity_gender.csv')
sex=rep(NA,length(animal_codes))
for(i in 1:length(animal_codes)){
  sex[i]=c('female','male')[gender$gender.1[which(gender$animal.number==animal_codes[i])]]
}


clinical_score_file=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv',stringsAsFactors = F)
for(i in 1:ncol(counts2)){
  clinical_score[i]=clinical_score_file$clinical.score[clinical_score_file$ID==animal_codes[i]&
                                                         clinical_score_file$dpi==dpis[i]]
  if(conds[i]=='SMC'|conds[i]=='SLC'|dpis[i]==0){
    cond_i[i]=conditions[1]
  }else{
    if(dpis[i]==1){
      cond_i[i]=conditions[2]
    }
    if(dpis[i]==3){
      cond_i[i]=conditions[3]
    }
    if(dpis[i]==7){
      cond_i[i]=conditions[4]
    }
  }
}


pal=alpha(brewer.pal(length(unique(conditions)),'Dark2'),0.7)

# ##
# c('dpi 0', ' dpi 1',  'dpi3-2006',  'dpi3-2013',  'dpi7-2006','dpi7-2013')
# c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3')


library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)
conditions=rep(0,nrow(PCA$x))
conditions[grepl('SMC',animal_codes)]='SMC dpi0,1,3'
conditions[dpis=='7'&grepl('SMC',animal_codes)]='SMC dpi7'
conditions[dpis=='7'&grepl('SLI6',animal_codes)]='SLI6dpi7'
conditions[dpis=='7'&grepl('SMI13',animal_codes)]='SMI13dpi7'
plot_ly(x=PCA$x[,1],y=PCA$x[,2],z=PCA$x[,3],color=factor(conditions),colors=brewer.pal(length(unique(conditions)),'Dark2'),alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))
########################


library(psych);library(ggplot2);library(caret);library(randomForest);library(RColorBrewer)
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
redundant_BTMS=c("Rho.GTPase.cycle..M4.14.","mitotic.cell.division..M6.",
                 "enriched.in.NK.cells..KIR.cluster...M61.1.","enriched.in.B.cells..VI...M69.")

dat=dat[,!colnames(dat)%in%redundant_BTMS]


for(i in 2:ncol(counts)){
  counts[,i]=counts[,i]/sum(counts[,i])
}

for(i in 1:nrow(counts)){
  counts[i,2:ncol(counts)]=counts[i,2:ncol(counts)]/sum(counts[i,2:ncol(counts)])
}

any(apply(counts[,2:ncol(counts)],2,function(x) sum(x))==0)

## melt BTMs
mat=matrix(ncol=2,nrow=500000)
total=1
for(i in 1:ncol(dat)){
  genes_i=unique(dat[2:nrow(dat),i])
  genes_i=as.character(genes_i[genes_i!=' '])
  len=length(genes_i)
  mat[total:(total+len-1),]=cbind(as.character(rep(names(dat)[i],len)),genes_i)
  total=total+len-1
}
mat=mat[!is.na(mat[,1]),]
mat=mat[mat[,2]%in%counts[,1],]
## get geometric means out for counts

#take only BTMs with at least 2 genes with counts
BTMs_work=unique(mat[duplicated(mat[,1]),1])
mat2=matrix(ncol=ncol(counts[2:ncol(counts)]),nrow=length(BTMs_work))
colnames(mat2)=colnames(counts[,2:ncol(counts)])
rownames(mat2)=BTMs_work
for(i in 1:length(BTMs_work)){
  for(j in 2:ncol(counts)){
    mat2[i,j-1]=geometric.mean(1+counts[counts$Geneid%in%
                                          mat[mat[,1]==BTMs_work[i],2]
                                        ,j])-1
  }
}
dim(dat)
which(is.na(mat2),arr.ind=T)
any(apply(mat2,1,function(x) sum(x))==0)

mat2[is.na(mat2[,1]),]
###########################
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
par(mfrow=c(3,2),mar=c(4,4,1,1))
mat3=mat2[rowSums(mat2>0)>1,]
pca_btm=prcomp(t(mat3),retx=T,scale.=T,center=T)

  pal=alpha(brewer.pal(length(unique(conditions)),'Dark2'),0.7)
  
  # ##
  # c('dpi 0', ' dpi 1',  'dpi3-2006',  'dpi3-2013',  'dpi7-2006','dpi7-2013')
  # c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3')
  
  
  library(plotly)
  counts3=counts2[rowSums(counts2>0)>7,]
  PCA=prcomp(t(counts3),scale=T)
  conditions=rep(0,nrow(PCA$x))
  conditions[grepl('SMC',animal_codes)]='SMC dpi0,1,3'
  conditions[dpis=='7'&grepl('SMC',animal_codes)]='SMC dpi7'
  conditions[dpis=='7'&grepl('SLI6',animal_codes)]='SLI6dpi7'
  conditions[dpis=='7'&grepl('SMI13',animal_codes)]='SMI13dpi7'
  plot_ly(
    x=pca_btm$x[,1],y=pca_btm$x[,2],
    z=pca_btm$x[,3],color=factor(conditions),colors=brewer.pal(length(unique(conditions)),'Dark2'),alpha=0.85,
    text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
    add_markers(marker=list(size=12))%>%
    layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*pca_btm$sdev[1]^2/sum(pca_btm$sdev^2),1),'%)')),
                      yaxis=list(title=paste("PCA2 (var prop=",round(100*pca_btm$sdev[2]^2/sum(pca_btm$sdev^2),1),'%)')),
                      zaxis=list(title=paste("PCA3 (var prop=",round(100*pca_btm$sdev[3]^2/sum(pca_btm$sdev^2),1),'%)'))),
           legend=list(font=list(size=20)))
  
  
  
  conditions=rep(0,nrow(PCA$x))
  conditions[dpis=='7'&grepl('SMC',animal_codes)]='SMC dpi7'

  rf=randomForest(y= as.factor(conditions),x= t(mat2),
               importance=T,do.trace = F)
names(rf$importance[order(rf$importance[,1],decreasing = T),1])

rf$importance
BTMs_to_plot=names(rf$importance[order(rf$importance[,4],decreasing = T),4])[1:11]

head(rf$importance[order(rf$importance[,4],decreasing = T),4])

mat_plot=mat2[match(BTMs_to_plot,rownames(mat2)),]

dpis=sapply(colnames(mat_plot),function(x) strsplit(x,'_')[[1]][3])
treatments=sapply(colnames(mat_plot),function(x) strsplit(x,'\\.')[[1]][1])
colnames(mat_plot)=paste(treatments,dpis,sep='_')

cond_list_order=c('SMC_0','SMC_1','SMC_3','SMC_7',
                  'SLC_0','SLC_1','SLC_3','SLC_7',
                  'SMI8_0','SMI13_0','SMI6_0',
                  'SLI13_0','SLI6_0',
                  'SMI8_1','SMI13_1','SMI6_1',
                  'SLI13_1','SLI6_1',
                  'SMI8_3','SMI13_3','SMI6_3',
                  'SLI13_3','SLI6_3',
                  'SMI8_7','SMI13_7','SMI6_7',
                  'SLI13_7','SLI6_7')
types_convert_plot=list()
types_convert_plot[['SMC']]=types_convert_plot[['SLC']]='control'
types_convert_plot[['SLI8']]=types_convert_plot[['SMI8']]=paste('BTV8',sep='')
types_convert_plot[['SMI13']]=types_convert_plot[['SLI13']]=paste('2013',sep='')
types_convert_plot[['SMI6']]=types_convert_plot[['SLI6']]=paste('2006',sep='')


plots=list()
for(i in 1:nrow(mat_plot)){
  BTM_to_plot=i
  to_plot=data.frame(vals=mat_plot[BTM_to_plot,],conds=colnames(mat_plot),treatment=treatments)
  to_plot$conds=factor(to_plot$conds,levels = cond_list_order)
  to_plot$treatment_col=factor(types_convert_plot[to_plot$treatment],levels=c('control','BTV8','2013','2006'))
  plots[[i]]=ggplot(to_plot,aes(x=conds,y=vals,fill=treatment_col))+
    geom_dotplot(binaxis='y',stackdir = 'center',alpha=0.3,dotsize = 0.9)+theme_bw()+
    ggtitle(rownames(mat_plot)[i])+
    geom_vline(xintercept=(13.5))+
    geom_vline(xintercept=(18.5))+
    geom_vline(xintercept=(23.5))
}
library(gridExtra)
png(filename=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn.new/test.SMC7.var.filt.png',sep=''),width=1100,height=200*length(plots))

grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]]
             ,plots[[5]], plots[[6]],plots[[7]],plots[[8]],plots[[9]]
             ,plots[[10]],plots[[11]]#,plots[[12]]
             ,ncol=1)
dev.off()    
  