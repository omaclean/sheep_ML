counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

cond_levels=c('SMC','SLC','SMI6','SLI6','SMI13','SLI13','SMI8')
dat=matrix(ncol=5,nrow=1)
for(dpi in c('0','1','3','7')){
  for(cond in cond_levels){
    counts_ij=counts[,grepl(cond,colnames(counts))&
                       sapply(colnames(counts),function(x) strsplit(x,'_')[[1]][3])==dpi]
    print(c(dpi,cond,ncol(counts_ij),length(which(rowSums(counts_ij>0)==ncol(counts_ij))),
            length(which(rowSums(counts_ij>0)>6))))
    dat=rbind(dat,c(dpi,cond,ncol(counts_ij),length(which(rowSums(counts_ij>0)==ncol(counts_ij))),
                          length(which(rowSums(counts_ij>0)>6))))
  }
}
write.csv(dat,file='~/read.info.csv')

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
# conditions=c(conditions,'BTV8_dpi1')
# cond_i[dpis==1&conds=='SMI8']=conditions[5]
# conditions=c(conditions,'BTV8_dpi3')
# cond_i[dpis==3&conds=='SMI8']=conditions[6]
# conditions=c(conditions,'BTV8_dpi7')
# cond_i[dpis==7&conds=='SMI8']=conditions[7]

  
cond_i=paste(sex,sapply(conds,function(x) substr(x,1,2)),cond_i)

pal=alpha(brewer.pal(length(unique(conditions)),'Dark2'),0.7)

# ##
# c('dpi 0', ' dpi 1',  'dpi3-2006',  'dpi3-2013',  'dpi7-2006','dpi7-2013')
# c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3')


library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)

plot_ly(x=PCA$x[,1],y=PCA$x[,2],z=PCA$x[,3],color=factor(cond_i),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=10))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))



plot_ly(x=PCA$x[,4],y=PCA$x[,5],z=PCA$x[,6],color=factor(cond_i),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA6 (var prop=",round(100*PCA$sdev[6]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))

