counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)

counts2=apply(counts[,2:ncol(counts)],2,function(x) x/sum(x))
library(plotly);library(RColorBrewer);library(scales)

cond_levels=c('SMC','SLC','SMI6','SLI6','SMI13','SLI13','SMI8')
conds=sapply(colnames(counts2),function(x) strsplit(x,'\\.')[[1]][1])
dpis=as.numeric(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][3]))
################
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
conditions=

pal=alpha(brewer.pal(length(unique(conditions)),'Dark2'),0.7)

# ##
# c('dpi 0', ' dpi 1',  'dpi3-2006',  'dpi3-2013',  'dpi7-2006','dpi7-2013')
# c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3')


library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)

plot_ly(x=PCA$x[,1],y=PCA$x[,2],z=PCA$x[,3],color=factor(cond_i,levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=16))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[,4],y=PCA$x[,5],color=factor(cond_i,levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[,2],y=PCA$x[,3],color=factor(cond_i,levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[,2],y=PCA$x[,3],color=factor(cond_i,levels=conditions),colors=brewer.pal(7,'Dark2')[c(1:3,5,4,6,7)],alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))



##

    round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1)
 
library(Rtsne)
tsne=Rtsne(t(counts3),dims=2, perplexity=30,verbose=T,max_iter=500)

plot_ly(x=tsne$Y[,1],y=tsne$Y[,2],color=factor(cond_i,levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=24))%>%
  layout(scene=list(xaxis=list(title="tSNE1"),yaxis=list(title="tSNE2")),
         legend=list(font=list(size=20)))
library(umap)
UMAP=umap(t(counts3))

plot_ly(x=UMAP$layout[,1],y=UMAP$layout[,2],color=factor(cond_i,levels=conditions),colors=brewer.pal(length(conditions),'Dark2'),alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=24))%>%
  layout(scene=list(xaxis=list(title="UMAP1"),yaxis=list(title="UMAP2")),
         legend=list(font=list(size=20)))
##sex

sex2=sex
sex2[grepl('SMI',animal_codes)&sex=='male']='male-Sassari'
plot_ly(x=UMAP$layout[,1],y=UMAP$layout[,2],color=factor(sex2,levels=unique(sex2)),colors=brewer.pal(length(conditions),'Accent')[1:3],alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=24))%>%
  layout(scene=list(xaxis=list(title="UMAP1"),yaxis=list(title="UMAP2")),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[,2],y=PCA$x[,3],color=factor(sex2,levels=unique(sex2)),colors=brewer.pal(length(conditions),'Accent')[1:3],alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex))%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))



#########
###################
###################
###################
###################
###################
###################
###################
###################
###################
###################
#############################


pal=alpha(brewer.pal(length(unique(conditions)),'Dark2'),0.7)
library(scales)
## more cols
conditions=c('uninfected', 'dpi 1 2006&13',  'dpi3-2006',  'dpi3-2013',  'dpi7-2006','dpi7-2013','dpi1,3,7-BTV8')
pal=alpha(c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3','#666666'),0.6)


clinical_score=rep(0,ncol(counts2))
animal_codes=sapply(colnames(counts2),function(x) gsub('\\.','-',
                                                       paste(strsplit(x,'_')[[1]][1])))

clinical_score_file=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv',stringsAsFactors = F)
cond_i=rep(NA,ncol(counts3))
for(i in 1:ncol(counts2)){
  clinical_score[i]=clinical_score_file$clinical.score[clinical_score_file$ID==animal_codes[i]&
                                                         clinical_score_file$dpi==dpis[i]]
  if(conds[i]=='SMC'|conds[i]=='SLC'|dpis[i]==0){
    cond_i[i]=conditions[1]
  }else{
    if(dpis[i]==1){
      cond_i[i]=conditions[2]
    }
    if(dpis[i]==3&(conds[i]=='SMI6'|conds[i]=='SLI6')){
      cond_i[i]=conditions[3]
    }    
    if(dpis[i]==3&(conds[i]=='SMI13'|conds[i]=='SLI13')){
      cond_i[i]=conditions[4]
    }
    if(dpis[i]==7&(conds[i]=='SMI6'|conds[i]=='SLI6')){
      cond_i[i]=conditions[5]
    }    
    if(dpis[i]==7&(conds[i]=='SMI13'|conds[i]=='SLI13')){
      cond_i[i]=conditions[6]
    }
    if(dpis[i]>0&conds[i]=='SMI8'){
      cond_i[i]=conditions[7]
    }
  }
}


cond_i


library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)

plot_ly(x=PCA$x[!is.na(cond_i),1],y=PCA$x[!is.na(cond_i),2],
        z=PCA$x[!is.na(cond_i),3],color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[!is.na(cond_i),4],y=PCA$x[!is.na(cond_i),5],
        color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))
###############


## more cols
conditions=c('uninfected', ' dpi 1',  'dpi3-2006-Teramo',
             'dpi3-2006-Sassari',  'dpi7-2006-Teramo','dpi7-2006-Sassari')
pal=alpha(c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3'),0.6)

cond_i=rep(NA,ncol(counts3))
for(i in 1:ncol(counts2)){
  clinical_score[i]=clinical_score_file$clinical.score[clinical_score_file$ID==animal_codes[i]&
                                                         clinical_score_file$dpi==dpis[i]]
  if(conds[i]=='SMC'|conds[i]=='SLC'|dpis[i]==0){
    cond_i[i]=conditions[1]
  }else{
    if(dpis[i]==1){
      cond_i[i]=conditions[2]
    }
    if(dpis[i]==3&(conds[i]=='SLI6')){
      cond_i[i]=conditions[3]
    }    
    if(dpis[i]==3&(conds[i]=='SMI6')){
      cond_i[i]=conditions[4]
    }
    if(dpis[i]==7&(conds[i]=='SLI6')){
      cond_i[i]=conditions[5]
    }    
    if(dpis[i]==7&(conds[i]=='SMI6')){
      cond_i[i]=conditions[6]
    }
  }
}


cond_i


library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)

plot_ly(x=PCA$x[!is.na(cond_i),1],y=PCA$x[!is.na(cond_i),2],
        z=PCA$x[!is.na(cond_i),3],color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=8))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[!is.na(cond_i),4],y=PCA$x[!is.na(cond_i),5],
        color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))

##############


## more cols
conditions=c('uninfected', ' dpi 1',  'dpi3-2013-Teramo',
             'dpi3-2013-Sassari',  'dpi7-2013-Teramo','dpi7-2013-Sassari')
pal=alpha(c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3'),0.6)

cond_i=rep(NA,ncol(counts3))
for(i in 1:ncol(counts2)){
  clinical_score[i]=clinical_score_file$clinical.score[clinical_score_file$ID==animal_codes[i]&
                                                         clinical_score_file$dpi==dpis[i]]
  if(conds[i]=='SMC'|conds[i]=='SLC'|dpis[i]==0){
    cond_i[i]=conditions[1]
  }else{
    if(dpis[i]==1){
      cond_i[i]=conditions[2]
    }
    if(dpis[i]==3&(conds[i]=='SLI13')){
      cond_i[i]=conditions[3]
    }    
    if(dpis[i]==3&(conds[i]=='SMI13')){
      cond_i[i]=conditions[4]
    }
    if(dpis[i]==7&(conds[i]=='SLI13')){
      cond_i[i]=conditions[5]
    }    
    if(dpis[i]==7&(conds[i]=='SMI13')){
      cond_i[i]=conditions[6]
    }
  }
}

library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)

plot_ly(x=PCA$x[!is.na(cond_i),1],y=PCA$x[!is.na(cond_i),2],
        z=PCA$x[!is.na(cond_i),3],color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=8))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[!is.na(cond_i),4],y=PCA$x[!is.na(cond_i),5],
        color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))
########


## more cols
conditions=c('uninfected', ' dpi 1',  'dpi3-Teramo',
             'dpi3-Sassari',  'dpi7-Teramo','dpi7-Sassari')
pal=alpha(c("#1B9E77","#D95F02",'#44487e'  ,'#9c8fbd' ,'#f94e87','#c322a3'),0.6)

cond_i=rep(NA,ncol(counts3))
for(i in 1:ncol(counts2)){
  clinical_score[i]=clinical_score_file$clinical.score[clinical_score_file$ID==animal_codes[i]&
                                                         clinical_score_file$dpi==dpis[i]]
  if(conds[i]=='SMC'|conds[i]=='SLC'|dpis[i]==0){
    cond_i[i]=conditions[1]
  }else{
    if(dpis[i]==1){
      cond_i[i]=conditions[2]
    }
    if(dpis[i]==3&(grepl('SLI',conds[i]))){
      cond_i[i]=conditions[3]
    }    
    if(dpis[i]==3&(grepl('SMI',conds[i]))){
      cond_i[i]=conditions[4]
    }
    if(dpis[i]==7&(grepl('SLI',conds[i]))){
      cond_i[i]=conditions[5]
    }    
    if(dpis[i]==7&(grepl('SMI',conds[i]))){
      cond_i[i]=conditions[6]
    }
  }
}

library(plotly)
counts3=counts2[rowSums(counts2>0)>7,]
PCA=prcomp(t(counts3),scale=T)

plot_ly(x=PCA$x[!is.na(cond_i),1],y=PCA$x[!is.na(cond_i),2],
        z=PCA$x[!is.na(cond_i),3],color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))

plot_ly(x=PCA$x[!is.na(cond_i),4],y=PCA$x[!is.na(cond_i),5],
        color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title=paste("PCA4 (var prop=",round(100*PCA$sdev[4]^2/sum(PCA$sdev^2),1),'%)')),
         yaxis=list(title=paste("PCA5 (var prop=",round(100*PCA$sdev[5]^2/sum(PCA$sdev^2),1),'%)')),
         legend=list(font=list(size=20)))

# 
# 
# 
# 
# plot(x=UMAP$layout[,1],y=UMAP$layout[,2],col=brewer.pal(length(loc_x_gen_cols),'Set2')[as.numeric(loc_x_gen_cols_i)],
#      pch=as.numeric(as.factor(cond_i)))
# legend(pch=unique(as.numeric(as.factor(cond_i))),x='bottomright',legend=levels(as.factor(cond_i)),bty='n')
# legend(col=brewer.pal(length(loc_x_gen_cols),'Dark2'),x='topleft',
#        legend=loc_x_gen_cols,bty='n',lwd=4)


