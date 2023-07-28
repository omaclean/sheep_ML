counts=read.table('/home/oscar/Documents/sheep_megadata/tru_seq_29.6.21/All_Count.txt',header=T)
#counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
cond_levels=c('SLI6','SMI8')
dim(counts[,grepl(paste(cond_levels,collapse='|'),colnames(counts))])
counts=counts[,grepl(paste(cond_levels,collapse='|'),colnames(counts))]


dim(counts)

counts2=apply(counts[,1:ncol(counts)],2,function(x) x/sum(x))
dpis=as.numeric(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][3]))

counts2=counts2[,dpis%in%c(0,1,3)]
dim(counts2)
library(plotly);library(RColorBrewer);library(scales)


conds=sapply(colnames(counts2),function(x) strsplit(x,'\\.')[[1]][1])
dpis=as.numeric(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][3]))

counts2=counts2[,dpis%in%c(0,1,3)]
animal_codes=sapply(colnames(counts2),function(x) gsub('\\.','-',
                                                       paste(strsplit(x,'_')[[1]][1])))
gender=read.csv('/home/oscar/Documents/sheep_megadata/temperature_humidity_gender.csv')
sex=rep(NA,length(animal_codes))
for(i in 1:length(animal_codes)){
  sex[i]=c('female','male')[gender$gender.1[which(gender$animal.number==animal_codes[i])]]
}


## more cols
conditions=c('BTV8-dpi0', 'BTV8-dpi1',  'BTV8-dpi3',  'T2006-dpi0',  'T2006-dpi1','T2006-dpi3')

pal=alpha(c('#666666','#9c8fbd','#44487e',"#D95F02","#1B9E77",'#f94e87'),0.6)
plot(1:7,col=pal,pch=19)

clinical_score=rep(0,ncol(counts2))
animal_codes=sapply(colnames(counts2),function(x) gsub('\\.','-',
                                                       paste(strsplit(x,'_')[[1]][1])))

clinical_score_file=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv',stringsAsFactors = F)
cond_i=rep(NA,ncol(counts2))
for(i in 1:ncol(counts2)){
  clinical_score[i]=clinical_score_file$clinical.score[clinical_score_file$ID==animal_codes[i]&
                                                         clinical_score_file$dpi==dpis[i]]
  if(dpis[i]==0&conds[i]=='SMI8'){
    cond_i[i]=conditions[1]
  }
  if(dpis[i]==1&conds[i]=='SMI8'){
    cond_i[i]=conditions[2]
  }    
  if(dpis[i]==3&conds[i]=='SMI8'){
    cond_i[i]=conditions[3]
  }
  if(dpis[i]==0&conds[i]=='SLI6'){
    cond_i[i]=conditions[4]
  }    
  if(dpis[i]==1&conds[i]=='SLI6'){
    cond_i[i]=conditions[5]
  }
  if(dpis[i]==3&conds[i]=='SLI6'){
    cond_i[i]=conditions[6]
  }
}


cond_i

counts2=counts2[rowSums(counts2>0)>6,]

PCA=prcomp(t(counts2),scale=T)


plot_ly(x=PCA$x[!is.na(cond_i),1],y=PCA$x[!is.na(cond_i),2],
        z=PCA$x[!is.na(cond_i),3],color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=12))%>%
  layout(scene=list(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)')),
                    yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)')),
                    zaxis=list(title=paste("PCA3 (var prop=",round(100*PCA$sdev[3]^2/sum(PCA$sdev^2),1),'%)'))),
         legend=list(font=list(size=20)))




plot_ly(x=PCA$x[!is.na(cond_i),1],y=PCA$x[!is.na(cond_i),2],
        color=factor(cond_i[!is.na(cond_i)],levels=conditions),colors=pal,alpha=0.85,
        text=~paste(sapply(colnames(counts2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis,sex)[!is.na(cond_i)])%>%
  add_markers(marker=list(size=15))%>%
  layout(xaxis=list(title=paste("PCA1 (var prop=",round(100*PCA$sdev[1]^2/sum(PCA$sdev^2),1),'%)'),
                    titlefont = list(size = 15)),
         yaxis=list(title=paste("PCA2 (var prop=",round(100*PCA$sdev[2]^2/sum(PCA$sdev^2),1),'%)'),
                    titlefont = list(size = 15)),
         legend=list(font=list(size=18)))

dim(counts2)
