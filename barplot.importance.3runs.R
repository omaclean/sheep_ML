###############################################

outdir="/home/oscar/scripts/github/sheep_ML/outdir/barplots"

N_clinical=100
clin=read.csv(paste('/home/oscar/Pictures/plots/Sheep_megadata/15.5.22/clinical_score.data_',
                          N_clinical,'.params.csv',sep=''))
four=read.csv('~/Pictures/plots/Sheep_megadata/15.5.22/final_four.data.17params.csv')
six=read.csv('~/Pictures/plots/Sheep_megadata/15.5.22/final_six_states.data50params.csv')


#################code
outdir=gsub("/$","",outdir)
shared_params=intersect(intersect(names(clin)[2:ncol(clin)],names(six)[2:ncol(six)]),
                        names(four)[2:ncol(four)])

shared_once_params=unique(c(intersect(names(clin)[2:ncol(clin)],names(six)[2:ncol(six)]),
                            intersect(names(six)[2:ncol(six)],names(four)[2:ncol(four)]),
                            intersect(names(clin)[2:ncol(clin)],names(four)[2:ncol(four)])
))

param_sets=list(clinical=names(clin)[2:ncol(clin)],six_states=names(six)[2:ncol(six)],
                four_states=names(four)[2:ncol(four)])

venn::venn(snames=c(paste('clinical score pred\nN =',ncol(clin)-1),
                    paste('six states pred\nN =',ncol(six)-1),
                    paste('four states pred\nN =',ncol(four)-1)),
           unname(param_sets),
           ilabels=T,
           zcolor = RColorBrewer::brewer.pal(7,'Set1'),ilcs=2,sncs=1.25 )


png(paste0(outdir,'/Venn.',N_clinical,
          'clinical_params.png'),width=700,height=700)
venn::venn(snames=c(paste('clinical score pred\nN =',ncol(clin)-1),
                    paste('six states pred\nN =',ncol(six)-1),
                    paste('four states pred\nN =',ncol(four)-1)),
           unname(param_sets),
           ilabels=T,
           zcolor = RColorBrewer::brewer.pal(7,'Set1'),ilcs=2,sncs=1.25 )
dev.off()

pdf(paste0(outdir,'/Venn.',N_clinical,
          'clinical_params.pdf'),width=10.00,height=10.00)
venn::venn(snames=c(paste('clinical score pred\nN =',ncol(clin)-1),
                    paste('six states pred\nN =',ncol(six)-1),
                    paste('four states pred\nN =',ncol(four)-1)),
           unname(param_sets),
           ilabels=T,
           zcolor = RColorBrewer::brewer.pal(7,'Set1'),ilcs=2,sncs=1.25 )
dev.off()

################################################
library(psych);library(ggplot2);library(caret);library(randomForest)
dpi_plot='7'

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
#filter counts & BTMs so that gene is expressed in at least 2 conditions (on average)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]


mat2=read.csv('/home/oscar/Documents/sheep_megadata/30.6.21/dat_plots_filt/BTMs.final.csv',row.names = 1)
mat2[1:2,1:2]
#########################################



############################
setwd("/home/oscar/Documents/sheep_megadata/30.6.21/")
library(ggplot2)
library(plotly)
library(gridExtra)
library(RColorBrewer)
library(Amelia)
library(reshape2)
data=read.csv('combined/combined_sheet.csv',stringsAsFactors = F)
#filter out SMI13, non dpi 7 data and <80% complete variables
for(i in 3:ncol(data)){
  if(length(which(is.na(data[data$dpi==7&grepl('SLC',data$ID),i])))==7){
    data[data$dpi==7&grepl('SLC',data$ID),i]=data[data$dpi==21&grepl('SLC',data$ID),i]
  }
  if(length(which(is.na(data[data$dpi==7&grepl('SLI6',data$ID),i])))==7){
    for(j in grep('SLI6',data$ID,value=T)){
      data[data$dpi==7&grepl(j,data$ID),i]=data[which(
        !is.na(data[,i])&
          data$dpi>7&
          grepl(j,data$ID))[1]
        ,i]
    }
    data[data$dpi==7&!grepl('SLC',data$ID),i]=data[,i]
  }
}

data=data[data$dpi==7&!grepl('SLI13',data$ID),]

cbind(colnames(data),colSums(is.na(data))/ncol(data))


data=data[,(colSums(is.na(data))<0.2*nrow(data))]

cols=c("#20B2AA",brewer.pal(6,"Reds")[1+c(1,3,5,4,2)])


# pca_dat=amelia(data[,3:ncol(data)],logs=colnames(data)[3:ncol(data)])$imputations$imp1
# 
# pca_dat=amelia(data[,3:ncol(data)])$imputations$imp1

## make empty cells equal to the column mean
pca_dat=data[,3:ncol(data)]
for(i in 1:ncol(pca_dat)){
  nas=is.na(pca_dat[,i])
  nas[pca_dat[,i]==""]=T
  if(any(nas)){
    pca_dat[nas,i]=mean(as.numeric((pca_dat[,i])[!nas]))
  }
  pca_dat[,i]=as.numeric(pca_dat[,i])
}

rownames(pca_dat)=paste(data$ID,data$dpi,sep='_')
dim(pca_dat)

metabol_dat_save=pca_dat
colnames(pca_dat)
###################################################################

colnames(mat2)=  sapply(colnames(mat2),function(x) 
  gsub('\\.','-',paste(strsplit(x,'_')[[1]][c(1,3)],collapse='_')))
mat3=t(mat2)

pca_dat=pca_dat[rownames(pca_dat)%in%rownames(mat3),]
mat3=mat3[match(rownames(pca_dat),rownames(mat3)),]
#write.csv(pca_dat,'/home/oscar/Documents/sheep_megadata/13.4.21/dat_plots_filt/metabolic.final.csv')

################################################
#test line: comb_dat_funct_in=comb_dat;types_funct_convert=types_convert;xvar=50;clinkey=NA

gettop150=function(comb_dat_funct_in,types_funct_convert,clinkey=NA){
  comb_dat_funct=comb_dat_funct_in
  comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
  
  if(!is.na(clinkey[1])){
    types_funct2=clinkey
  }else{
    types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
    types_funct2=as.character(types_funct_convert[types_funct])
  }
  
  xvar2=150
  RF_1=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
  hits=rownames(RF_1$importance)[order(RF_1$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
  
  RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
  
  gini= RF$importance[,'MeanDecreaseGini']
  head(names(gini))
  top150=gini[order(gini,decreasing = T)[1:150]]
  head((top150))
  return(top150)
}
###################################################

library(smotefamily)
library(randomForest)
library(caret)
colnames(mat3)=paste('BTM:',colnames(mat3),sep='')
#drop highly correlated
#pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('Conteggio.Globuli.Rossi..RBC.','GOT.AST.')==x))]#'TNF.alpha'

comb_dat=as.data.frame(cbind(pca_dat,mat3))
comb_dat=comb_dat[,!grepl('rectal|Albumin',colnames(comb_dat))]
comb_dat=comb_dat[!grepl('SMI6-131|SMI6-132|SMI6-133',rownames(comb_dat)),]
colnames(comb_dat)=gsub("[[:blank:]]", "", colnames(comb_dat))
colnames(comb_dat)=gsub('\\.\\.','\\.',colnames(comb_dat))
colnames(comb_dat)=gsub('\\.\\.','\\.',colnames(comb_dat))
colnames(comb_dat)=gsub('\\)','',colnames(comb_dat))
library(smotefamily)
library(randomForest)



types=unlist(strsplit(rownames(comb_dat),'-'))[(1:nrow(data))*2-1]
types_convert=list()
types_convert[['SMC']]='control_Sass'
types_convert[['SLC']]='control_Ter'
types_convert[['SMI8']]='BTV8_Sass'
types_convert[['SMI13']]='2013_Sass'
types_convert[['SLI13']]='2013_Ter'
types_convert[['SMI6']]='2006_Sass'
types_convert[['SLI6']]='2006_Ter'



top150_6=gettop150(comb_dat,types_convert)
######################################

types=unlist(strsplit(rownames(comb_dat),'-'))[(1:nrow(data))*2-1]
types_convert=list()
types_convert[['SMI6']]=types_convert[['SLI6']]='2006'
types_convert[['SMI13']]=types_convert[['SLI13']]='2013'
types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8'
types_convert[['SMC']]=types_convert[['SLC']]='control'

top150_4=gettop150(comb_dat,types_convert)


##################################################

clinicals_in=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv')
clinicals_in=clinicals_in[clinicals_in$dpi==dpi_plot,]
#classes_N=30
animals=unlist(strsplit(rownames(comb_dat),'_'))[(1:nrow(comb_dat))*2-1]


clinicals_in=clinicals_in[clinicals_in$ID%in%animals,]
comb_dat=comb_dat[order(clinicals_in$clinical.score[match(animals,clinicals_in$ID)]),]
clinicals=clinicals_in$clinical.score[match(animals,clinicals_in$ID)]
names(clinicals)=animals # bad coding makes this necessary for function

convert_clinicals=list()
convert_clinicals[['0']]=convert_clinicals[['1']]=convert_clinicals[['2']]='inf clin_score 0-2'
convert_clinicals[['3']]=convert_clinicals[['4']]=convert_clinicals[['5']]='inf clin_score 3-5'
convert_clinicals[['6']]=convert_clinicals[['7']]=convert_clinicals[['8']]='inf clin_score 6-8'
clinicals_discrete=rep('control',length(animals))
clinicals_discrete[!grepl('SMC|SLC',animals)]=
  unlist(convert_clinicals[as.character(clinicals[!grepl('SMC|SLC',animals)])])

top150_clin=gettop150(comb_dat,types_convert,clinicals_discrete)


##
shared_params=gsub('BTM\\.','BTM:',shared_params)
shared_params=gsub('combine\\.\\.','combine(',shared_params)
shared_params=gsub('\\.\\.','\\.',shared_params)
shared_params=gsub('\\.\\.','\\.',shared_params)
shared_params=substr(shared_params,1,50)

#######################################
#######################################
#######################################
#######################################
ggdata=data.frame(parameters=shared_params)


ggdata$parameters%in%names(top150_6)
ggdata$parameters
names(top150_6)



ggdata$six_states=top150_6[match(ggdata$parameters,substr(names(top150_6),1,50))]
ggdata$four_states=top150_4[match(ggdata$parameters,substr(names(top150_4),1,50))]
ggdata$clinical=top150_clin[match(ggdata$parameters,substr(names(top150_clin),1,50))]
ggdata$parameters[7]
match("BTM:combine(TBA.M137.TBA.M180..",substr(names(top150_clin),1,50))
grep('TBA.M137',substr(names(top150_clin),1,50),value=T)
############
new=melt(ggdata,id='parameters',value.name = 'importance')
grep('enriched',names(top150_4),value=T)
unique(new$parameters[is.na(new$importance)])
dim(new)

new$parameters=substr(new$parameters,1,30)
ggdata$sum=ggdata$clinical+ggdata$four_states+ggdata$six_states
new$parameters=factor(new$parameters,levels=new$parameters[order(ggdata$sum,decreasing = T)])

png(paste0(outdir,'/importance_fully_shared_params_using',N_clinical,
          '_clinical_params.png'),width=500,height=350)
ggplot(new,aes(y=importance,fill=variable,x=parameters))+
  geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
  ylab('gini importance')+
  ggtitle(paste('importance of',length(unique(new$parameters)),'parameters shared across all three'))
dev.off()


pdf(paste0(outdir,'/importance_fully_shared_params_using',N_clinical,
          '_clinical_params.pdf'),width=7.00,height=6.00)
ggplot(new,aes(y=importance,fill=variable,x=parameters))+
  geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
  ylab('gini importance')+
  ggtitle(paste('importance of',length(unique(new$parameters)),'parameters shared across all three'))
dev.off()


##
shared_once_params=gsub('BTM\\.','BTM:',shared_once_params)
shared_once_params=gsub('combine\\.\\.','combine(',shared_once_params)
shared_once_params=gsub('\\.\\.','\\.',shared_once_params)
shared_once_params=gsub('\\.\\.','\\.',shared_once_params)
shared_once_params=gsub('\\)','',shared_once_params)
shared_once_params=substr(shared_once_params,1,50)

#######################################
#######################################
#######################################
#######################################
ggdata=data.frame(parameters=shared_once_params)


ggdata$parameters%in%names(top150_6)
ggdata$parameters
names(top150_6)



ggdata$six_states=top150_6[match(ggdata$parameters,substr(names(top150_6),1,50))]
ggdata$four_states=top150_4[match(ggdata$parameters,substr(names(top150_4),1,50))]
ggdata$clinical=top150_clin[match(ggdata$parameters,substr(names(top150_clin),1,50))]
############
new=melt(ggdata,id='parameters',value.name = 'importance')

unique(new$parameters[is.na(new$importance)])
dim(new)

new$parameters=substr(new$parameters,1,30)
ggdata$sum=rep(0,nrow(ggdata))
for(i in 1:nrow(ggdata)){
  if(!is.na(ggdata$clinical[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$clinical[i]}
  if(!is.na(ggdata$six_states[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$six_states[i]}
  if(!is.na(ggdata$four_states[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$four_states[i]}
  
}




new$parameters=factor(new$parameters,levels=substr(ggdata$parameters,1,30)[order(ggdata$sum,decreasing = T)])

ggdata[order(ggdata$sum,decreasing = T),]


png(paste0(outdir,'/importance_any_pair_shared_params_using',N_clinical,
          'clinical_params.png'),width=700,height=400)
ggplot(new,aes(y=importance,fill=variable,x=parameters))+
  geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
  ylab('gini importance')+
  ggtitle(paste('importance of',length(unique(new$parameters)),'parameters shared in any combination'))
dev.off()

pdf(paste0(outdir,'/importance_any_pair_shared_params_using',N_clinical,
          'clinical_params.pdf'),width=10.00,height=7.00)

ggplot(new,aes(y=importance,fill=variable,x=parameters))+
  geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
  ylab('gini importance')+
  ggtitle(paste('importance of',length(unique(new$parameters)),'parameters shared in any combination'))
dev.off()

# 
# 
# 
# new2=new[new$variable=='four_states',]
# 
# 
# new2$parameters=factor(new2$parameters,levels=new2$parameters[rev(order(new2$importance))])
# 
# new2=new2[!is.na(new2$importance),]
# 
# ggplot(head(new2,10),aes(y=importance,fill=variable,x=parameters))+
#   geom_bar(stat='identity')+theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
#   ylab('gini importance')+
#   ggtitle('importance of 40 shared in any combo parameters')
# 
# 
