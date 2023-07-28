library(psych);library(ggplot2);library(caret);library(randomForest)
dpi_plot='7'

runs_of_RF=500;classes_N=30
shuffle=FALSE

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
#filter counts & BTMs so that gene is expressed in at least 2 conditions (on average)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]


mat2=read.csv('/home/oscar/Documents/sheep_megadata/14.9.21/dat_plots_filt/BTMs.final.csv',row.names = 1)
mat2[1:2,1:2]
#########################################

outdir="/home/oscar/scripts/github/sheep_ML/outdir"

############################
setwd("/home/oscar/Documents/sheep_megadata/14.9.21/")

outdir=gsub("/$","",outdir)
library(ggplot2)
library(plotly)
library(gridExtra)
library(RColorBrewer)
library(Amelia)
library(reshape2)
library(wesanderson)
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

data=data[data$dpi==dpi_plot&!grepl('SLI13',data$ID),]

data=data[,(colSums(is.na(data))<0.19*nrow(data))]

cols=c("#20B2AA",brewer.pal(6,"Reds")[1+c(1,3,5,4,2)])


#if data to be shuffled to show performance of RF not random:

  if(shuffle==T){
      out_file_extra="SHUFFLED_"
  }else{out_file_extra=""}


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
mat3=mat3[sapply(rownames(mat3),function(x) strsplit(x,'_')[[1]][2])==dpi_plot,]
pca_dat=pca_dat[rownames(pca_dat)%in%rownames(mat3),]
mat3=mat3[match(rownames(pca_dat),rownames(mat3)),]
# write.csv(pca_dat,'/home/oscar/Documents/sheep_megadata/13.4.21/dat_plots_filt/metabolic.final.csv')

################################################
###################################################

library(smotefamily)
library(randomForest)
library(caret)
colnames(mat3)=paste('BTM:',colnames(mat3),sep='')

comb_dat=as.data.frame(cbind(pca_dat,mat3))
comb_dat=comb_dat[,!grepl('rectal|Albumin',colnames(comb_dat))]
comb_dat=comb_dat[!grepl('SMI6-131|SMI6-132|SMI6-133',rownames(comb_dat)),]



#drop highly correlated
pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('MIP.1.alpha','IL.1alpha','GOT.AST.')==x))]

library(smotefamily)
library(randomForest)


RUN_FILTRATION_AND_PREDICTION_clin_imbalanced=function(comb_dat_funct_in,types_funct_in,runs_of_RF_funct,xvar){
  comb_dat_funct=comb_dat_funct_in
  comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
  types_funct2=types_funct_in
  xvar2=150
  if(shuffle==TRUE){
      comb_dat_funct=comb_dat_funct[sample(1:nrow(comb_dat_funct),replace = FALSE),]
  }else{out_file_extra=""}

  RF_1=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
  hits=rownames(RF_1$importance)[order(RF_1$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
  
  rfcont=rfeControl(functions=rfFuncs,repeats=5,verbose=F,rerank=F,method="repeatedcv",number=5)
  RFE=rfe(x=comb_dat_funct[,hits], y=as.factor(types_funct2),rfeControl=rfcont,sizes=c(1:10,seq(12,xvar2,2)))
  
  ##

  plot_loop_fun=function(){
    par(mfrow=c(1,2),mar=c(5,4,1,1))
  plot(RFE$results$Variables,RFE$results$Kappa,ylim=c(0,1),xlab='parameters',ylab='prediction accurracy',
       col=scales::alpha(wesanderson::wes_palette('Darjeeling1')[1],.6),pch=19)
  par(new=T)
  plot(RFE$results$Variables,RFE$results$Accuracy,ylim=c(0,1),
       col=wesanderson::wes_palette('Darjeeling1')[2],type='l',xlab='',ylab='')
  abline(v=RFE$bestSubset,col=wesanderson::wes_palette('Darjeeling1')[3],lty=2)
  if(xvar!=0){
    abline(v=xvar,col=wesanderson::wes_palette('Darjeeling1')[5],lty=2)
  }
  plot(RFE$results$Variables,RFE$results$Kappa,ylim=c(0.35,.8),xlab='',ylab='',
       col=scales::alpha(wesanderson::wes_palette('Darjeeling1')[1],.6),pch=19)
  par(new=T)
  plot(RFE$results$Variables,RFE$results$Accuracy,ylim=c(0.35,.8),
       col=wesanderson::wes_palette('Darjeeling1')[2],type='l',
       xlab='parameters',ylab='prediction accuracy ZOOM')
  abline(v=RFE$bestSubset,col=wesanderson::wes_palette('Darjeeling1')[3],lty=2)
  if(xvar!=0){
    abline(v=xvar,col=wesanderson::wes_palette('Darjeeling1')[5],lty=2)
  }
  if(xvar==0){
    legend(legend=c('kappa','accuracy','"best" subset'),col=wesanderson::wes_palette('Darjeeling1')[c(1:3)],
           lwd=2,bty='n',x='topleft',lty=c(1,1,2))
  }
  else{
    legend(legend=c('kappa','accuracy','"best" subset #','used subset #' ),
           col=wesanderson::wes_palette('Darjeeling1')[c(1:3,5)],lwd=2,bty='n',x='topleft',
           lty=c(1,1,2,2))
  }
  }
 
  png(paste0(outdir,'/plots/clinical_',out_file_extra,
            paste(unique(types_funct2),collapse='_'),'.png'),width=700,height=400)
  plot_loop_fun()
  dev.off()
  ##
  pdf(paste0(outdir,'/plots/clinical_',out_file_extra,
            paste(unique(types_funct2),collapse='_'),'.pdf'),width=12.0,height=7.0)
  
  if(xvar==0){
    xvar=RFE$optsize
  }
  scores=matrix(nrow=3,ncol=0)
  param_hits=c()
  ###############################################################
  comb_dat_funct=comb_dat_funct_in
  comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
  
  types_funct2=types_funct_in
  RF=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F)
  hits_xvar2=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
  
  ##
  for(test in 1:runs_of_RF){
    comb_dat_funct=comb_dat_funct_in
    comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
    #if set to shuffle mode rotate sample order data each loop
  if(shuffle==TRUE){
      comb_dat_funct=comb_dat_funct[sample(1:nrow(comb_dat_funct),replace = FALSE),]
  }


    types_funct2=types_funct_in
    #drop highly correlated
    #pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('MIP.1.alpha','IL.1alpha','GOT.AST.')==x))]
    
    library(caret);library(randomForest)
    
    holdouts=c()
    for(i in unique(types_funct2)){
      holdouts=c(holdouts,sample(which(types_funct2==i),2))
    }
    comb_dat_funct_test=comb_dat_funct[holdouts,]
    comb_dat_funct=comb_dat_funct[-holdouts,]
    types_funct2_test=types_funct2[holdouts]
    types_funct2=types_funct2[-holdouts]
    
    
    
    
    if(length(unique(table(types_funct2)))>1){
      
      ########### add synthetic
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=3,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    if(length(unique(table(types_funct2)))>1){
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=3,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F)
    RF$confusion
    ######################################
    hits=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar]
    
    param_hits=c(param_hits,hits)

    RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F)
    scores=cbind(scores,rbind(names(clinicals)[holdouts],types_funct2_test,
                              as.character(predict(RF,newdata=comb_dat_funct_test))))
  }
  # hits_final_funct=(table(param_hits)[table(param_hits)>0.90*max(table(param_hits))])[
  #   order(as.numeric((table(param_hits)[table(param_hits)>0.90*max(table(param_hits))])),decreasing=T)]
  hits_final_funct=table(param_hits)[order(table(param_hits),decreasing=T)[1:xvar]]
  
  #############################################################################################################################################
  #############################################################################################################################################
  scores=matrix(nrow=3,ncol=0)
  
  for(test in 1:runs_of_RF){
    comb_dat_funct=comb_dat_funct_in
    comb_dat_funct=comb_dat_funct[,names(hits_final_funct)]
    
    types_funct2=types_funct_in
    
    library(caret);library(randomForest)
    holdouts=c()
    for(i in unique(types_funct2)){
      holdouts=c(holdouts,sample(which(types_funct2==i),2))
    }
    comb_dat_funct_test=comb_dat_funct[holdouts,]
    comb_dat_funct=comb_dat_funct[-holdouts,]
    types_funct2_test=types_funct2[holdouts]
    types_funct2=types_funct2[-holdouts]
    
    if(length(unique(table(types_funct2)))>1){
      ########### add synthetic
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=3,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    if(length(unique(table(types_funct2)))>1){
      ########### add synthetic
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=3,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
        # generate new random forest model each loop with new holdouts
    RF=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F)
    scores=cbind(scores,rbind(names(clinicals)[holdouts],
                              types_funct2_test,as.character(predict(RF,newdata=comb_dat_funct_test))))
  }
  
  
  confusion=matrix(nrow=length(unique(types_funct2)),ncol=length(unique(types_funct2)),dat=0)
  rownames(confusion)=colnames(confusion)=score_vals=unique(unique(scores[2,]))
  for(scores_i in 1:ncol(confusion)){
    test_i=which(scores[2,]==score_vals[scores_i])
    for(scores_j in 1:ncol(confusion)){
      confusion[scores_i,scores_j]=length(which(scores[3,test_i]==score_vals[scores_j]))
    }
  }
  return(list(table(param_hits),confusion,xvar,RF_1$importance[,ncol(RF_1$importance)],RFE$results$Kappa,scores,names(hits_final_funct)))
}



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
groups_hierarchy=c("control","inf clin_score 0-2","inf clin_score 3-5","inf clin_score 6-8")
clinicals_discrete=rep('control',length(animals))
clinicals_discrete[!grepl('SMC|SLC',animals)]=
  unlist(convert_clinicals[as.character(clinicals[!grepl('SMC|SLC',animals)])])

comb_dat_funct_in=bin1=comb_dat;types_funct_in=clinicals_discrete
############## RUN it!

#run with defined optimal parameter numbers first
funct_outclinclass=RUN_FILTRATION_AND_PREDICTION_clin_imbalanced(bin1,clinicals_discrete,runs_of_RF,classes_N)

# bin2=bin1[!grepl('SMC|SLC',rownames(bin1)),];clinicals_discrete2=clinicals_discrete[!grepl('SMC|SLC',rownames(bin1))]
# #run again without SMC and SLC, with no set number of parameters
# funct_outclinclass2=RUN_FILTRATION_AND_PREDICTION_clin_imbalanced(bin2,clinicals_discrete2,runs_of_RF,0)
# print(funct_outclinclass2[[3]])
# substr(names(funct_outclinclass2[[1]])[order(funct_outclinclass2[[1]],decreasing=T)[1:30]],1,30)

plot((funct_outclinclass[[5]]))

print(funct_outclinclass[[2]])

confusion_out=funct_outclinclass[[2]]
confusion_out=confusion_out[,groups_hierarchy]
confusion_out=confusion_out[groups_hierarchy,]
write.csv(confusion_out,file=paste0(outdir,'/clinical_RF_performance_table_',out_file_extra,classes_N,'params.csv'))

knitr::kable(confusion_out)
#return(list(table(param_hits),confusion,xvar,RF_1$importance[,ncol(RF_1$importance)],RFE$results$Kappa,scores))

par(mfrow=c(1,1))
plot((funct_outclinclass[[1]]))
print(funct_outclinclass[[2]])

substr(names(funct_outclinclass[[1]])[order(funct_outclinclass[[1]],decreasing=T)[1:print(funct_outclinclass[[3]])]],1,30)

















if(shuffle==TRUE){stop("no point making heatmaps for shuffled datasets")}


################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
################ heatmap#####################################
data=read.csv('combined/combined_sheet.csv',stringsAsFactors = F)
#filter out SMI13, non dpi 7 data and <80% complete variables
for(i in 3:ncol(data)){
  if(length(which(is.na(data[data$dpi==dpi_plot&grepl('SLC',data$ID),i])))==7){
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

data=data[data$dpi==dpi_plot,]
rownames(data)=data$ID
pca_dat=data[,3:ncol(data)]


mat3=t(mat2)
mat3=mat3[sapply(rownames(mat3),function(x) strsplit(x,'_')[[1]][2])==dpi_plot,]

rownames(mat3)=sapply(rownames(mat3),function(x) strsplit(x,'_')[[1]][1])
pca_dat=pca_dat[rownames(pca_dat)%in%rownames(mat3),]
mat3=mat3[match(rownames(pca_dat),rownames(mat3)),]
colnames(mat3)=paste('BTM:',colnames(mat3),sep='')

###################################################################
comb_dat=as.data.frame(cbind(pca_dat,mat3))
comb_dat=comb_dat[!grepl('SMI6-131|SMI6-132|SMI6-133',rownames(comb_dat)),]

clinicals_in=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv')
clinicals_in=clinicals_in[clinicals_in$dpi==dpi_plot,]
clinicals_in=clinicals_in[match(rownames(comb_dat),clinicals_in$ID),]
clinicals_plot=clinicals_in$clinical.score[match(rownames(comb_dat),clinicals_in$ID)]

comb_dat=comb_dat[order(clinicals_plot),]

comb_dat=comb_dat[c(grep('SLC|SMC',rownames(comb_dat)),grep('SLC|SMC',rownames(comb_dat),invert=T)),]
clinicals_plot=clinicals_in$clinical.score[match(rownames(comb_dat),clinicals_in$ID)]

plot(clinicals_plot)
# 
# comb_dat=comb_dat[,c('IFN.y','Proteine.Totali')]
# 
# 

topotopparams= names(funct_outclinclass[[1]])[order(funct_outclinclass[[1]],decreasing=T)[1:funct_outclinclass[[3]]]]
# topotopparams= names(funct_outclinclass2[[1]])[order(funct_outclinclass2[[1]],decreasing=T)[1:30]]

write.csv(comb_dat,file=paste0(outdir,'/all_RF.data_clin',out_file_extra,'.csv'))
write.csv(comb_dat[,topotopparams],file=paste0(outdir,'/clinical_score.data_',out_file_extra,funct_outclinclass[[3]],
                                              '.params.csv'))


#comb_dat=read.csv(file='/home/oscar/Documents/sheep_megadata/14.9.21/RF_input_data/clinical_score.data.csv',row.names=TRUE)

#pheaty=apply(comb_dat,2,function(x) (log(x+1) - mean(log(x[grepl('SMC|SLC',rownames(comb_dat))]+1)))/sd(log(x[!is.na(x)]+1)))[,topotopparams]
pheaty=apply(comb_dat[,topotopparams],2,function(x) ((x) - mean((x[grepl('SMC|SLC',rownames(comb_dat))])))/sd((x[!is.na(x)])))#[,topotopparams]
rownames(pheaty)
storage=rownames(pheaty)
colnames(pheaty)=substr(colnames(pheaty),1,32)
#rownames(pheaty)=paste(rownames(pheaty),clinicals_plot,sep=' #')
clinicals_pheaty_plot=clinicals_in$clinical.score[match(storage,clinicals_in$ID)]

rownames(pheaty)[grep('SMC',rownames(pheaty))]='Mock S'
rownames(pheaty)[grep('SLC',rownames(pheaty))]='Mock T'
rownames(pheaty)[grep('SMI6',rownames(pheaty))]='BTV-1 2006 S'
rownames(pheaty)[grep('SLI6',rownames(pheaty))]='BTV-1 2006 T'
rownames(pheaty)[grep('SLI13',rownames(pheaty))]='BTV-1 2013 S'
rownames(pheaty)[grep('SLI13',rownames(pheaty))]='BTV-1 2013 T'
rownames(pheaty)[grep('SMI8',rownames(pheaty))]='BTV-8 S'
rownames(pheaty)[grep('SLI8',rownames(pheaty))]='BTV-8 T'

rownames(pheaty)=paste(rownames(pheaty),' CS:',clinicals_pheaty_plot,sep='')

colnames(pheaty)=substr(colnames(pheaty),1,32)
library(grid);library(pheatmap)

pal=c(colorRampPalette(c('#1842a5',"#3A5DB1",RColorBrewer::brewer.pal(7,'RdBu')[c(4,2)],
                         "#C4393B","#B2182B"
))(60))

pal=pal[6:60]
pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= pal,na_col='#000000')

grid.abline(194.3,0,range='x')
#grid.abline(255.5,0,range='x')
#grid.abline(282.5,0,range='x')
grid.abline(315.3,0,range='x')
#grid.abline(347.5,0,range='x')
#grid.abline(452.5,0,range='x')
#grid.abline(531.5,0,range='x')
grid.abline(530.3,0,range='x')
#grid.abline(649.5,0,range='x')



#########################################################################################
#export at 900 height
# 
# #pheaty=apply(comb_dat,2,function(x) (log(x+1) - mean(log(x[grepl('SMC|SLC',rownames(comb_dat))]+1)))/sd(log(x[!is.na(x)]+1)))[,topotopparams]

# rownames(pheaty)
# 
# colnames(pheaty)=substr(colnames(pheaty),1,32)
# rownames(pheaty)=paste(rownames(pheaty),clinicals_plot,sep=' #')
# library(grid);library(pheatmap)
# 
# pal=c(colorRampPalette(c('#1842a5',"#3A5DB1",RColorBrewer::brewer.pal(7,'RdBu')[c(4,2)],
#                          "#C4393B","#B2182B"
# ))(60))
# 
# pal=pal[6:60]
# pheatmap((pheaty),
#          cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
#          color= pal,na_col='#000000')
# 
# grid.abline(194.3,0,range='x')
# #grid.abline(255.5,0,range='x')
# #grid.abline(282.5,0,range='x')
# grid.abline(315.3,0,range='x')
# #grid.abline(347.5,0,range='x')
# #grid.abline(452.5,0,range='x')
# #grid.abline(531.5,0,range='x')
# grid.abline(530.3,0,range='x')
# #grid.abline(649.5,0,range='x')
# 
# #export at 900 height

dev.off()

png(paste0(outdir,'/plots/clinical_pheatmap_',classes_N,'_params_',
          paste(unique(clinicals_discrete),collapse='_'),'.png'),width=1400,height=850)

pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= pal,na_col='#000000')

grid.abline(194.3,0,range='x')
#grid.abline(255.5,0,range='x')
#grid.abline(282.5,0,range='x')
grid.abline(316.3,0,range='x')
#grid.abline(347.5,0,range='x')
#grid.abline(452.5,0,range='x')
#grid.abline(531.5,0,range='x')
grid.abline(534.3,0,range='x')
#grid.abline(649.5,0,range='x')

#export at 900 height
dev.off()



png(paste0(outdir,'/plots/clinical_pheatmap_no_lines_',classes_N,'_params_',
          paste(unique(clinicals_discrete),collapse='_'),'.png'),width=1400,height=850)

pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= pal,na_col='#000000')
# 
# grid.abline(194.3,0,range='x')
# #grid.abline(255.5,0,range='x')
# #grid.abline(282.5,0,range='x')
# grid.abline(316.3,0,range='x')
# #grid.abline(347.5,0,range='x')
# #grid.abline(452.5,0,range='x')
# #grid.abline(531.5,0,range='x')
# grid.abline(534.3,0,range='x')
# #grid.abline(649.5,0,range='x')

#export at 900 height
dev.off()



pdf(paste0(outdir,'/plots/clinical_pheatmap_',classes_N,'_params_',
          paste(unique(clinicals_discrete),collapse='_'),'.pdf'),width=14.00,height=9.00)


pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= pal,na_col='#000000')

grid.abline(272.3,0,range='x')
#grid.abline(255.5,0,range='x')
#grid.abline(282.5,0,range='x')
grid.abline(424.6,0,range='x')
#grid.abline(347.5,0,range='x')
#grid.abline(452.5,0,range='x')
#grid.abline(531.5,0,range='x')
grid.abline(510.3,0,range='x')
#grid.abline(649.5,0,range='x')
#export at 900 height
dev.off()


pdf(paste0(outdir,'/plots/clinical_pheatmap_no_lines_',classes_N,'_params_',
          paste(unique(clinicals_discrete),collapse='_'),'.pdf'),width=14.00,height=9.00)

pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= pal,na_col='#000000')

# grid.abline(272.3,0,range='x')
# #grid.abline(255.5,0,range='x')
# #grid.abline(282.5,0,range='x')
# grid.abline(424.6,0,range='x')
# #grid.abline(347.5,0,range='x')
# #grid.abline(452.5,0,range='x')
# #grid.abline(531.5,0,range='x')
# grid.abline(510.3,0,range='x')
# #grid.abline(649.5,0,range='x')

#export at 900 height
dev.off()

