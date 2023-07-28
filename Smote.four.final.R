library(psych);library(ggplot2);library(caret);library(randomForest)
dpi_plot='7'

#N_classes gives number of parameters in RF model
#runs_of_RF gives number of RF models to run for confusion matrix
runs_of_RF=500;N_classes=17

shuffle=FALSE

outdir="/home/oscar/scripts/github/sheep_ML/outdir"
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
#filter counts & BTMs so that gene is expressed in at least 2 conditions (on average)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]


mat2=read.csv('/home/oscar/Documents/sheep_megadata/14.9.21/dat_plots_filt/BTMs.final.csv',row.names = 1)
mat2[1:2,1:2]
#########################################



############################
setwd("/home/oscar/Documents/sheep_megadata/14.9.21/")
library(ggplot2)
library(plotly)
library(gridExtra)
library(RColorBrewer)
library(Amelia)
library(reshape2)
library(caret)
library(randomForest)
library(smotefamily)
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


outdir=gsub("/$","",outdir)
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

pca_dat=pca_dat[rownames(pca_dat)%in%rownames(mat3),]
mat3=mat3[match(rownames(pca_dat),rownames(mat3)),]
#write.csv(pca_dat,'/home/oscar/Documents/sheep_megadata/13.4.21/dat_plots_filt/metabolic.final.csv')

################################################
###################################################


colnames(mat3)=paste('BTM:',colnames(mat3),sep='')
comb_dat=as.data.frame(cbind(pca_dat,mat3))
comb_dat=comb_dat[,!grepl('rectal|Albumin',colnames(comb_dat))]
comb_dat=comb_dat[!grepl('SMI6-131|SMI6-132|SMI6-133',rownames(comb_dat)),]



#drop highly correlated
pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('Conteggio.Globuli.Rossi..RBC.','GOT.AST.')==x))]#'TNF.alpha'

library(smotefamily)
library(randomForest)



#comb_dat_funct_in=bin1;types_funct_convert=types_convert
# RUN_FILTRATION_AND_PREDICTION=function(comb_dat_funct_in,types_funct_convert,runs_of_RF_funct,xvar){
#   comb_dat_funct=comb_dat_funct_in
#   comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
#   types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
#   types_funct2=as.character(types_funct_convert[types_funct])
#   xvar2=150
#   RF_1=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
#   hits=rownames(RF_1$importance)[order(RF_1$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
#   xvar2=50
#   RF_1=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
#   hits=rownames(RF_1$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
#   
#   rfcont=rfeControl(functions=rfFuncs,repeats=5,verbose=F,rerank=F,method="repeatedcv",number=5)
#   RFE=rfe(x=comb_dat_funct[,hits], y=as.factor(types_funct2),rfeControl=rfcont,sizes=c(1:10,seq(12,xvar2,2)))
#   
#   ##
#   png(paste('/home/oscar/Pictures/plots/Sheep_megadata/RF_11.5.21/binary_classifier',substr(types_funct[1],1,2),
#             paste(unique(types_funct2),collapse='_'),'.png',sep='_'),width=700,height=400)
#   par(mfrow=c(1,2),mar=c(5,4,1,1))
#   plot(RFE$results$Variables,RFE$results$Kappa,ylim=c(0,1),xlab='parameters',ylab='prediction accur')
#   par(new=T)
#   plot(RFE$results$Variables,RFE$results$Accuracy,ylim=c(0,1),col=2,type='l',xlab='',ylab='')
#   abline(v=RFE$bestSubset,col=3)
#   plot(RFE$results$Variables,RFE$results$Kappa,ylim=c(0.8,1),xlab='',ylab='')
#   par(new=T)
#   plot(RFE$results$Variables,RFE$results$Accuracy,ylim=c(0.8,1),col=2,type='l',xlab='parameters',ylab='prediction accur ZOOM')
#   abline(v=RFE$bestSubset,col=3)
#   legend(legend=c('kappa','accuracy','"best" subset'),col=1:3,lwd=2,bty='n',x='bottom')
#   dev.off()
#   if(xvar==0){
#     xvar=RFE$optsize
#   }
#   scores=matrix(nrow=2,ncol=0)
#   param_hits=c()
#   ###############################################################
#   for(test in 1:runs_of_RF){
#     comb_dat_funct=comb_dat_funct_in
#     comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
#     types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
#     types_funct2=as.character(types_funct_convert[types_funct])
#     #drop highly correlated
#     #pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('MIP.1.alpha','IL.1alpha','GOT.AST.')==x))]
#     
#     library(caret);library(randomForest)
#     
#     holdouts=c()
#     for(i in unique(types_funct2)){
#       holdouts=c(holdouts,sample(which(types_funct2==i),2))
#     }
#     comb_dat_funct_test=comb_dat_funct[holdouts,]
#     comb_dat_funct=comb_dat_funct[-holdouts,]
#     types_funct2_test=types_funct2[holdouts]
#     types_funct2=types_funct2[-holdouts]
#     
#     
#     #
#     
#     RF=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F)
#     hits=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
#     RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F)
#     RF$confusion
#     
#     hits=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar]
#     
#     param_hits=c(param_hits,hits)
#     # if(length(unique(table(types_funct2)))>1){
#     # 
#     # ########### add synthetic
#     #   smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=3,dup_size=2)
#     #   smote_test$syn_data$class
#     #   comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
#     #   types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
#     # }
#     # if(length(unique(table(types_funct2)))>1){
#     #   smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=4,dup_size=2)
#     #   smote_test$syn_data$class
#     #   comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
#     #   types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
#     # }
#     #######################################
#     
#     RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F)
#     scores=cbind(scores,rbind(types_funct2_test,as.character(predict(RF,newdata=comb_dat_funct_test))))
#   }
#   
#   
#   hits_final_funct=(table(param_hits)[table(param_hits)>0.90*max(table(param_hits))])[
#     order(as.numeric((table(param_hits)[table(param_hits)>0.90*max(table(param_hits))])),decreasing=T)]
#   
#   
#   
#   #############################################################################################################################################
#   #############################################################################################################################################
#   scores=matrix(nrow=2,ncol=0)
#   
#   for(test in 1:runs_of_RF){
#     comb_dat_funct=comb_dat_funct_in
#     comb_dat_funct=comb_dat_funct[,names(hits_final_funct)]
#     types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
#     types_funct2=as.character(types_funct_convert[types_funct])
#     #drop highly correlated
#     #pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('MIP.1.alpha','IL.1alpha','GOT.AST.')==x))]
#     
#     library(caret);library(randomForest)
#     holdouts=c()
#     for(i in unique(types_funct2)){
#       holdouts=c(holdouts,sample(which(types_funct2==i),2))
#     }
#     comb_dat_funct_test=comb_dat_funct[holdouts,]
#     comb_dat_funct=comb_dat_funct[-holdouts,]
#     types_funct2_test=types_funct2[holdouts]
#     types_funct2=types_funct2[-holdouts]
#     
#     # if(length(unique(table(types_funct2)))>1){
#     #   ########### add synthetic
#     # smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=4,dup_size=2)
#     # smote_test$syn_data$class
#     # comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
#     # types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
#     # }
#     
#     RF=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F)
#     scores=cbind(scores,rbind(types_funct2_test,as.character(predict(RF,newdata=comb_dat_funct_test))))
#   }
#   
#   
#   confusion=matrix(nrow=length(unique(types_funct2)),ncol=length(unique(types_funct2)),dat=0)
#   rownames(confusion)=colnames(confusion)=score_vals=unique(unique(scores[1,]))
#   for(scores_i in 1:ncol(confusion)){
#     test_i=which(scores[1,]==score_vals[scores_i])
#     for(scores_j in 1:ncol(confusion)){
#       confusion[scores_i,scores_j]=length(which(scores[2,test_i]==score_vals[scores_j]))
#     }
#   }
#   return(list(table(param_hits),confusion,xvar,RF_1$importance[,ncol(RF_1$importance)],RFE$results$Kappa))
# }

RUN_FILTRATION_AND_PREDICTION_imbalanced=function(comb_dat_funct_in,types_funct_convert,runs_of_RF_funct,xvar){
  comb_dat_funct=comb_dat_funct_in
  comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
  types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
  types_funct2=as.character(types_funct_convert[types_funct])
  
  xvar2=150

  if(shuffle==TRUE){
      comb_dat_funct=comb_dat_funct[sample(1:nrow(comb_dat_funct),replace = FALSE),]
  }


  RF_1=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
  hits=rownames(RF_1$importance)[order(RF_1$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
  
  rfcont=rfeControl(functions=rfFuncs,repeats=5,verbose=F,rerank=F,method="repeatedcv",number=5)
  #rfe with top 10 params then every other param up to xvar2 default 150 params.
  RFE=rfe(x=comb_dat_funct[,hits], y=as.factor(types_funct2),rfeControl=rfcont,sizes=c(1:10,seq(12,xvar2,2)))
  
  plot_loop_fun=function(){
    par(mfrow=c(1,2),mar=c(5,4,1,1))
    plot(RFE$results$Variables,RFE$results$Kappa,ylim=c(0,1),xlab='parameters',ylab='prediction accurracy',
        col=scales::alpha(wesanderson::wes_palette('Darjeeling1')[1],.6),pch=19)
    par(new=T)
    plot(RFE$results$Variables,RFE$results$Accuracy,ylim=c(0,1),col=wesanderson::wes_palette('Darjeeling1')[2],type='l',xlab='',ylab='')
    abline(v=RFE$bestSubset,col=wesanderson::wes_palette('Darjeeling1')[3],lty=2)
    if(xvar!=0){
      abline(v=xvar,col=wesanderson::wes_palette('Darjeeling1')[5],lty=2)
    }
    plot(RFE$results$Variables,RFE$results$Kappa,ylim=c(0.8,1),xlab='',ylab='',
        col=scales::alpha(wesanderson::wes_palette('Darjeeling1')[1],.6),pch=19)
    par(new=T)
    plot(RFE$results$Variables,RFE$results$Accuracy,ylim=c(0.8,1),col=wesanderson::wes_palette('Darjeeling1')[2],type='l',xlab='parameters',ylab='prediction accuracy ZOOM')
    abline(v=RFE$bestSubset,col=wesanderson::wes_palette('Darjeeling1')[3],lty=2)
    if(xvar!=0){
      abline(v=xvar,col=wesanderson::wes_palette('Darjeeling1')[5],lty=2)
    }
    if(xvar==0){
      legend(legend=c('kappa','accuracy','"best" subset'),col=wesanderson::wes_palette('Darjeeling1')[1:3],lwd=2,bty='n',x='bottom',lty=c(1,1,2))
    }
    else{
      legend(legend=c('kappa','accuracy','"best" subset #','used subset #' ),col=wesanderson::wes_palette('Darjeeling1')[c(1:3,5)],
            lwd=2,bty='n',x='topleft',lty=c(1,1,2,2))
    }
  }

  ##
  png(paste0(outdir,'/plots/four_states_',out_file_extra,
            paste(unique(types_funct2),collapse='_'),'.png'),width=700,height=400)
    plot_loop_fun()
  dev.off()
  ##
  pdf(paste0(outdir,'/plots/four_states_',out_file_extra,
            paste(unique(types_funct2),collapse='_'),'.pdf'),width=12.0,height=7.0)
    plot_loop_fun()
  dev.off()
  if(xvar==0){
    xvar=RFE$optsize
  }
  scores=matrix(nrow=2,ncol=0)
  param_hits=c()
  #do first round to get top 150 or whatever xvar2 is set to params
  types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
  types_funct2=as.character(types_funct_convert[types_funct])
  RF=randomForest(x=comb_dat_funct_in[,!grepl('rectal|Albumin',colnames(comb_dat_funct))], y=as.factor(types_funct2),importance=T,do.trace = F)
  hits_xvar2=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
  
  
  #RUN lots of random forests to find top performing params
  ###############################################################
  for(test in 1:runs_of_RF){
    comb_dat_funct=comb_dat_funct_in
    comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
    types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
    types_funct2=as.character(types_funct_convert[types_funct])
    #drop highly correlated
    #pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('MIP.1.alpha','IL.1alpha','GOT.AST.')==x))]
    if(shuffle==TRUE){
        comb_dat_funct=comb_dat_funct[sample(1:nrow(comb_dat_funct),replace = FALSE),]
    }
    library(caret);library(randomForest)
    
    holdouts=c()
    for(i in unique(types_funct2)){
      holdouts=c(holdouts,sample(which(types_funct2==i),2))
    }
    comb_dat_funct_test=comb_dat_funct[holdouts,]
    comb_dat_funct=comb_dat_funct[-holdouts,]
    types_funct2_test=types_funct2[holdouts]
    types_funct2=types_funct2[-holdouts]
    
    
    #
    
    # RF=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F)
    # hits=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
    # RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F)
    # RF$confusion
    
  
    if(length(unique(table(types_funct2)))>1){
      
      ########### add synthetic
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=3,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    if(length(unique(table(types_funct2)))>1){
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=4,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    ######################################
    #TAKE TOP xvar params into list 
    RF=randomForest(x=comb_dat_funct[,hits_xvar2], y=as.factor(types_funct2),importance=T,do.trace = F)
    hits=rownames(RF$importance)[order(RF$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar]
    
    param_hits=c(param_hits,hits)
    scores=cbind(scores,rbind(types_funct2_test,as.character(predict(RF,newdata=comb_dat_funct_test))))
  }
  
  
  # hits_final_funct=(table(param_hits)[table(param_hits)>0.90*max(table(param_hits))])[
  #   order(as.numeric((table(param_hits)[table(param_hits)>0.90*max(table(param_hits))])),decreasing=T)]
  #tabulate which paremeters came up top most often in looped random forest runs
  hits_final_funct=((table(param_hits))[
    order(as.numeric(table(param_hits)),decreasing=T)])[1:xvar]
  
  
  #############################################################################################################################################
  #############################################################################################################################################
  scores=matrix(nrow=2,ncol=0)
  
  for(test in 1:runs_of_RF){
    comb_dat_funct=comb_dat_funct_in
    #only use top parameters selected by above model
    comb_dat_funct=comb_dat_funct[,names(hits_final_funct)]
    

    types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
    types_funct2=as.character(types_funct_convert[types_funct])
    #drop highly correlated
    #pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('MIP.1.alpha','IL.1alpha','GOT.AST.')==x))]
    
    #if set to shuffle mode rotate sample order data each loop
      if(shuffle==TRUE){
          comb_dat_funct=comb_dat_funct[sample(1:nrow(comb_dat_funct),replace = FALSE),]
      }
    #test correct parameters chosen
    if(ncol(comb_dat_funct)!=xvar){
      print('error param choice not working for confusion matrix')
    }

    
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
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=4,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    if(length(unique(table(types_funct2)))>1){
      ########### add synthetic
      smote_test=smotefamily::SMOTE(comb_dat_funct,as.factor(types_funct2),K=4,dup_size=2)
      smote_test$syn_data$class
      comb_dat_funct=rbind(comb_dat_funct,smote_test$syn_data[1:7,1:(ncol(smote_test$syn_data)-1)])
      types_funct2=c(types_funct2,smote_test$syn_data$class[1:7])
    }
    
    RF=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F)
    scores=cbind(scores,rbind(types_funct2_test,as.character(predict(RF,newdata=comb_dat_funct_test))))
  }
  
  
  confusion=matrix(nrow=length(unique(types_funct2)),ncol=length(unique(types_funct2)),dat=0)
  rownames(confusion)=colnames(confusion)=score_vals=unique(unique(scores[1,]))
  for(scores_i in 1:ncol(confusion)){
    test_i=which(scores[1,]==score_vals[scores_i])
    for(scores_j in 1:ncol(confusion)){
      confusion[scores_i,scores_j]=length(which(scores[2,test_i]==score_vals[scores_j]))
    }
  }
  return(list(table(param_hits),confusion,xvar,RF_1$importance[,ncol(RF_1$importance)],
              RFE$results$Kappa,hits_final_funct))
}



###########################################
#comb_dat=cbind(mat3)
#comb_dat=comb_dat[,!grepl('rectal|Albumin',colnames(comb_dat))]
types=unlist(strsplit(rownames(comb_dat),'-'))[(1:nrow(data))*2-1]
types_convert=list()
types_convert[['SMI6']]=types_convert[['SLI6']]='2006'
types_convert[['SMI13']]=types_convert[['SLI13']]='2013'
types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8'
types_convert[['SMC']]=types_convert[['SLC']]='control'

bin1=comb_dat

funct_out_4_cats_comb=RUN_FILTRATION_AND_PREDICTION_imbalanced(bin1,types_convert,runs_of_RF,N_classes)
plot(funct_out_4_cats_comb[[1]])

write.csv(funct_out_4_cats_comb[[2]],
    file=paste0(outdir,'/four_states_performance_table_Nparams',out_file_extra,
            N_classes,'.csv'))

paste(substr(names(funct_out_4_cats_comb[[1]][order((funct_out_4_cats_comb[[1]]),decreasing = T)])[1:funct_out_4_cats_comb[[3]]],1,30),sep='\t')
print(funct_out_4_cats_comb[[3]])
print(funct_out_4_cats_comb[[2]])


#  return(list(table(param_hits),confusion,xvar,RF_1$importance[,ncol(RF_1$importance)],RFE$results$Kappa,hits_final_funct))
funct_out_4_cats_comb[[5]][1:40]
funct_out_4_cats_comb[[6]]
hits_final_4_cats_comb=names(funct_out_4_cats_comb[[1]])[
  order(funct_out_4_cats_comb[[1]],decreasing=T)[1:funct_out_4_cats_comb[[3]]]]

########################




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
rownames(comb_dat)
comb_dat=comb_dat[c(grep('SLC|SMC',rownames(comb_dat)),grep('SLC|SMC',rownames(comb_dat),invert=T)),]
clinicals_plot=clinicals_in$clinical.score[match(rownames(comb_dat),clinicals_in$ID)]



topotopparams=hits_final_4_cats_comb


write.csv(comb_dat,file=paste0(outdir,'/all_RF.data_fourstates.csv'))
write.csv(comb_dat[,topotopparams],file=paste0(outdir,'/final_four.data.'
                                              ,N_classes,'params.csv'))

#pheaty=apply(comb_dat,2,function(x) (log(x+1) - mean(log(x[grepl('SMC|SLC',rownames(comb_dat))]+1)))/sd(log(x[!is.na(x)]+1)))[,topotopparams]
#Z scores:
pheaty=apply(comb_dat,2,function(x) ((x) - mean((x[grepl('SMC|SLC',rownames(comb_dat))])))/sd((x[!is.na(x)])))[,topotopparams]


pheaty=pheaty[c(grep('SMC',rownames(pheaty)),grep('SLC',rownames(pheaty)),grep('SMI8',rownames(pheaty)),
                grep('SMI13',rownames(pheaty)),
                grep('SMI6',rownames(pheaty)),grep('SLI13',rownames(pheaty)),grep('SLI6',rownames(pheaty))),]


storage=rownames(pheaty)
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
#950 long
pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= c(colorRampPalette(c('#1842a5',rev(RColorBrewer::brewer.pal(7,'RdBu')[c(1,2,4)])))(60)))
grid.abline(107,0,range='x')
grid.abline(207,0,range='x')
#grid.abline(255.5,0,range='x')
grid.abline(307,0,range='x')
#grid.abline(307.5*1.1,0,range='x')
grid.abline(407.5,0,range='x')
grid.abline(507.5,0,range='x')
#grid.abline(531.5,0,range='x')
grid.abline(607.5,0,range='x')
#grid.abline(649.5,0,range='x')

dev.off()
png(paste0(outdir,'/plots/four_states_heatmap',
          N_classes,'params.png',sep=''),width=700,height=871.5)
pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= c(colorRampPalette(c('#1842a5',rev(RColorBrewer::brewer.pal(7,'RdBu')[c(1,2,4)])))(60)))
grid.abline(105,0,range='x')
grid.abline(205,0,range='x')
#grid.abline(255.5,0,range='x')
grid.abline(306,0,range='x')
#grid.abline(307.5*1.1,0,range='x')
grid.abline(407.5,0,range='x')
grid.abline(507.5,0,range='x')
#grid.abline(531.5,0,range='x')
grid.abline(607.5,0,range='x')
#grid.abline(649.5,0,range='x')
dev.off()

png(paste0(outdir,'/four_states_heatmap_no_lines_',
          N_classes,'params.png',),width=700,height=871.5)
pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= c(colorRampPalette(c('#1842a5',rev(RColorBrewer::brewer.pal(7,'RdBu')[c(1,2,4)])))(60)))
# grid.abline(105,0,range='x')
# grid.abline(205,0,range='x')
# grid.abline(306,0,range='x')
# grid.abline(407.5,0,range='x')
# grid.abline(507.5,0,range='x')
# grid.abline(607.5,0,range='x')

dev.off()
# 
# svg(paste('~/Pictures/plots/Sheep_megadata/14.9.21/RF_plots/four_states_heatmap',
#           N_classes,'params.svg',sep=''),width=7.00,height=12.05)
# #950 long
# pheatmap((pheaty),
#          cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
#          color= c(colorRampPalette(c('#1842a5',rev(RColorBrewer::brewer.pal(7,'RdBu')[c(1,2,4)])))(60)))
# grid.abline(105,0,range='x')
# grid.abline(205,0,range='x')
# #grid.abline(255.5,0,range='x')
# grid.abline(305,0,range='x')
# #grid.abline(307.5*1.1,0,range='x')
# grid.abline(405,0,range='x')
# grid.abline(505,0,range='x')
# #grid.abline(531.5,0,range='x')
# grid.abline(605,0,range='x')
# #grid.abline(649.5,0,range='x')
# 
# dev.off()

pdf(paste0(outdir,'/plots/four_states_heatmap',
          N_classes,'params.pdf'),width=7.00,height=10.05
    )
#950 long
pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= c(colorRampPalette(c('#1842a5',rev(RColorBrewer::brewer.pal(7,'RdBu')[c(1,2,4)])))(60)))
grid.abline(638.5-79.6*5,0,range='x')
grid.abline(638.5-79.6*4,0,range='x')
#grid.abline(255.5,0,range='x')
grid.abline(638.5-79.5*3,0,range='x')
#grid.abline(307.5*1.1,0,range='x')
grid.abline(638.5-79.4*2,0,range='x')
grid.abline(638.5-79.5,0,range='x')
#grid.abline(531.5,0,range='x')
grid.abline(638.5,0,range='x')
#grid.abline(649.5,0,range='x')

dev.off()


pdf(paste0(outdir,'/plots/four_states_heatmap_no_lines',
          N_classes,'params.pdf'),width=7.00,height=10.05
)
#950 long
pheatmap((pheaty),
         cluster_rows=F,cluster_cols=T,treeheight_row = 0, treeheight_col = 0,
         color= c(colorRampPalette(c('#1842a5',rev(RColorBrewer::brewer.pal(7,'RdBu')[c(1,2,4)])))(60)))
# grid.abline(638.5-79.6*5,0,range='x')
# grid.abline(638.5-79.6*4,0,range='x')
# #grid.abline(255.5,0,range='x')
# grid.abline(638.5-79.5*3,0,range='x')
# #grid.abline(307.5*1.1,0,range='x')
# grid.abline(638.5-79.4*2,0,range='x')
# grid.abline(638.5-79.5,0,range='x')
# #grid.abline(531.5,0,range='x')
# grid.abline(638.5,0,range='x')
#grid.abline(649.5,0,range='x')

dev.off()



