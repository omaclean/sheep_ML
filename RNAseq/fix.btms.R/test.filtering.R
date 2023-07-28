library(psych);library(ggplot2);library(caret);library(randomForest)

if(!exists('run')){
  dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
  counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
  counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
  ####################
  #remove BTMs which are contained in another (i.e. full subsets)
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
  
  
  
  
  
  
  #################
  mat2_copy=mat2
  
  
  mat2_cor=cor(t(mat2))
  cors=rep(0,nrow(mat2_cor))
  
  for(i in 1:(ncol(mat2_cor)-1)){
    count=0
    for(j in (i+1):ncol(mat2_cor)){
      if(mat2_cor[i,j]>0.95){
        cors[i]=cors[i]+1
        rownames(mat2)[j]=paste('mean(',rownames(mat2)[i],rownames(mat2)[j],')')
        mat2[j,]=sapply(1:ncol(mat2),function(x) mean(mat2[j,x],mat2[i,x]))
      }
    }
  }
  table(cors)
  mat2=mat2[cors==0,]
  
  mat3=mat2
  
  ####
  # convert colnames from counts into 'types' &extract DPI
  types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
  types_convert=list()
  types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
  types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8_dpi'
  types_convert[['SMI13']]=types_convert[['SLI13']]='2013_dpi'
  types_convert[['SMI6']]=types_convert[['SLI6']]='2006_dpi'
  types2=as.character(types_convert[types])
  dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
  types3=paste(types2,'_dpi',dpis)
  
}

run=1
############################################################
###############################################################################################
################################################################################################
############################################################
###############################################################################################
################################################################################################
############################################################
###############################################################################################
################################################################################################
############################################################
###############################################################################################
################################################################################################
############################################################
###############################################################################################
################################################################################################
dpi_plot='1'
cond='SMI8'#'SLI6|SMI6|SLI13|SMI13'

to_plot=(dpis==dpi_plot)|(types2=='uninfected'|dpis=='0')
animal_states=rep(1,length(which(to_plot)))
rfedat=t(mat3[,to_plot])
animal_states[grepl(cond,rownames(rfedat))&dpis[to_plot]==dpi_plot]=2
BTV8DPI1=factor(c('no','yes')[animal_states],levels=c('no','yes'))

#######################################################################
infected=!(types2=='uninfected'|dpis=='0')[to_plot]
excess_var=rep(0,ncol(rfedat))
for(i in 1:ncol(rfedat)){
  if(var(rfedat[!infected,i])>var(rfedat[,i])){
    excess_var[i]=1
  }
}
sum(excess_var)/length(excess_var)
rfedat=rfedat[,as.logical(excess_var)]


library(caret)
rfcont=rfeControl(functions=rfFuncs,repeats=5,verbose=F,rerank=F,method="repeatedcv",number=5)


cbind(BTV8DPI1,rownames(rfedat))




RF=randomForest(y= BTV8DPI1,x= rfedat,
                importance=T,do.trace = F)
RF=randomForest(y= BTV8DPI1,x= rfedat[,colnames(rfedat)%in%order(RF$importance[,4],decreasing=T)[1:50]],
                importance=T,do.trace = F)
RF$importance[1:10,]
########################################
# 
# RFE3=rfe(rfedat[,colnames(rfedat)%in%rownames(RF$importance)[order(RF$importance[,4],decreasing=T)[1:50]]], y=BTV8DPI1,
#          rfeControl = rfcont,sizes=(1:50))#(nrow(mat3))))
# #######################
# #plot optimum BTM set
# par(mfrow=c(1,2))
# plot(RFE3$results$Variables,RFE3$results$Accuracy,ylim=c(0,1),ylab='',xlab='')
# par(new=T)
# plot(RFE3$results$Variables,RFE3$results$Kappa,ylim=c(0,1),col=2,xlab='number of parameters in run',
#      ylab='Accuracy metric',main='random forest performance')
# abline(v=RFE3$bestSubset,col=3)
# legend(x='bottomright',legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),pch=1,col=1:2,bty='n')
# ############################################# zoomed in plot
# plot(RFE3$results$Variables,RFE3$results$Accuracy,ylim=c(0.2,1),xlim=c(0,20),type='l',ylab='',xlab='')
# par(new=T)
# plot(RFE3$results$Variables,RFE3$results$Kappa,ylim=c(0.2,1),col=2,xlim=c(0,20),type='l',xlab='number of parameters in run',
#      ylab='Accuracy metric',main='random forest performance')
# abline(v=RFE3$bestSubset,col=3)
# abline(v=4,col=4)
# legend(x='bottomright',legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),lty=1,col=1:2,bty='n')
# ############ get confusion matrix
# RF=randomForest(y= BTV8DPI1,x= rfedat[,colnames(rfedat)%in%RFE3$optVariables[1:4]],
#                 importance=T,do.trace = F)
# RF

#########################
# plot
##

BTMs_to_plot=rownames(RF$importance)[order(RF$importance[,4],decreasing=T)[1:11]]

#BTMs_to_plot=RFE3$optVariables[1:11]
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
png(filename=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn.new/'
                   ,cond,'_dpi',dpi_plot,'.var.antifilt.png',sep=''),width=1100,height=200*length(plots))

grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]]
             ,plots[[5]], plots[[6]],plots[[7]],plots[[8]],plots[[9]]
             ,plots[[10]],plots[[11]]#,plots[[12]]
             ,ncol=1)
dev.off()
print(length(BTMs_to_plot))



