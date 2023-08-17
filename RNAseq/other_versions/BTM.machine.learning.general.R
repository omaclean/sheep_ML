library(psych);library(ggplot2);library(caret);library(randomForest)
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
for(i in 2:ncol(counts)){
  counts[,i]=counts[,i]/sum(counts[,i])
}
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

## get geometric means out for counts
mat2=matrix(ncol=ncol(counts[2:ncol(counts)]),nrow=ncol(dat))
colnames(mat2)=colnames(counts[,2:ncol(counts)])
rownames(mat2)=names(dat)
for(i in 1:ncol(dat)){
  for(j in 2:ncol(counts)){
    mat2[i,j-1]=geometric.mean(1+counts[counts$Geneid%in%dat[2:nrow(dat),i],j])-1
  }
}

# convert colnames from counts into 'types' &extract DPI
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
types_convert=list()
types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8_dpi',dpi_plot,sep='')
types_convert[['SMI13']]=types_convert[['SLI13']]=paste('2013_dpi',dpi_plot,sep='')
types_convert[['SMI6']]=types_convert[['SLI6']]=paste('2006_dpi',dpi_plot,sep='')
types2=as.character(types_convert[types])
dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])


############################################################
################################################################
################################################################################################
dpi_plot='7'
cond='SLI13'

to_plot=(dpis==dpi_plot)|(types2=='uninfected')
animal_states=rep(1,length(which(to_plot)))
rfedat=t(mat3[,to_plot])
animal_states[grepl(cond,rownames(rfedat))&dpis[to_plot]==dpi_plot]=2
BTV8DPI1=factor(c('no','yes')[animal_states],levels=c('no','yes'))

RFE3=rfe(rfedat, y=BTV8DPI1,
         rfeControl = rfcont,sizes=(1:(nrow(mat3))))

#plot optimum BTM set
par(mfrow=c(1,2))
plot(RFE3$results$Variables,RFE3$results$Accuracy,ylim=c(0,1),ylab='',xlab='')
par(new=T)
plot(RFE3$results$Variables,RFE3$results$Kappa,ylim=c(0,1),col=2,xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=RFE3$bestSubset,col=3)
legend(x=30,y=0.8,legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),pch=1,col=1:2,bty='n')
############################################# zoomed in plot
plot(RFE3$results$Variables,RFE3$results$Accuracy,ylim=c(0.2,1),xlim=c(0,20),type='l',ylab='',xlab='')
par(new=T)
plot(RFE3$results$Variables,RFE3$results$Kappa,ylim=c(0.2,1),col=2,xlim=c(0,20),type='l',xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=RFE3$bestSubset,col=3)
abline(v=9,col=4)
legend(x=c(0.4,0.9),legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),lty=1,col=1:2,bty='n')
############ get confusion matrix
RF=randomForest(y= BTV8DPI1,x= rfedat[,colnames(rfedat)%in%RFE3$optVariables[1:11]],
                importance=T,do.trace = F)
RF
RFE3$optVariables[1:9]


#########################
# plot
##

BTMs_to_plot=RFE3$bestSubset

BTMs_to_plot=RFE3$optVariables[1:9]
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
png(filename=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn/'
                   ,cond,'_dpi',dpi_plot,'.png',sep=''),width=1100,height=800)

grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]]
             #,plots[[5]], plots[[6]],plots[[7]],plots[[8]],plots[[9]]
             ,ncol=1)
dev.off()




#########################
# plot
# ##
# 
# BTMs_to_plot=grep('T.cell',unique(mat[,1]),value=T)
# 
# mat_plot=mat2[match(BTMs_to_plot,rownames(mat2)),]
# 
# dpis=sapply(colnames(mat_plot),function(x) strsplit(x,'_')[[1]][3])
# treatments=sapply(colnames(mat_plot),function(x) strsplit(x,'\\.')[[1]][1])
# colnames(mat_plot)=paste(treatments,dpis,sep='_')
# 
# cond_list_order=c('SMC_0','SMC_1','SMC_3','SMC_7',
#                   'SLC_0','SLC_1','SLC_3','SLC_7',
#                   'SMI8_0','SMI13_0','SMI6_0',
#                   'SLI13_0','SLI6_0',
#                   'SMI8_1','SMI13_1','SMI6_1',
#                   'SLI13_1','SLI6_1',
#                   'SMI8_3','SMI13_3','SMI6_3',
#                   'SLI13_3','SLI6_3',
#                   'SMI8_7','SMI13_7','SMI6_7',
#                   'SLI13_7','SLI6_7')
# types_convert_plot=list()
# types_convert_plot[['SMC']]=types_convert_plot[['SLC']]='control'
# types_convert_plot[['SLI8']]=types_convert_plot[['SMI8']]=paste('BTV8',sep='')
# types_convert_plot[['SMI13']]=types_convert_plot[['SLI13']]=paste('2013',sep='')
# types_convert_plot[['SMI6']]=types_convert_plot[['SLI6']]=paste('2006',sep='')
# 
# 
# plots=list()
# for(i in 1:nrow(mat_plot)){
#   BTM_to_plot=i
#   to_plot=data.frame(vals=mat_plot[BTM_to_plot,],conds=colnames(mat_plot),treatment=treatments)
#   to_plot$conds=factor(to_plot$conds,levels = cond_list_order)
#   to_plot$treatment_col=factor(types_convert_plot[to_plot$treatment],levels=c('control','BTV8','2013','2006'))
#   plots[[i]]=ggplot(to_plot,aes(x=conds,y=vals,fill=treatment_col))+
#     geom_dotplot(binaxis='y',stackdir = 'center',alpha=0.3,dotsize = 1.8)+theme_bw()+
#     ggtitle(rownames(mat_plot)[i])+
#     geom_vline(xintercept=(13.5))+
#     geom_vline(xintercept=(18.5))+
#     geom_vline(xintercept=(23.5))
# }
# library(gridExtra)
# png(filename=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn/Tcells.allconds'
#                    ,'.png',sep=''),width=1100,height=3000)
# 
# grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],
#              plots[[6]],plots[[7]],plots[[8]],plots[[9]], plots[[10]],plots[[11]],plots[[12]],plots[[13]],
#              plots[[14]],plots[[15]],plots[[16]],plots[[17]],plots[[18]],plots[[19]],plots[[20]],
#              plots[[21]],plots[[22]]
#              ,ncol=1)
# dev.off()
# 
# 
# 
# 
