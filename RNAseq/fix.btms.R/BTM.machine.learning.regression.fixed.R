library(psych);library(ggplot2);library(caret);library(randomForest)
dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
clinical_score=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv')
####################
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]
#remove BTMs which are contained in another (i.e. full subsets)
redundant_BTMS=c("Rho.GTPase.cycle..M4.14.","mitotic.cell.division..M6.",
                 "enriched.in.NK.cells..KIR.cluster...M61.1.","enriched.in.B.cells..VI...M69.")

dat=dat[,!colnames(dat)%in%redundant_BTMS]

for(i in 2:ncol(counts)){#adjust counts for #reads/run
  counts[,i]=counts[,i]/sum(counts[,i])
}
for(i in 1:nrow(counts)){#adjust counts for #reads/gene (increases variance but avoids overweighting highly expressed genes- assumes genes vary in potency)
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

dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
animal_names=gsub('\\.','-',sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][1]))
clinical_score_i=sapply(1:ncol(mat2),function(x) clinical_score$clinical.score[
        clinical_score$dpi==as.numeric(dpis[x])&clinical_score$ID==animal_names[x]])
############################################################
################################################################
################################################################################################
rfedat=t(mat2)
library(caret)
rfcont=rfeControl(functions=rfFuncs,repeats=5,verbose=F,rerank=F,method="repeatedcv",number=5)
# RFE_reg_save=RFE_reg
# 
# ########################################
# RFE_reg=rfe(rfedat[,RFE_reg_save$variables$var[1:50]], y=clinical_score_i,
#          rfeControl = rfcont,sizes=(1:(nrow(mat3))))
#######################
#plot optimum BTM set
par(mfrow=c(1,2))
plot(RFE_reg$results$Variables,RFE_reg$results$Rsquared,ylim=c(0,1),ylab='',xlab='')
par(new=T)
plot(RFE_reg$results$Variables,RFE_reg$results$MAE,ylim=c(0,1),col=2,xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=RFE_reg$bestSubset,col=3)
legend(x='bottomright',legend=c('Rsquared','Mean absolute error'),pch=1,col=1:2,bty='n')
############################################# zoomed in plot
plot(RFE_reg$results$Variables,RFE_reg$results$Rsquared,ylim=c(0.4,1),xlim=c(0,30),type='l',ylab='',xlab='')
par(new=T)
plot(RFE_reg$results$Variables,RFE_reg$results$MAE,ylim=c(0.4,1),col=2,xlim=c(0,30),type='l',xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=RFE_reg$bestSubset,col=3)
abline(v=13,col=4)
legend(x='bottomright',legend=c('Rsquared','Mean absolute error'),lty=1,col=1:2,bty='n')
######## get confusion matrix

RF=randomForest(rfedat[,RFE_reg$optVariables], y=clinical_score_i,
                importance=T,do.trace = F)
par(mfrow=c(1,1))
col_names=c('uninfected','dpi1','dpi3','dpi7')
cols=rep(1,length(dpis))
cols[dpis=='1']=2;cols[dpis=='3']=3;cols[dpis=='7']=4;cols[grepl('SMC|SLC',animal_names)]=1
plot(clinical_score_i,RF$predicted,pch=19,col=alpha(brewer.pal(4,'Dark2')[cols],.5))
legend(x='topleft',bty='n',legend=col_names,col=brewer.pal(4,'Dark2'),pch=19)
RFE_reg$optVariables[1:13]
library(plotly)

plot_ly(x=clinical_score_i,y=RF$predicted,color=factor(col_names[cols],levels=col_names),colors=brewer.pal(4,'Dark2'),alpha=0.85,
        text=~paste(sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][1]),'dpi=',dpis))%>%
  add_markers(marker=list(size=12))%>%
  layout(xaxis=list(title='actual clinical score'),
         yaxis=list(title='predicted clinical score'),
         legend=list(font=list(size=20)))



#########################
# plot
##

BTMs_to_plot=RFE_reg$bestSubset

BTMs_to_plot=RFE_reg$optVariables[1:13]
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
png(filename=paste('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn.new/regression_clinical_score2.png',sep=''),width=1100,height=200*length(plots))

grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]]
             ,plots[[5]], plots[[6]],plots[[7]],plots[[8]],plots[[9]]
             ,plots[[10]],plots[[11]],plots[[12]],plots[[13]]
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



