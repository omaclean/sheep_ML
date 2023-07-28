types_convert=list()
types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8_dpi',dpi_plot,sep='')
types_convert[['SMI13']]=types_convert[['SLI13']]=paste('2013_dpi',dpi_plot,sep='')
types_convert[['SMI6']]=types_convert[['SLI6']]=paste('2006_dpi',dpi_plot,sep='')
types2=as.character(types_convert[types])


##fix up BTMS- remove duplicates and melt 
for(i in 1:ncol(dat)){
  genes_i=unique(dat[2:nrow(dat),i])
  genes_i=as.character(genes_i[genes_i!=' '])
  len=length(genes_i)
  
  mat[total:(total+len-1),]=cbind(as.character(rep(names(dat)[i],len)),genes_i)
  if(any(duplicated(genes_i))){
    print(c(i,length(which(duplicated(genes_i)))))}
  total=total+len-1
}



mat=mat[!is.na(mat[,1]),]

for(i in unique(mat[,1])){
  if(any(duplicated((mat[mat[,1]==i,2])))){
    print(i)
  }
}
###########################
library(psych)
mat2=matrix(ncol=ncol(counts[2:ncol(counts)]),nrow=ncol(dat))
colnames(mat2)=colnames(counts[,2:ncol(counts)])
rownames(mat2)=names(dat)
for(i in 1:ncol(dat)){
  for(j in 2:ncol(counts)){
    mat2[i,j-1]=geometric.mean(1+counts[counts$Geneid%in%dat[2:nrow(dat),i],j])-1
  }
}
plot(mat2[,1],mat2[,31])




###########################
types=sapply(colnames(mat2),function(x) strsplit(x,'\\.')[[1]][1])
par(mfrow=c(3,2),mar=c(4,4,1,1))
for(dpi_plot in c('1','3','7')){
  
  types_convert=list()
  types_convert[['SMC']]=types_convert[['SLC']]='uninfected'
  types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8_dpi',dpi_plot,sep='')
  types_convert[['SMI13']]=types_convert[['SLI13']]=paste('2013_dpi',dpi_plot,sep='')
  types_convert[['SMI6']]=types_convert[['SLI6']]=paste('2006_dpi',dpi_plot,sep='')
  types2=as.character(types_convert[types])
  
  dpis=sapply(colnames(mat2),function(x) strsplit(x,'_')[[1]][3])
  
  types2[dpis==0]='uninfected'
  
  
  mat3=mat2[rowSums(mat2>0)>1,]
  
  pca_btm=prcomp(t(mat3),retx=T,scale.=T,center=T)
  
  to_plot=(dpis==dpi_plot)|(types2=='uninfected')
  types2_fact=factor(types2,levels=unique(unlist(types_convert)))
  
  plot(pca_btm$x[to_plot,1],pca_btm$x[to_plot,2],pch=19,
       col=(brewer.pal(length(unique(types2)),'Dark2')[as.numeric(types2_fact)])[to_plot],
       xlab=paste('PCA1 prop dev =',100*round(pca_btm$sdev[1]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA2 prop dev =',100*round(pca_btm$sdev[2]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-17,30),ylim=c(-25,20))
  legend(x='topright',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
  
  plot(pca_btm$x[to_plot,3],pca_btm$x[to_plot,4],pch=19,
       col=(brewer.pal(length(unique(types2)),'Dark2')[as.numeric(types2_fact)])[to_plot],
       xlab=paste('PCA3 prop dev =',100*round(pca_btm$sdev[3]^2/sum(pca_btm$sdev^2),4),'%'),
       ylab=paste('PCA4 prop dev =',100*round(pca_btm$sdev[4]^2/sum(pca_btm$sdev^2),4),'%'),
       main=paste('BTM PCA dpi ',dpi_plot),xlim=c(-17,30),ylim=c(-25,20))
  legend(x='topright',bty='n',legend=levels(types2_fact),pch=19,
         col=brewer.pal(length(unique(types2)),'Dark2'))
}
#################################
#################################
#################################
#################################
dpi_plot='3'
to_plot=(dpis==dpi_plot)|(types2=='uninfected')
types2_fact=factor(types2,levels=unique(unlist(types_convert)))


library(caret)
rfcont=rfeControl(functions=rfFuncs,repeats=5,verbose=F,rerank=F,method="repeatedcv",number=5)


mat2=matrix(ncol=ncol(counts[2:ncol(counts)]),nrow=ncol(dat))
colnames(mat2)=colnames(counts[,2:ncol(counts)])
for(i in 1:ncol(dat)){
  for(j in 2:ncol(counts)){
    mat2[i,j-1]=geometric.mean(1+counts[counts$Geneid%in%dat[2:nrow(dat),i],j])-1
  }
}
rownames(mat2)=names(dat)


##########################
mat3=mat2[rowSums(mat2>0)>1,]
rfedat=t(mat3[,to_plot])

RFE1=rfe(rfedat, y=types2_fact[to_plot],
    rfeControl = rfcont,sizes=(1:nrow(mat3)))

test_set=RFE1$optVariables[1:17]
library(randomForest)


randomForest(y= types2_fact[to_plot],x= rfedat[,colnames(rfedat)%in%test_set],
             importance=T,do.trace = F)

par(mfrow=c(1,2))
plot(RFE1$results$Variables,RFE1$results$Accuracy,ylim=c(0.5,1),ylab='',xlab='')
par(new=T)
plot(RFE1$results$Variables,RFE1$results$Kappa,ylim=c(0.5,1),col=2,xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=30,col=3)
abline(v=RFE1$optsize,col=4)
legend(x='bottomright',legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),pch=1,col=1:2,bty='n')


plot(RFE1$results$Variables,RFE1$results$Accuracy,ylim=c(0.5,1),xlim=c(0,80),type='l',ylab='',xlab='')
par(new=T)
plot(RFE1$results$Variables,RFE1$results$Kappa,ylim=c(0.5,1),col=2,xlim=c(0,80),type='l',xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=17,col=3)
abline(v=RFE1$optsize,col=4)
legend(x='right',legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),lty=1,col=1:2,bty='n')
################################################################
################################################################
################################################################
################################################################
################################################################################################

animals_due_to_die=rep(1,length(which(to_plot)))
animals_due_to_die[grepl('SLI6',rownames(rfedat))&dpis[to_plot]=='3']=2
death=factor(c('yes','no')[animals_due_to_die],levels=c('yes','no'))

RFE2=rfe(rfedat, y=death,
         rfeControl = rfcont,sizes=(1:nrow(mat3)))
RFE2

par(mfrow=c(1,2))
plot(RFE2$results$Variables,RFE2$results$Accuracy,ylim=c(0.5,1),ylab='',xlab='')
par(new=T)
plot(RFE2$results$Variables,RFE2$results$Kappa,ylim=c(0.5,1),col=2,xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=5,col=3)

legend(x='right',legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),pch=1,col=1:2,bty='n')

plot(RFE2$results$Variables,RFE2$results$Accuracy,ylim=c(0.5,1),xlim=c(0,20),type='l',ylab='',xlab='')
par(new=T)
plot(RFE2$results$Variables,RFE2$results$Kappa,ylim=c(0.5,1),col=2,xlim=c(0,20),type='l',xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=5,col=3)
legend(x='right',legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),lty=1,col=1:2,bty='n')

RF=randomForest(y= death,x= rfedat[,colnames(rfedat)%in%RFE2$optVariables],
             importance=T,do.trace = F)
RF$predicted[RF$predicted=='no']
RF$predicted[RF$predicted=='yes'&death=='no']


################################################################
################################################################
################################################################
################################################################################################
dpi_plot='3'
cond='SLI6'

to_plot=(dpis==dpi_plot)|(types2=='uninfected')
animal_states=rep(1,length(which(to_plot)))
rfedat=t(mat3[,to_plot])
animal_states[grepl(cond,rownames(rfedat))&dpis[to_plot]=='1']=2
BTV8DPI1=factor(c('no','yes')[animal_states],levels=c('no','yes'))

RFE3=rfe(rfedat, y=BTV8DPI1,
         rfeControl = rfcont,sizes=(1:nrow(mat3)))
RFE3

par(mfrow=c(1,2))
plot(RFE3$results$Variables,RFE3$results$Accuracy,ylim=c(0,1),ylab='',xlab='')
par(new=T)
plot(RFE3$results$Variables,RFE3$results$Kappa,ylim=c(0,1),col=2,xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=RFE3$bestSubset,col=3)

legend(x=20,y=0.8,legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),pch=1,col=1:2,bty='n')

plot(RFE3$results$Variables,RFE3$results$Accuracy,ylim=c(0.5,1),xlim=c(0,20),type='l',ylab='',xlab='')
par(new=T)
plot(RFE3$results$Variables,RFE3$results$Kappa,ylim=c(0.5,1),col=2,xlim=c(0,20),type='l',xlab='number of parameters in run',
     ylab='Accuracy metric',main='random forest performance')
abline(v=RFE3$bestSubset,col=3)
abline(v=11,col=4)
legend(x=c(0.4,0.9),legend=c('accuracy','kappa (biased class \n adjusted accuracy)'),lty=1,col=1:2,bty='n')
RF=randomForest(y= BTV8DPI1,x= rfedat[,colnames(rfedat)%in%RFE3$optVariables[1:11]],
                importance=T,do.trace = F)
RF
RFE3$optVariables[1:11]








# #### [1] 
# dpi7_params=c("Activated..LPS..dendritic.cell.surface.signature..S11."       
# , "type.I.interferon.response..M127."                            
# ,"proinflammatory.dendritic.cell..myeloid.cell.response..M86.1."
# ,"DC.surface.signature..S5."                                    
# ,"RIG.1.like.receptor.signaling..M68."                          
# ,"TBA..M131."                                                   
# , "viral.sensing...immunity..IRF2.targets.network..II...M111.1." 
# ,"signaling.in.T.cells..II...M35.1."                            
# ,"innate.antiviral.response..M150."                             
# , "antiviral.IFN.signature..M75."                                
# , "signaling.in.T.cells..I...M35.0."                             
# ,"myeloid..dendritic.cell.activation.via.NFkB..II...M43.1."     
# ,"TBA..M66."                                                    
# ,"inflammasome.receptors.and.signaling..M53."                   
# , "proteasome..M226."                                            
# , "innate.activation.by.cytosolic.D..sensing..M13."              
# ,"enriched.in.activated.dendritic.cells..II...M165."            
# ,"enriched.in.NK.cells..I...M7.2."                              
# , "golgi.membrane..I...M113."                                    
# , "T.cell.activation..III...M7.4."                               
# , "activated.dendritic.cells..M67."                              
# , "CORO1A.DEF6.network..I...M32.2."                              
# , "enriched.in.nuclear.pore.complex.interacting.proteins..M247." 
# ,"T.cell.differentiation.via.ITK.and.PKC..M18."                 
# ,"collagen..TGFB.family.et.al..M77."                            
# ,"chaperonin.mediated.protein.folding..I...M204.0."             
# , "chaperonin.mediated.protein.folding..II...M204.1."            
# , "T.cell.differentiation..Th2...M19."                           
# , "mismatch.repair..II...M22.1."                                 
# , "TBA..M121."  )

