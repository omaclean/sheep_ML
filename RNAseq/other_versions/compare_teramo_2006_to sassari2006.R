library(edgeR)
genes_hits=list()
dpi_levels=c('0'#,'1','3'
             ,'3')
location1='Sassari'
types_list1=c("BTV8")
types_headers1=c('SMI8')

location2='Sassari'
types_list2=c("2006")
types_headers2=c('SMI6')


###########################


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
###########
types_list=types_list1
types_headers=types_headers1

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]


dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))


counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]

dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

dpis=factor(dpis,levels=dpi_levels)
animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))

animal_names2=as.factor(animal_names)

################################################
design <- model.matrix(~dpis+animal_names2)
##################
d3=DGEList(counts=counts2)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=paste(animal_names,'dpi',as.character(dpis))
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results=matrix(ncol=3,nrow=2)
results[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )


contrast=rep(0,ncol(design))
contrast[2]=1
de2=glmQLFTest(d3,contrast=contrast)
test= topTags(de2, n=nrow(counts2))
FDR_threshold=0.05


genes_hits[[1]]=rownames(test$table)
genes_hits[[2]]=test$table$FDR
genes_hits[[3]]=test$table$logFC

#2)
###########
types_list=types_list2

types_headers=types_headers2

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]


dpis=sapply(colnames(counts2),function(x) 
  tail(strsplit(strsplit(x,'_dpi')[[1]][1],'')[[1]],1))


counts2=counts2[,dpis==dpi_levels[1]|dpis==dpi_levels[2]]

dpis=dpis[dpis==dpi_levels[1]|dpis==dpi_levels[2]]

dpis=factor(dpis,levels=dpi_levels)
animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))

animal_names2=as.factor(animal_names)

################################################
design <- model.matrix(~dpis+animal_names2)
##################
d3=DGEList(counts=counts2)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=paste(animal_names,'dpi',as.character(dpis))
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results=matrix(ncol=3,nrow=2)
results[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )


contrast=rep(0,ncol(design))
contrast[2]=1
de2=glmQLFTest(d3,contrast=contrast)
test= topTags(de2, n=nrow(counts2))
FDR_threshold=0.05


genes_hits[[4]]=rownames(test$table)
genes_hits[[5]]=test$table$FDR
genes_hits[[6]]=test$table$logFC
results[2,]=c(length(which(test$table$FDR<FDR_threshold)),colnames(design)[1],colnames(design)[2])
###########################################################################################################
################################################################################################
################################################################################################
################################################
max_FDR=ceiling(max(c(-log(genes_hits[[2]]),-log(genes_hits[[5]])))/50)*50
max_logFC=ceiling(max(abs(c(genes_hits[[3]],genes_hits[[6]]))))


names_plot=gene_names$gene_name[match(genes_hits[[1]],gene_names$sheep_ID)]
names_plot[(log(genes_hits[[2]])>(-max_FDR*0.5))&
             (log(genes_hits[[2]])>(-max_FDR*0.2)| abs(genes_hits[[3]])<2.5)]=''

par(mfrow=c(2,2),mar=c(3.8, 4, 2, 0) )
log(genes_hits[[2]])

c(range(-log(genes_hits[[2]])),range(-log(genes_hits[[5]])))
c(range((genes_hits[[3]])),range((genes_hits[[6]])))

plot((genes_hits[[3]]),-log(genes_hits[[2]]),
     pch=19,col=alpha('#0b895b',0.25),
     xlab=paste('logFC',types_list1,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1],location1),
     ylab=paste('log(FDR)',types_list1,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1],location1),
     main='model includes individual variation',yaxt='n',
     ylim=c(0,max_FDR),xlim=c(-max_logFC,max_logFC),
     col=alpha('#000000',.4))
axis(2,at=max_FDR*(0:5)/5,labels=paste('1e-',max_FDR*(0:5)/5,sep=''))
abline(h=c(-log(0.05)),v=c(-2.5,2.5),col=alpha('#2bA97b',0.5))

text((genes_hits[[3]]),-log(genes_hits[[2]]),col=alpha('#000000',0.5),
     names_plot,srt=-45)


############################################################

names_plot=gene_names$gene_name[match(genes_hits[[4]],gene_names$sheep_ID)]

names_plot[(log(genes_hits[[5]])>(-max_FDR*0.5))&
             (log(genes_hits[[5]])>(-max_FDR*0.2)| abs(genes_hits[[6]])<2.5)]=''


plot((genes_hits[[6]]),-log(genes_hits[[5]]),
     pch=19,col=alpha('#0b895b',0.25),
     xlab=paste('logFC',types_list1,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1],location2),
     ylab=paste('log(FDR)',types_list1,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1],location2),
     main='model includes individual variation',yaxt='n',
     ylim=c(0,max_FDR),xlim=c(-max_logFC,max_logFC))
axis(2,at=max_FDR*(0:5)/5,labels=paste('1e-',max_FDR*(0:5)/5,sep=''))
abline(h=c(-log(0.05)),v=c(-2.5,2.5),col=alpha('#2bA97b',0.5))

text((genes_hits[[6]]),-log(genes_hits[[5]]),col=alpha('#000000',0.5),
     names_plot,srt=-45)

############################################################
print(results)

names_plot=gene_names$gene_name[match(genes_hits[[1]],gene_names$sheep_ID)]
names_plot[log(genes_hits[[2]])>(-max_FDR*0.3)&
             log(genes_hits[[5]])[match(genes_hits[[1]],genes_hits[[4]])]>(-max_FDR*0.3) ]=''

log(genes_hits[[2]])

c(range(-log(genes_hits[[2]])),range(-log(genes_hits[[5]])))
plot(-log(genes_hits[[2]]),-log(genes_hits[[5]])[match(genes_hits[[1]],genes_hits[[4]])],
     pch=19,col=alpha('#0b895b',0.25),
     xlab=paste('log(FDR)',types_list1,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1] ,location1),
     ylab=paste('log(FDR)',types_list2,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1],location2),
     main='model includes individual variation',yaxt='n',xaxt='n',
     ylim=c(0,max_FDR),xlim=c(0,max_FDR))
axis(2,at=max_FDR*(0:5)/5,labels=paste('1e-',max_FDR*(0:5)/5,sep=''))
axis(1,at=max_FDR*(0:5)/5,labels=paste('1e-',max_FDR*(0:5)/5,sep=''))

text(-log(genes_hits[[2]]),
     -log(genes_hits[[5]])[match(genes_hits[[1]],genes_hits[[4]])],
     names_plot,srt=-45)
abline(h=-log(0.05),v=-log(0.05),col=alpha('#2bA97b',0.8))


############


names_plot1=names_plot=gene_names$gene_name[match(genes_hits[[1]],gene_names$sheep_ID)]
names_plot1[abs(genes_hits[[3]])<2.5|
            abs(genes_hits[[6]][match(genes_hits[[1]],genes_hits[[4]])])<2.5]=''

c(range(genes_hits[[3]]),range(genes_hits[[6]]))
plot(genes_hits[[3]],genes_hits[[6]][match(genes_hits[[1]],genes_hits[[4]])],
     pch=19,col=alpha('#0b895b',0.25),
 
     xlab=paste('logFC',types_list1,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1] ,location1),
     ylab=paste('logFC',types_list2,'dpi',dpi_levels[2],'vs dpi',
                dpi_levels[1],location2),
     main='model includes individual variation',
     ylim=c(-max_logFC,max_logFC),xlim=c(-max_logFC,max_logFC))

text(genes_hits[[3]],
     genes_hits[[6]][match(genes_hits[[1]],genes_hits[[4]])],
     names_plot1,srt=-45,col=alpha('#000000',.4))
abline(h=c(-2.5,2.5),v=c(-2.5,2.5),col=alpha('#2bA97b',0.8))


# 
# cbind(genes_hits[[3]],genes_hits[[6]][match(genes_hits[[1]],genes_hits[[4]])])[which(names_plot1=='ALAS2'),]
# cbind(genes_hits[[2]],genes_hits[[5]][match(genes_hits[[1]],genes_hits[[4]])])[which(names_plot1=='ALAS2'),]



#######################################


########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
# ########################################################
# ########################################################
# plot(test$table$logFC,-log(test$table$FDR),pch=19,col=alpha('#0b895b',0.25),
#      yaxt='n',ylim=c(0,-log(10^-125)),ylab='FDR',xlab='logFC',
#      main=paste('baseline vs ',colnames(design)[3],'with individual effect terms'))
# axis(2,at=-log(10^(-10*c(2*(0:6)))),labels=10^(-10*c(2*(0:6))))
# 
# 
# 
# dim(counts2)
# plot(genes_hits[[3]],genes_hits[[6]][match(genes_hits[[1]],genes_hits[[4]])]
#      ,pch=19,col=alpha('#0b895b',0.25),
#      xlab='logFC 2006 dpi 7 vs dpi 0 Sassari',
#      ylab='logFC 2006 dpi 7 vs dpi 0 Teramo',
#      main='model includes individual variation'
# )
# 
# head(cbind(genes_hits[[1]],genes_hits[[4]][match(genes_hits[[1]],genes_hits[[4]])]))
# 
# 
# 
# 
# 
# 
# 
# test$table$FDR[test$table$FDR<10^-125]=10^-125
# names_plot=gene_names$gene_name[match(rownames(test$table),gene_names$sheep_ID)]
# names_plot[(abs(test$table$logFC)<3|is.na(names_plot)|test$table$FDR>0.01)&test$table$FDR>10^-30]=''
# plot(test$table$logFC,-log(test$table$FDR),pch=19,col=alpha('#0b895b',0.25),
#      yaxt='n',ylim=c(0,-log(10^-125)),ylab='FDR',xlab='logFC',
#      main=paste('baseline vs ',colnames(design)[3],'with individual effect terms'))
# axis(2,at=-log(10^(-10*c(2*(0:6)))),labels=10^(-10*c(2*(0:6))))
# 
# text(test$table$logFC,-log(test$table$FDR),names_plot,srt=-45)
# abline(h=-log(0.05),v=c(-3,3),col=alpha('#2bA97b',0.8))

