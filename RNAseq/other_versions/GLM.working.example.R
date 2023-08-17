##############
gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

dpi=7
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
rownames(counts)=counts$Geneid
counts=counts[,c(1,grep(paste(dpi,'_dpi|0_dpi',sep=''),colnames(counts)))]

counts2=counts[,2:ncol(counts)]
counts2=counts2[,grepl(paste(types_headers,collapse='|',sep=''),colnames(counts2))]
counts2=counts2[rowSums(counts2>0)>7,]
###########
types_list=c("control","2013","2006","BTV8")
types_headers=c('SLC|SMC','SLI13|SMI13','SLI6|SMI6','SMI8')
#types_array=
animal_names=as.character(sapply(colnames(counts2),function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
inf_status=as.character(sapply(colnames(counts2),function(x) 
  types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))
inf_status[grep('0_dpi',colnames(counts2))]="control"

inf_status=factor(inf_status,levels=types_list)
locations=c("short","long")

location=factor(locations[sapply(colnames(counts2),function(x) 
  grep(substr(strsplit(x,'\\.')[[1]][1],1,2) ,c('SM','SL') ))],levels=locations)


animal_names2=as.factor(animal_names)
################################################
design <- model.matrix(~inf_status+animal_names2)
##################
d3=DGEList(counts=counts2,group= inf_status)
d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)
rownames(design)=animal_names
library(edgeR)
print(paste('number DEGs , condition 1, condition 2 on dpi',dpi))

results=matrix(ncol=3,nrow=choose(length(unique(types_list)),2)+2)
results[1,]=c(paste('number DEGs on dpi',dpi), 'condition 1', 'condition 2' )

cont=rep(0,length(design[1,]))
cont[3]=1
de2=glmQLFTest(d3,contrast=cont)
test= topTags(de2, n=nrow(counts2))
test2=test
##############
test$table$FDR[test$table$FDR<10^-125]=10^-125
names_plot=gene_names$gene_name[match(rownames(test$table),gene_names$sheep_ID)]
names_plot[(abs(test$table$logFC)<3|is.na(names_plot)|test$table$FDR>0.01)&test$table$FDR>10^-30]=''
plot(test$table$logFC,-log(test$table$FDR),pch=19,col=alpha('#0b895b',0.25),
     yaxt='n',ylim=c(0,-log(10^-125)),ylab='FDR',xlab='logFC',
     main=paste('baseline vs ',colnames(design)[3],'with individual effect terms'))
axis(2,at=-log(10^(-10*c(2*(0:6)))),labels=10^(-10*c(2*(0:6))))

text(test$table$logFC,-log(test$table$FDR),names_plot,srt=-45)
abline(h=-log(0.05),v=c(-3,3),col=alpha('#2bA97b',0.8))

###############################
results[2,]=c(length(which(test$table$FDR<0.05)),'baseline',paste(colnames(design)[cont]))
count=3
for(i in 1:(length(types_list)-1)){
  for(j in (i+1):length(types_list)){
    contrast=as.numeric(types_list%in%types_list[j])
    if(types_list[i]!='control'){
      contrast[which(types_list%in%types_list[i])]=-1
    }
    de2=glmQLFTest(d3,contrast=c(contrast,rep(0,length(design[1,])-length(contrast))))
    test= topTags(de2, n=nrow(counts))
    results[count,]=(c(length(which(test$table$FDR<0.05)),colnames(design)[i],colnames(design)[j]))
    count=count+1
  }
}
print(results)

