
counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
sex=read.csv('~/Documents/sheep_megadata/temperature_humidity_gender.csv')

dpi=7
types_headers=c('SMC','SMI6')


counts=counts[,c(1,grep(paste(types_headers,collapse ='|'),colnames(counts)))]
dpis=sapply(colnames(counts[,2:ncol(counts)]),function(x) strsplit(x,'_')[[1]][3])
counts=counts[,c(1,1+which(dpis=='0'|dpis=='7'))]
dpis=dpis[which(dpis=='0'|dpis=='7')]


types_list=c("male","infected","infected:male")

#types_array=
animal_names=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) 
  paste(strsplit(x,'\\.|_')[[1]][1:2],sep='',collapse='-')))
animal_sex=sex$gender.1[match(animal_names,sex$animal.number)]



inf_status=types_list[animal_sex+2*as.numeric(grepl('SMI6',animal_names)&dpis=='7')]
cbind(inf_status,animal_names,dpis)
# inf_status=as.character(sapply(colnames(counts)[2:ncol(counts)],function(x) 
#   types_list[grep(strsplit(x,'\\.')[[1]][1],types_headers)]))
library(edgeR)
inf_status=as.numeric(grepl('SMI6',animal_names)&dpis=='7')
male=animal_sex-1
design <- model.matrix(~male+inf_status+male:inf_status)

counts2=counts[,2:ncol(counts)]
rownames(counts2)=counts[,1]
counts2=counts2[rowSums(counts2>0)>7,]

d3=DGEList(counts=counts2)


d3=calcNormFactors(d3)
d3 <- estimateDisp(d3,design)
d3 <- glmQLFit(d3,design)

head(design)
contrast=c(0,0,0,0)
contrast[4]=1

de2 =glmQLFTest(d3, contrast =contrast  )
test= topTags(de2, n=nrow(counts))
print(c(length(which(test$table$FDR<0.05))))


print(test$table[test$table$FDR<0.05,])

orths=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')

orths$gene_names_new[orths$ens_IDs%in%rownames(test$table[test$table$FDR<0.05,])]
tab=cbind(orths$gene_names_new[orths$ens_IDs%in%rownames(test$table[test$table$FDR<0.05,])],
      test$table[test$table$FDR<0.05,c(1,5)])

colnames(tab)[1]='genename'
tab





#### sex overall
head(design)
contrast=c(0,0,0,1)


de2 =glmQLFTest(d3, contrast =contrast  )
test= topTags(de2, n=nrow(counts))
print(c(length(which(test$table$FDR<0.05))))

orths$gene_names_new[orths$ens_IDs%in%rownames(test$table[test$table$FDR<0.05,])]
tab=cbind(orths$gene_names_new[orths$ens_IDs%in%rownames(test$table[test$table$FDR<0.05,])],
          test$table[test$table$FDR<0.05,c(1,5)])

colnames(tab)[1]='genename'
tab
test[test$table$FDR<0.4,]

