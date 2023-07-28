dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
dat=gsub(' ','',dat)
dim(dat)

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
mat2_copy=mat2
mat2[1:4,1:4]

mat2_cor=cor(t(mat2))
hist(mat2_cor)

hits=which(abs(mat2_cor)>0.95,arr.ind=T)

for(i in 1:nrow(hits)){
  if(hits[i,1]>hits[i,2]){
    print(c(mat2_cor[hits[i,1],hits[i,2]],rownames(mat2_cor)[hits[i,2]],rownames(mat2_cor)[hits[i,1]]))
    mat2[hits[i,1],]=rep(NA,ncol(mat2))
    rownames(mat2)[hits[i,1]]=paste(rownames(mat2)[hits[i,1]],rownames(mat2)[hits[i,2]])
  }
}
mat2=mat2[!is.na(mat2[,1]),]

######################

BTMs_to_plot=c( "innate.antiviral.response..M150."                       
                , "antiviral.IFN.signature..M75."                          
                , "type.I.interferon.response..M127."                      
                , "enriched.in.activated.dendritic.cells..II...M165."      
                ,"Activated..LPS..dendritic.cell.surface.signature..S11." 
                , "myeloid..dendritic.cell.activation.via.NFkB..I...M43.0."
                ,"Ran.mediated.mitosis..M15."   )

mat2_cor[colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))],
         colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))]]

###
RSAD2_IFIT1=c("ENSOARG00000014648","ENSOARG00000015177")
counts[counts[,1]%in%RSAD2_IFIT1,"SMI8.292_b3_1_dpi_S36" ]
counts[match((dat[2:nrow(dat),"innate.antiviral.response..M150."])[!is.na(dat[,"innate.antiviral.response..M150."])],counts[,1]),"SMI8.292_b3_1_dpi_S36" ]

counts[match((dat[2:nrow(dat),"type.I.interferon.response..M127."])[!is.na(dat[,"type.I.interferon.response..M127."])],counts[,1]),"SMI8.292_b3_1_dpi_S36" ]

 

library(pheatmap)
pheatmap(mat2_cor[colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))],
                  colnames(mat2_cor)[match(BTMs_to_plot,colnames(mat2_cor))]]
         ,cluster_rows=F,cluster_cols=F)

########################################################################
BTMs_to_plot=c("TBA..M186."                                                           
               ,"enriched.in.monocytes..IV...M118.0."                                  
               ,"immune.activation...generic.cluster..M37.0."                          
               ,"growth.factor.induced..enriched.in.nuclear.receptor.subfamily.4..M94."
               ,"Naive.B.cell.surface.signature..S8." )


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

plots=list()
for(i in 1:5){
  BTM_to_plot=i
  to_plot=data.frame(vals=mat_plot[BTM_to_plot,],conds=colnames(mat_plot),treatment=treatments)
  to_plot$conds=factor(to_plot$conds,levels = cond_list_order)
  plots[[i]]=ggplot(to_plot,aes(x=conds,y=vals,fill=treatments))+
    geom_dotplot(binaxis='y',stackdir = 'center',alpha=0.3,dotsize = 1.8)+theme_bw()+
    ggtitle(rownames(mat_plot)[BTM_to_plot])+
    geom_vline(xintercept=(13.5))+
    geom_vline(xintercept=(18.5))+
    geom_vline(xintercept=(23.5))
    
}
library(gridExtra)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],ncol=1)



######################################
######################################
######################################
######################################


counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]

dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs.csv')

#####################
library(psych)
mat2=matrix(ncol=ncol(counts[2:ncol(counts)]),nrow=ncol(dat))
colnames(mat2)=colnames(counts[,2:ncol(counts)])
rownames(mat2)=names(dat)
for(i in 1:ncol(dat)){
  for(j in 2:ncol(counts)){
    mat2[i,j-1]=geometric.mean(1+counts[counts$Geneid%in%dat[2:nrow(dat),i],j])-1
  }
}


BTMs_to_plot=c( "innate.antiviral.response..M150."                       
               , "antiviral.IFN.signature..M75."                          
               , "type.I.interferon.response..M127."                      
               , "enriched.in.activated.dendritic.cells..II...M165."      
               ,"Activated..LPS..dendritic.cell.surface.signature..S11." 
               , "myeloid..dendritic.cell.activation.via.NFkB..I...M43.0."
               ,"Ran.mediated.mitosis..M15."   )


colnames(mat2_cor)

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
types_convert=list()
types_convert[['SMC']]=types_convert[['SLC']]='control'
types_convert[['SLI8']]=types_convert[['SMI8']]=paste('BTV8',sep='')
types_convert[['SMI13']]=types_convert[['SLI13']]=paste('2013',sep='')
types_convert[['SMI6']]=types_convert[['SLI6']]=paste('2006',sep='')


plots=list()
for(i in 1:length(BTMs_to_plot)){
  BTM_to_plot=i
  to_plot=data.frame(vals=mat_plot[BTM_to_plot,],conds=colnames(mat_plot),treatment=treatments)
  to_plot$conds=factor(to_plot$conds,levels = cond_list_order)
  to_plot$treatment_col=factor(types_convert[to_plot$treatment],levels=c('control','BTV8','2013','2006'))
  plots[[i]]=ggplot(to_plot,aes(x=conds,y=vals,fill=treatment_col))+
    geom_dotplot(binaxis='y',stackdir = 'center',alpha=0.3,dotsize = 1.3)+theme_bw()+
    ggtitle(rownames(mat_plot)[BTM_to_plot])+
    geom_vline(xintercept=(13.5))+
    geom_vline(xintercept=(18.5))+
    geom_vline(xintercept=(23.5))
  
}
library(gridExtra)
png(filename='/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/BTM.mach.learn/BTV8.dpi1.png',width=1100,height=1000)
  grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],ncol=1)
dev.off()



orths=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
overlap=matrix(ncol=7,nrow=7,dat=NA)
rownames(overlap)=rep(NA,7)
#colnames(overlap)=rep(NA,7)
for(i in 1:6){
  rownames(overlap)[i]=BTMs_to_plot[i]
  #colnames(overlap)[i]=BTMs_to_plot[i]
  for(j in (i+1):7){
  
    intersect_ij=intersect(mat[mat[,1]==BTMs_to_plot[i],2],mat[mat[,1]==BTMs_to_plot[j],2])
    overlap[i,j]=length(intersect_ij)
    
    print(c(i,j))
    print(intersect_ij)
    print(orths$gene_names_new[match(intersect_ij,orths$ens_IDs)])
  }
}
overlap

