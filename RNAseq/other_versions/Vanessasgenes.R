
#####
bigloves=readLines('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/big_love-genes/new.explanation.csv')

bigloves2=read.csv('/home/oscar/Downloads/for_Vanessa2.csv')
bigloves_table=data.frame(gene_name=NULL,ens_ID=NULL)
for(i in bigloves){
  genes=strsplit(i,',')[[1]]
  for(j in 2:length(genes)){
    bigloves_table=rbind(bigloves_table,c(paste(genes[1],'-',letters[j-1]),genes[j]))  
  }
}
colnames(bigloves_table)=c('gene_name','Ens_ID')
for(i in 1:nrow(bigloves2)){
  if(!bigloves2[i,1]%in%bigloves_table$Ens_ID){
    bigloves_table=rbind(bigloves_table,c(bigloves2[i,2],bigloves2[i,1]))
  }
}
head(bigloves_table)
dir='~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all'
setwd(dir)
files=c(files[grepl('Sassari.control',files)],
        files[grepl('Teramo.control',files)],
        files[grepl('BTV8',files)],
        files[grepl('Sassari.2013',files)],
        files[grepl('Teramo.2013',files)],
        files[grepl('Sassari.2006',files)],
        files[grepl('Teramo.2006',files)]
)
logFCs=matrix(ncol=length(files),nrow=nrow(bigloves_table))
for(i in 1:length(files)) {
  dat=read.csv(files[i],stringsAsFactors = F)
  logFCs[,i]=dat[match(bigloves_table$Ens_ID,dat$ensembl_ID),4]
}
rownames(logFCs)=bigloves_table$gene_name
colnames(logFCs)=sapply(files,function(x) gsub('.dpi','',strsplit(x,'.csv')[[1]][1]))
colnames(logFCs)=sapply(colnames(logFCs),function(x) paste(substr((strsplit(x,'\\.')[[1]]),1,4),collapse='.'))
pheatmap(logFCs,
         cluster_rows=F,cluster_cols=F,treeheight_row = 0, treeheight_col = 0)

# length(dat$ensembl_ID)
# print(bigloves_table$Ens_ID%in%dat$ensembl_ID)
# print(bigloves_table$Ens_ID%in%
#         read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt')[,1])
# bigloves_table$in_RNA_Seq_data=(bigloves_table$Ens_ID%in%
#                                   read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt')[,1])
 write.csv(bigloves_table,file='/home/oscar/Downloads/All_comb_van.table.csv',row.names = F)
# ###########