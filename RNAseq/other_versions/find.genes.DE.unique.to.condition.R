files=list.files('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE',pattern='.csv')
setwd('/home/oscar/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.DE')

ISGs=read.csv('/home/oscar/Documents/Alex_RNA_seq/Ov_ISGlist.csv')

dpi_comp='0.vs.1'
cond='Sassari.BTV8'

mat=matrix(nrow=100000,ncol=4,dat=NA)
count=1
for(i in files){
  DE=read.csv(i)
  mat[count:(count+nrow(DE)-1),]=cbind(DE$human_gene_names,
                                       DE$ensembl_ID,rep(i,nrow(DE)),DE[,4])
  count=count+nrow(DE)
}
mat=mat[1:which(is.na(mat[,3]))[1],]

matdpi1=mat[grepl(dpi_comp,mat[,3]),]
matdpi1_BTV8=matdpi1[grepl(cond,matdpi1[,3]),]
matdpi1_nonBTV8=matdpi1[!grepl(cond,matdpi1[,3]),]


hits=unique(matdpi1_BTV8[,2])[!matdpi1_BTV8[,2]%in%matdpi1_nonBTV8[,2]]

matdedup=mat[!duplicated(mat[,2]),]
matdedup[matdedup[,2]%in%hits,1]

orths=read.csv('/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
orths$gene_names_new[orths$ens_IDs%in%hits]

is.ISG=hits%in%ISGs$Gene.ID
output=cbind(orths$gene_names_new[match(hits,orths$ens_IDs)],
      hits,
      matdpi1_BTV8[match(hits,matdpi1_BTV8[,2]),4],is.ISG)

write.csv(output[order(output[,3],decreasing=T),],row.names=F,file=paste(
  '~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/unique_sets/',cond,'_',dpi_comp,'.csv',sep=''))
