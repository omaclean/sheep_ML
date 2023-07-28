library(biomaRt)

counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
gene_names=read.csv('/home/oscar/Documents/orthologue_sets/sheep_cow_and_human_esmbl_entrez_cow_esmbl.csv',stringsAsFactors = F)

gene_list=counts$Geneid
gene_names_new=gene_names$gene_name[match(gene_list,gene_names$sheep_ID)]
#httr::set_config(httr::config(ssl_verifypeer = FALSE))
sheepgenes=useMart("ensembl",dataset="oaries_gene_ensembl")

test_orthologues=getBM(attributes=c("ensembl_gene_id",
                                    "hsapiens_homolog_associated_gene_name","hsapiens_homolog_orthology_confidence","hsapiens_homolog_ensembl_gene"),
                       filters = 'ensembl_gene_id', values = gene_list,mart=sheepgenes)
test_orthologues$hsapiens_homolog_associated_gene_name
for(i in 1:length(gene_list)){
  if(is.na(gene_names_new[i])){
    if(gene_list[i]%in%test_orthologues$ensembl_gene_id){
      trials=test_orthologues$hsapiens_homolog_associated_gene_name[test_orthologues$ensembl_gene_id==gene_list[i]]
      if(any(trials!='')){
        gene_names_new[i]=trials[which(trials!='')[1]]
      }
    }
  }
}
extra_genes=read.csv('/home/oscar/Documents/orthologue_sets/Alex/Ov_fullmerge_final.csv')
for(i in 1:length(gene_list)){
  if(is.na(gene_names_new[i])){
    if(gene_list[i]%in%extra_genes$Gene.ID_ov){
      trials=c(extra_genes$Gene_name[extra_genes$Gene.ID_ov==gene_list[i]],
               extra_genes$Gene_name_human[extra_genes$Gene.ID_ov==gene_list[i]])
      if(any(trials!='')){
        gene_names_new[i]=trials[which(trials!='')[1]]
      }
    }
  }
}
extra_extra_genes=read.csv('/home/oscar/Documents/orthologue_sets/Extra.sheep.names.csv')
extra_extra_genes=extra_extra_genes[extra_extra_genes[,2]!='',]
for(i in 1:length(gene_list)){
  if(is.na(gene_names_new[i])){
    if(gene_list[i]%in%extra_genes$Ensembl.ID){
      gene_names_new[i]=extra_extra_genes[which(extra_extra_genes$Ensembl.ID==gene_list[i]),2]
    }}}
ens_IDs=gene_list
write.csv(cbind(gene_names_new,ens_IDs),'/home/oscar/Documents/orthologue_sets/Sheep_mega_data.list.csv')
