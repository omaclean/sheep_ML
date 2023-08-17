library('biomaRt')

SLI6.tru=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_5.7.21/TRUSEQSLI6dpi3.csv')
SLI6.3prim=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SLI6dpi3.csv')

sheepgenes=useMart("ensembl",dataset="oaries_gene_ensembl",host='www.ensembl.org')

genes=getBM(attributes=c("ensembl_gene_id",'entrezgene_id',"external_gene_name","refseq_mrna","uniprot_gn_symbol"),
            filters = 'ensembl_gene_id', values =SLI6.tru$X ,mart=sheepgenes)
genes=genes[!duplicated(genes[,1]),]

write(unique(genes[,5]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/tru.all.names.txt',ncol=1)
write(unique(genes[,2]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/tru.all.entrez.txt',ncol=1)


genes=getBM(attributes=c("ensembl_gene_id",'entrezgene_id',"external_gene_name","refseq_mrna","uniprot_gn_symbol"),
            filters = 'ensembl_gene_id', values =SLI6.tru$X[SLI6.tru$adj.P.Val<0.05&
                                                              SLI6.tru$logFC>0] ,mart=sheepgenes)
genes=genes[!duplicated(genes[,1]),]

write(unique(genes[,5]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/tru.sigup.names.txt',ncol=1)
write(unique(genes[,2]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/tru.sigup.entrez.txt',ncol=1)

#################


genes=getBM(attributes=c("ensembl_gene_id",'entrezgene_id',"external_gene_name","refseq_mrna","uniprot_gn_symbol"),
            filters = 'ensembl_gene_id', values =SLI6.3prim$X ,mart=sheepgenes)
genes=genes[!duplicated(genes[,1]),]

write(unique(genes[,5]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/prim.all.names.txt',ncol=1)
write(unique(genes[,2]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/prim.all.entrez.txt',ncol=1)


genes=getBM(attributes=c("ensembl_gene_id",'entrezgene_id',"external_gene_name","refseq_mrna","uniprot_gn_symbol"),
            filters = 'ensembl_gene_id', values =SLI6.3prim$X[SLI6.3prim$adj.P.Val<0.05&
                                                               SLI6.3prim$logFC>0] ,mart=sheepgenes)
genes=genes[!duplicated(genes[,1]),]

write(unique(genes[,5]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/prim.sigup.names.txt',ncol=1)
write(unique(genes[,2]),file='/home/oscar/Documents/sheep_megadata/test.pathway.analysis/prim.sigup.entrez.txt',ncol=1)

head(genes[,2])
