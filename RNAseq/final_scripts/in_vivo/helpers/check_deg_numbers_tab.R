#check_deg_numbers

dir="/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/"

files=list.files(dir,pattern=".csv")
outdir="/home/oscar/scripts/github/sheep_ML/RNAseq/outputs"


#################
dir=paste0(gsub("/$","",dir),"/")
days=gsub("\\.csv","",unique(sapply(files,function(x) strsplit(x,'dpi')[[1]][2])))

conditions=unique(sapply(files,function(x) strsplit(x,'dpi')[[1]][1]))
days
conditions
files
tab=matrix(nrow=length(conditions),ncol=length(days))
rownames(tab)=conditions
colnames(tab)=days
for(condition_i in conditions){
    for(day_i in days){
        file=read.csv(paste0(dir,condition_i,'dpi',day_i,'.csv'))
        head(file)
        DEG=sum(file$adj.P.Val<0.05)
        tab[condition_i,day_i]=DEG
    }
}
knitr::kable(tab)

write.csv(tab,file=paste0(outdir,'/','in_vivo_dreams_deg_numbers.csv',sep=''))
