dat=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs.csv')
for(i in 1:ncol(dat)){
  
  dat[,i]=gsub(' ','',dat[,i])  
  dat[duplicated(dat[,i]),i]=NA
}

write.csv(dat,row.names=F,file='/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/BTMs_from_Artur_deduplicated_with_NAs2.csv')
