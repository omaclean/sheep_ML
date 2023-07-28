dat1=read.csv(file="~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/3.tier.GLM/Sassari.2013_BTV8.dpi.0.vs.7.csv")


 dat2=read.csv('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SMI8dpi7.csv')
dat1=dat1[dat1$ensembl_ID%in%dat2$X,]
dat2=dat2[match(dat1$ensembl_ID,dat2$X),]
par(mfrow=c(1,1))
  plot(dat2$logFC[match(dat1$ensembl_ID,dat2$X)],dat1$Sassari_SMI8_logFC)
  