dream=read.csv('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/3.tier.GLM/Teramo.2013_2006.dpi.0.vs.7.csv')
old=read.csv('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/Teramo.2006.dpi.0.vs.7.csv')
old=read.csv('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all/Teramo.2013.dpi.0.vs.7.csv')
dream=read.csv('~/Documents/sheep_megadata/RNA_Seq_1.12.20/mixed_effects/controls_infect_dpi0_other/SLI6dpi7.csv')

old=old[old[,2]%in%dream[,1],]

dream=dream[match(old[,2],dream[,1]),]
head(dream)
head(old)
plot(dream[,2],old[,4])

length(which(dream$adj.P.Val<0.05))
length(which(old[,3]<0.05))

length(which(dream$adj.P.Val<0.05&old[,3]<0.05))

colnames(dream)
colnames(old)
