setwd('~/Pictures/plots/Sheep_megadata/RNA_seq/init_analyses/big_loop/sep.all')

cont=read.csv('Teramo.control.dpi.0.vs.3.csv')
thirt=read.csv('Teramo.2013.dpi.0.vs.3.csv')
six=read.csv('Teramo.2006.dpi.0.vs.3.csv')

common=intersect(intersect(thirt$ensembl_ID,six$ensembl_ID),cont$ensembl_ID)
vals=matrix(ncol=4,nrow=length(common))
vals[,1]=common
vals[,2]=as.numeric(cont[match(common,cont$ensembl_ID),4])
vals[,3]=as.numeric(thirt[match(common,thirt$ensembl_ID),4])
vals[,4]=as.numeric(six[match(common,six$ensembl_ID),4])

par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(vals[,2],vals[,3])
plot(vals[,2],vals[,4])
plot(vals[,3],vals[,4])

cor(as.numeric(vals[,2]),as.numeric(vals[,3]))
cor(as.numeric(vals[,2]),as.numeric(vals[,4]))
cor(as.numeric(vals[,3]),as.numeric(vals[,4]))


####################


cont=read.csv('Sassari.control.dpi.0.vs.3.csv')
thirt=read.csv('Sassari.2013.dpi.0.vs.3.csv')
six=read.csv('Sassari.2006.dpi.0.vs.3.csv')

common=intersect(intersect(thirt$ensembl_ID,six$ensembl_ID),cont$ensembl_ID)
vals=matrix(ncol=4,nrow=length(common))
vals[,1]=common
vals[,2]=as.numeric(cont[match(common,cont$ensembl_ID),4])
vals[,3]=as.numeric(thirt[match(common,thirt$ensembl_ID),4])
vals[,4]=as.numeric(six[match(common,six$ensembl_ID),4])

par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(vals[,2],vals[,3])
plot(vals[,2],vals[,4])
plot(vals[,3],vals[,4])

cor(as.numeric(vals[,2]),as.numeric(vals[,3]))
cor(as.numeric(vals[,2]),as.numeric(vals[,4]))
cor(as.numeric(vals[,3]),as.numeric(vals[,4]))



thirt=read.csv('Sassari.2013.dpi.0.vs.7.csv')
quan=read.csv('/home/oscar/Downloads/SMI13dpi0-VS-SMI13dpi7_siggenes_edgeR_final.csv',header=F)
head(quan)


common2=intersect(quan[,1],thirt$ensembl_ID)
vals=matrix(ncol=3,nrow=length(common2))
vals[,1]=common2
vals[,2]=as.numeric(thirt[match(common2,thirt$ensembl_ID),4])
vals[,3]=as.numeric(quan$V2[match(common2,quan[,1])])

par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
plot(vals[,2],vals[,3])
head(vals)
summary(lm(as.numeric(vals[,2])~as.numeric(vals[,3])))
