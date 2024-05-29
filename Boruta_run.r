


library(randomForest)
library(Boruta)

test_data=read.csv("/home/oscar/scripts/github/sheep_ML/outdir/six_states_RF_data.csv")

classes=sapply(test_data$X,function(x) strsplit(x,"-")[[1]][1])


types_convert=list()
types_convert[['SMC']]='control_Sass'
types_convert[['SLC']]='control_Ter'
types_convert[['SMI8']]='BTV8_Sass'
types_convert[['SMI13']]='2013_Sass'
types_convert[['SLI13']]='2013_Ter'
types_convert[['SMI6']]='2006_Sass'
types_convert[['SLI6']]='2006_Ter'
outdir="/home/oscar/scripts/github/sheep_ML/outdir/Boruta/sixstates"

types=factor(as.character(types_convert[classes]),levels=unlist(types_convert)[unlist(types_convert)%in%as.character(types_convert[classes])])

test_data_reg=as.data.frame(test_data[,2:ncol(test_data)])

set.seed(42);test_dat=Boruta(x=test_data_reg , y=types,doTrace=T,maxRuns=ncol(test_data_reg))
# results=(attStats(test_dat))

# knitr::kable(results)


confirmed=rownames(results)[results$decision=="Confirmed"]
tentative=rownames(results)[results$decision=="Tentative"]

six_states_comb=c(confirmed,tentative)
six_states_conf=confirmed

write.csv(confirmed,paste0(outdir,"/confirmed_params.csv"))
write.csv(tentative,paste0(outdir,"/tentative_params.csv"))

print(paste("length(confirmed)",length(confirmed),"length(tentative)",length(tentative)))
six_state_50params=read.csv("/home/oscar/scripts/github/sheep_ML/outdir/final_six_states.data50params.csv")
sum(colnames(six_state_50params)%in%confirmed)
sum(colnames(six_state_50params)%in%c(confirmed,tentative))


################################################################################################################################################
################################################################################################
################################################four states
classes=sapply(test_data$X,function(x) strsplit(x,"-")[[1]][1])

types_convert=list()
types_convert[['SMI6']]=types_convert[['SLI6']]='2006'
types_convert[['SMI13']]=types_convert[['SLI13']]='2013'
types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8'
types_convert[['SMC']]=types_convert[['SLC']]='control'


outdir="/home/oscar/scripts/github/sheep_ML/outdir/Boruta/fourstates"
types=factor(as.character(types_convert[classes]),levels=unique(unlist(types_convert)[unlist(types_convert)%in%as.character(types_convert[classes])]))


test_data_reg=as.data.frame(test_data[,2:ncol(test_data)])

set.seed(42);test_dat=Boruta(x=test_data_reg , y=types,doTrace=T,maxRuns=ncol(test_data_reg))
results=(attStats(test_dat))



confirmed=rownames(results)[results$decision=="Confirmed"]
tentative=rownames(results)[results$decision=="Tentative"]

four_states_comb=c(confirmed,tentative)
four_states_conf=c(confirmed)

print(paste("length(confirmed)",length(confirmed),"length(tentative)",length(tentative)))
write.csv(confirmed,paste0(outdir,"/confirmed_params.csv"))
write.csv(tentative,paste0(outdir,"/tentative_params.csv"))

four_state_17params=read.csv("/home/oscar/scripts/github/sheep_ML/outdir/final_four.data.17params.csv")
sum(colnames(four_state_17params)%in%confirmed)
sum(colnames(four_state_17params)%in%c(confirmed,tentative))


################################################################################################################################################
################################################################################################
################################################clinical
outdir="/home/oscar/scripts/github/sheep_ML/outdir/Boruta/clinical"

clinicals_in=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv')
clinicals_in=clinicals_in[clinicals_in$dpi==dpi_plot,]
#classes_N=30

animals=unlist(strsplit(rownames(comb_dat),'_'))[(1:nrow(comb_dat))*2-1]


clinicals_in=clinicals_in[clinicals_in$ID%in%animals,]
comb_dat=comb_dat[order(clinicals_in$clinical.score[match(animals,clinicals_in$ID)]),]
clinicals=clinicals_in$clinical.score[match(animals,clinicals_in$ID)]
names(clinicals)=animals # bad coding makes this necessary for function


convert_clinicals=list()
convert_clinicals[['0']]=convert_clinicals[['1']]=convert_clinicals[['2']]='inf clin_score 0-2'
convert_clinicals[['3']]=convert_clinicals[['4']]=convert_clinicals[['5']]='inf clin_score 3-5'
convert_clinicals[['6']]=convert_clinicals[['7']]=convert_clinicals[['8']]='inf clin_score 6-8'
groups_hierarchy=c("control","inf clin_score 0-2","inf clin_score 3-5","inf clin_score 6-8")
clinicals_discrete=rep('control',length(animals))
clinicals_discrete[!grepl('SMC|SLC',animals)]=
  unlist(convert_clinicals[as.character(clinicals[!grepl('SMC|SLC',animals)])])

clinicals_discrete=factor(clinicals_discrete,levels=c("control",(unique(convert_clinicals))))
set.seed(42);test_dat=Boruta(x=test_data_reg , y=clinicals_discrete,doTrace=T,maxRuns=ncol(test_data_reg))
results=(attStats(test_dat))



confirmed=rownames(results)[results$decision=="Confirmed"]
tentative=rownames(results)[results$decision=="Tentative"]


print(paste("length(confirmed)",length(confirmed),"length(tentative)",length(tentative)))

clin_comb=c(confirmed,tentative)
clin_conf=c(confirmed)
write.csv(confirmed,paste0(outdir,"/confirmed_params.csv"))
write.csv(tentative,paste0(outdir,"/tentative_params.csv"))


clin_100params=read.csv("/home/oscar/scripts/github/sheep_ML/outdir/clinical_score.data_100.params.csv",row.names=1)
sum(colnames(clin_100params)%in%confirmed)
sum(colnames(clin_100params)%in%c(confirmed,tentative))



clin_30params=read.csv("/home/oscar/scripts/github/sheep_ML/outdir/clinical_score.data_30.params.csv",row.names=1)
sum(colnames(clin_30params)%in%confirmed)
sum(colnames(clin_30params)%in%c(confirmed,tentative))


intersect(clin_comb,intersect(six_states_comb,four_states_comb))
intersect(clin_conf,intersect(six_states_conf,four_states_conf))
####total overlap


