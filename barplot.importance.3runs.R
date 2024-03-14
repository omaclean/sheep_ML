###############################################

library(venn)
library(smotefamily)
library(randomForest)
library(caret)
library(ggplot2)
library(plotly)
library(gridExtra)
library(RColorBrewer)
library(Amelia)
library(reshape2)
library(psych)
library(pheatmap)

outdir="/home/oscar/scripts/github/sheep_ML/outdir/barplots"


indir="/home/oscar/scripts/github/sheep_ML/outdir"
N_clinical=100
N_fourstates=17
N_sixstates=50
clin=read.csv(paste0(indir,'/clinical_score.data_',N_clinical,'.params.csv',sep=''))
four=read.csv(paste0(indir,"/final_four.data.",N_fourstates,'params.csv'))
six=read.csv(paste0(indir,"/final_six_states.data",N_sixstates,'params.csv'))




#functions


#tidy up names for plotting and consistency from sources
clean_names=function(x){
  x=gsub('BTM\\.','BTM:',x)
  x=gsub('combine\\.\\.','combine(',x)
  #x=gsub('combine\\(','comb(',x)
  x=gsub('\\.\\.','\\.',x)
  x=gsub('\\.\\.','\\.',x)
  x=substr(x,1,50)
  x[x=="IFN_y"|x=="IFN.y"]="IFN_y"
  return(x)
}

#################code
outdir=gsub("/$","",outdir)
#find params shared in all three models
shared_params=intersect(intersect(names(clin)[2:ncol(clin)],names(six)[2:ncol(six)]),
                        names(four)[2:ncol(four)])

write.csv(shared_params,paste0(outdir,'/shared_three_times_using',N_clinical,
          'params.csv'))

#find params shared in any combo of three models
shared_once_params=unique(c(intersect(names(clin)[2:ncol(clin)],names(six)[2:ncol(six)]),
                            intersect(names(six)[2:ncol(six)],names(four)[2:ncol(four)]),
                            intersect(names(clin)[2:ncol(clin)],names(four)[2:ncol(four)])
))

write.csv(shared_once_params,paste0(outdir,'/shared_at_least_once_params_using',N_clinical,
          '.csv'))
param_sets=list(clinical=names(clin)[2:ncol(clin)],six_states=names(six)[2:ncol(six)],
                four_states=names(four)[2:ncol(four)])



png(paste0(outdir,'/Venn.',N_clinical,
          'clinical_params.png'),width=700,height=700)
venn::venn(snames=c(paste('clinical score pred\nN =',ncol(clin)-1),
                    paste('six states pred\nN =',ncol(six)-1),
                    paste('four states pred\nN =',ncol(four)-1)),
           unname(param_sets),
           ilabels=T,
           zcolor = RColorBrewer::brewer.pal(7,'Set1'),ilcs=2,sncs=1.25 )
dev.off()

pdf(paste0(outdir,'/Venn.',N_clinical,
          'clinical_params.pdf'),width=10.00,height=10.00)
venn::venn(snames=c(paste('clinical score pred\nN =',ncol(clin)-1),
                    paste('six states pred\nN =',ncol(six)-1),
                    paste('four states pred\nN =',ncol(four)-1)),
           unname(param_sets),
           ilabels=T,
           zcolor = RColorBrewer::brewer.pal(7,'Set1'),ilcs=2,sncs=1.25 )
dev.off()

################################################
dpi_plot='7' #what day's clinical scores to use

# #read in RNA seq raw counts -not needed in this script
# counts=read.table('/home/oscar/Documents/sheep_megadata/RNA_Seq_1.12.20/All_Count.txt',header=T)
# #filter counts & BTMs so that gene is expressed in at least 2 conditions (on average)
# counts=counts[rowSums(counts[,2:ncol(counts)]>0)>14,]

#read in raw BTM params
mat2=read.csv('/home/oscar/Documents/sheep_megadata/30.6.21/dat_plots_filt/BTMs.final.csv',row.names = 1)
mat2[1:2,1:2]
#########################################



############################
setwd("/home/oscar/Documents/sheep_megadata/30.6.21/")

data=read.csv('combined/combined_sheet.csv',stringsAsFactors = F)
#filter out SMI13, non dpi 7 data and <80% complete variables
for(i in 3:ncol(data)){
  if(length(which(is.na(data[data$dpi==dpi_plot&grepl('SLC',data$ID),i])))==dpi_plot){
    data[data$dpi==dpi_plot&grepl('SLC',data$ID),i]=data[data$dpi==21&grepl('SLC',data$ID),i]
  }
  #if dpi 7 data missing fill in from neighbouring days
  if(length(which(is.na(data[data$dpi==dpi_plot&grepl('SLI6',data$ID),i])))==dpi_plot){
    for(j in grep('SLI6',data$ID,value=T)){
      data[data$dpi==dpi_plot&grepl(j,data$ID),i]=data[which(
        !is.na(data[,i])&
          data$dpi>dpi_plot&
          grepl(j,data$ID))[1]
        ,i]
    }
    data[data$dpi==dpi_plot&!grepl('SLC',data$ID),i]=data[,i]
  }
}

data=data[data$dpi==dpi_plot&!grepl('SLI13',data$ID),]


#remove variables with >20% missing data
data=data[,(colSums(is.na(data))<0.2*nrow(data))]

#set colours for plotting
cols=c("#20B2AA",brewer.pal(6,"Reds")[1+c(1,3,5,4,2)])


# pca_dat=amelia(data[,3:ncol(data)],logs=colnames(data)[3:ncol(data)])$imputations$imp1
# 
# pca_dat=amelia(data[,3:ncol(data)])$imputations$imp1

## make empty cells equal to the column mean
pca_dat=data[,3:ncol(data)]
for(i in 1:ncol(pca_dat)){
  nas=is.na(pca_dat[,i])
  nas[pca_dat[,i]==""]=T
  if(any(nas)){
    pca_dat[nas,i]=mean(as.numeric((pca_dat[,i])[!nas]))
  }
  pca_dat[,i]=as.numeric(pca_dat[,i])
}
#make rownames the animal ID and day
rownames(pca_dat)=paste(data$ID,data$dpi,sep='_')

metabol_dat_save=pca_dat
###################################################################

colnames(mat2)=  sapply(colnames(mat2),function(x) 
  gsub('\\.','-',paste(strsplit(x,'_')[[1]][c(1,3)],collapse='_')))
mat3=t(mat2)

pca_dat=pca_dat[rownames(pca_dat)%in%rownames(mat3),]
mat3=mat3[match(rownames(pca_dat),rownames(mat3)),]
#write.csv(pca_dat,'/home/oscar/Documents/sheep_megadata/13.4.21/dat_plots_filt/metabolic.final.csv')

################################################
#test line: comb_dat_funct_in=comb_dat;types_funct_convert=types_convert;xvar=50;clinkey=NA

gettop150=function(comb_dat_funct_in,types_funct_convert,clinkey=NA){
  comb_dat_funct=comb_dat_funct_in
  comb_dat_funct=comb_dat_funct[,!grepl('rectal|Albumin',colnames(comb_dat_funct))]
  
  if(!is.na(clinkey[1])){
    types_funct2=clinkey
  }else{
    types_funct=unlist(strsplit(rownames(comb_dat_funct),'-'))[(1:nrow(comb_dat_funct))*2-1]
    types_funct2=as.character(types_funct_convert[types_funct])
  }
  
  xvar2=150
  RF_1=randomForest(x=comb_dat_funct, y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
  hits=rownames(RF_1$importance)[order(RF_1$importance[,'MeanDecreaseGini'],decreasing = T)][1:xvar2]
  
  RF=randomForest(x=comb_dat_funct[,hits], y=as.factor(types_funct2),importance=T,do.trace = F,ntree=5000)
  
  gini= RF$importance[,'MeanDecreaseGini']
  
  print(c("SUM of all params gini importances::",sum(gini)))
  head(names(gini))
  top150=gini[order(gini,decreasing = T)[1:150]]
  return(top150)
}
###################################################


colnames(mat3)=paste('BTM:',colnames(mat3),sep='')
#drop highly correlated
#pca_dat=pca_dat[,sapply(colnames(pca_dat),function(x)!any(c('Conteggio.Globuli.Rossi..RBC.','GOT.AST.')==x))]#'TNF.alpha'



### make heatmaps

# #convert to Z-scores
# pheatdata=apply(mat3,2,function(x)(x-mean(x,na.rm=T))/sd(x,na.rm=T))
# write.csv(pheatdata,paste0(outdir,'/pheatdata_transformed_params.csv'))

# tail(colnames(comb_dat_heat))
# png(paste0(outdir,'/heatmap_shared_once_params_using',N_clinical,
#           'clinical_params.png'),width=700,height=700)
# colnames(pheatdata)=gsub('BTM\\.\\.','BTM:',colnames(pheatdata))
# colnames(pheatdata)=gsub('\\.\\.','\\.',colnames(pheatdata))
# colnames(pheatdata)=substr(colnames(pheatdata),1,30)

# pheatmap(pheatdata,cluster_cols=FALSE,cluster_rows=FALSE,NA_col='black')
# dev.off()

# rm(pheatdata)
# png(paste0(outdir,'/heatmap_shared_three_times_using',N_clinical,
#           'clinical_params.png'),width=700,height=700)
# pheatdata=comb_dat_heat[,match(clean_names(shared_params),clean_names(colnames(comb_dat_heat)))]

# pheatdata=apply(pheatdata,2,function(x)(x-mean(x,na.rm=T))/sd(x,na.rm=T))

# colnames(pheatdata)=gsub('BTM\\.\\.','BTM:',colnames(pheatdata))
# colnames(pheatdata)=gsub('\\.\\.','\\.',colnames(pheatdata))
# colnames(pheatdata)=substr(colnames(pheatdata),1,30)
# dim(pheatdata)
# pheatmap(pheatdata,cluster_cols=FALSE,cluster_rows=FALSE,NA_col='black')
# dev.off()






###processing for plotting:


comb_dat=as.data.frame(cbind(pca_dat,mat3))
comb_dat=comb_dat[,!grepl('rectal|Albumin',colnames(comb_dat))]
comb_dat=comb_dat[!grepl('SMI6-131|SMI6-132|SMI6-133',rownames(comb_dat)),]
colnames(comb_dat)=gsub("[[:blank:]]", "", colnames(comb_dat))
colnames(comb_dat)=gsub('\\.\\.','\\.',colnames(comb_dat))
colnames(comb_dat)=gsub('\\.\\.','\\.',colnames(comb_dat))
colnames(comb_dat)=gsub('\\)','',colnames(comb_dat))



types=unlist(strsplit(rownames(comb_dat),'-'))[(1:nrow(data))*2-1]
types_convert=list()
types_convert[['SMC']]='control_Sass'
types_convert[['SLC']]='control_Ter'
types_convert[['SMI8']]='BTV8_Sass'
types_convert[['SMI13']]='2013_Sass'
types_convert[['SLI13']]='2013_Ter'
types_convert[['SMI6']]='2006_Sass'
types_convert[['SLI6']]='2006_Ter'



top150_6=gettop150(comb_dat,types_convert)
######################################

types=unlist(strsplit(rownames(comb_dat),'-'))[(1:nrow(data))*2-1]
types_convert=list()
types_convert[['SMI6']]=types_convert[['SLI6']]='2006'
types_convert[['SMI13']]=types_convert[['SLI13']]='2013'
types_convert[['SLI8']]=types_convert[['SMI8']]='BTV8'
types_convert[['SMC']]=types_convert[['SLC']]='control'

top150_4=gettop150(comb_dat,types_convert)


##################################################

clinicals_in=read.csv('/home/oscar/Documents/sheep_megadata/clinical_score_Oscar.csv')
clinicals_in$ID_dpi=paste(clinicals_in$ID,clinicals_in$dpi,sep='_')

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
clinicals_discrete=rep('control',length(animals))
clinicals_discrete[!grepl('SMC|SLC',animals)]=
  unlist(convert_clinicals[as.character(clinicals[!grepl('SMC|SLC',animals)])])

top150_clin=gettop150(comb_dat,types_convert,clinicals_discrete)


##


#######################################
#######################################
#######################################
#######################################
for( run_nature in c("shared_once","fully_shared")){
  if(run_nature=="shared_once"){
    shared_param_set=shared_params
    ggdata=data.frame(parameters=shared_params)
  }else{
    shared_param_set=shared_once_params
    ggdata=data.frame(parameters=shared_once_params)
  }






  ggdata$parameters=clean_names(ggdata$parameters)

  #test if all parameters are present in all top150 lists
  if(sum(ggdata$parameters%in%clean_names(names(top150_6)))!=length(ggdata$parameters)&
    sum(ggdata$parameters%in%clean_names(names(top150_4)))!=length(ggdata$parameters)&
    sum(ggdata$parameters%in%clean_names(names(top150_clin)))!=length(ggdata$parameters)){
    stop('not all overlapping parameters matching in all 3, line 253')
  }


  ggdata$six_states=top150_6[match(ggdata$parameters,
    clean_names(names(top150_6)))]
  ggdata$four_states=top150_4[match(ggdata$parameters,
    clean_names(names(top150_4)))]
  ggdata$clinical=top150_clin[match(ggdata$parameters,
    clean_names(names(top150_clin)))]

  grep("RIG",clean_names(names(top150_clin)),value=T)
  grep("RIG",ggdata$parameters,value=T)

  knitr::kable(ggdata)

  #test presence of one BTM:
  match("BTM:combine(TBA.M137.TBA.M180..",substr(names(top150_clin),1,50))
  grep('TBA.M137',substr(names(top150_clin),1,50),value=T)
  ############
  new=melt(ggdata,id='parameters',value.name = 'importance')


  new$parameters=substr(new$parameters,1,30)
  #get summed importance, allowing for missing values from certain runs
  ggdata$sum=rep(0,nrow(ggdata))
  for(i in 1:nrow(ggdata)){
    if(!is.na(ggdata$clinical[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$clinical[i]}
    if(!is.na(ggdata$six_states[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$six_states[i]}
    if(!is.na(ggdata$four_states[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$four_states[i]}
    
  }
  new$parameters=factor(new$parameters,levels=new$parameters[order(ggdata$sum,decreasing = T)])
  write.csv(new,file=paste0(outdir,'/importance_',run_nature,'params_using',N_clinical,
          '_clinical_params.csv'))

  p=ggplot(new,aes(y=importance,fill=variable,x=parameters))+
    geom_bar(stat='identity')+theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
    ylab('gini importance')+
    ggtitle(paste('importance of',length(unique(new$parameters)),'parameters shared across all three'))

  ggsave(p,filename=paste0(outdir,'/importance_',run_nature,'params_using',N_clinical,
            '_clinical_params.png'),width=7,height=6)
  ggsave(p,filename=paste0(outdir,'/importance_',run_nature,'_params_using',N_clinical,
            '_clinical_params.pdf'),width=7,height=6)


  ## MAKE ordered heatmap

  param_order_for_plot=ggdata$parameters[order(ggdata$sum,decreasing = T)]



  mat4=mat3
  rownames(data)=paste(data$ID,data$dpi,sep='_')
  data2=data[rownames(data)%in%rownames(mat4),]
  mat4=mat4[match(rownames(data2),rownames(mat4)),]

  comb_dat_heat=as.data.frame(cbind(data2,mat4))
  print(dim(comb_dat_heat))

  #tidy up names
  colnames(comb_dat_heat)=gsub("[[:blank:]]", "", colnames(comb_dat_heat))
  colnames(comb_dat_heat)=gsub('\\.\\.','\\.',colnames(comb_dat_heat))
  colnames(comb_dat_heat)=gsub('\\.\\.','\\.',colnames(comb_dat_heat))
  colnames(comb_dat_heat)=gsub('\\)','',colnames(comb_dat_heat))


  #clean up names for plotting
  colnames(comb_dat_heat)=clean_names(colnames(comb_dat_heat))
  
  comb_dat_heat=comb_dat_heat[,param_order_for_plot[param_order_for_plot%in%clean_names(shared_param_set)]]
  comb_dat_heat=apply(comb_dat_heat,2,function(x) (log10(x+1) - mean(log10(x[grepl('SMC|SLC',rownames(comb_dat_heat))]+1)))/sd(log10(x[!is.na(x)]+1)))

  pheaty=comb_dat_heat
  write.csv(pheaty,paste0(outdir,'/pheaty_ordered.',run_nature,'params_using',N_clinical,
          '_clinical_params.csv'))

  pheaty=pheaty[!rownames(pheaty)%in% c("SMI6-131_7","SMI6-132_7","SMI6-133_7"),]


  rownames(pheaty)
  storage=rownames(pheaty)
  colnames(pheaty)=substr(colnames(pheaty),1,32)
  #rownames(pheaty)=paste(rownames(pheaty),clinicals_plot,sep=' #')
  clinicals_pheaty_plot=clinicals_in$clinical.score[match(storage,clinicals_in$ID_dpi)]
  
  #convert plotting names to more readable format & order them by spectrum of clinical score
  order_to_plot=c("Mock S","Mock T","BTV-8 S","BTV-8 T","BTV-1 2013 S","BTV-1 2013 T","BTV-1 2006 S","BTV-1 2006 T")

  rownames(pheaty)[grep('SMC',rownames(pheaty))]='Mock S'
  rownames(pheaty)[grep('SLC',rownames(pheaty))]='Mock T'
  rownames(pheaty)[grep('SMI6',rownames(pheaty))]='BTV-1 2006 S'
  rownames(pheaty)[grep('SLI6',rownames(pheaty))]='BTV-1 2006 T'
  rownames(pheaty)[grep('SMI13',rownames(pheaty))]='BTV-1 2013 S'
  rownames(pheaty)[grep('SLI13',rownames(pheaty))]='BTV-1 2013 T'
  rownames(pheaty)[grep('SMI8',rownames(pheaty))]='BTV-8 S'
  rownames(pheaty)[grep('SLI8',rownames(pheaty))]='BTV-8 T'
  
  store=rownames(pheaty)

  rownames(pheaty)=paste(rownames(pheaty),' CS:',clinicals_pheaty_plot,sep='')

  print(c(unlist(sapply(order_to_plot,function(x) which(grepl(x,store))))))
  
  pheaty=pheaty[c(unlist(sapply(order_to_plot,function(x) which(grepl(x,store))))),]


  colnames(pheaty)=substr(colnames(pheaty),1,32)
  library(grid);library(pheatmap)

  pal=c(colorRampPalette(c('#1842a5',"#3A5DB1",RColorBrewer::brewer.pal(7,'RdBu')[c(4,2)],
                          "#C4393B","#B2182B"
  ))(60))

  pal=pal[4:60]
  plot_fun=function(){
    pheatmap(t(pheaty),
            cluster_rows=F,cluster_cols=F,treeheight_row = 0, treeheight_col = 0,
            color= pal,na_col='#000000')
  }
  #print(plot_fun)
  png(paste0(outdir,'/heatmap_ordered.',run_nature,'pheatmap.png'),width=10.00,height=10.00,units='in',res=200)
  plot_fun()
  dev.off()
  pdf(paste0(outdir,'/heatmap_ordered.',run_nature,'pheatmap.pdf'),width=10.00,height=10.00)
  plot_fun()
  dev.off()

  
}


#######################################
#######################################
#######################################
#######################################
ggdata=data.frame(parameters=shared_once_params)


ggdata$parameters=clean_names(ggdata$parameters)

if(sum(ggdata$parameters%in%c(clean_names(names(top150_6)),clean_names(names(top150_4)),
                               clean_names(names(top150_clin))))!=length(ggdata$parameters)){
  stop('not all parameters matching in files, line 307')
}


ggdata$six_states=top150_6[match(ggdata$parameters,
  clean_names(names(top150_6)))]
ggdata$four_states=top150_4[match(ggdata$parameters,
  clean_names(names(top150_4)))]
ggdata$clinical=top150_clin[match(ggdata$parameters,
  clean_names(names(top150_clin)))]

############
new=melt(ggdata,id='parameters',value.name = 'importance')


new$parameters=substr(new$parameters,1,30)
ggdata$sum=rep(0,nrow(ggdata))
for(i in 1:nrow(ggdata)){
  if(!is.na(ggdata$clinical[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$clinical[i]}
  if(!is.na(ggdata$six_states[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$six_states[i]}
  if(!is.na(ggdata$four_states[i])){ggdata$sum[i]=ggdata$sum[i]+ggdata$four_states[i]}
  
}




new$parameters=factor(new$parameters,levels=substr(ggdata$parameters,1,30)[order(ggdata$sum,decreasing = T)])

#plot importance of parameters shared in any combination
write.csv(new,file=paste0(outdir,'/importance_allruns_params.csv'))

p2=ggplot(new,aes(y=importance,fill=variable,x=parameters))+
  geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
  ylab('gini importance')+
  ggtitle(paste('importance of',length(unique(new$parameters)),'parameters shared in any combination'))

ggsave(p2,filename=paste0(outdir,'/importance_any_pair_shared_params_using',N_clinical,
          'clinical_params.png'),width=10.00,height=7.00)
ggsave(p2,filename=paste0(outdir,'/importance_any_pair_shared_params_using',N_clinical,
          'clinical_params.pdf'),width=10.00,height=7.00)


# 
# 
# 
# new2=new[new$variable=='four_states',]
# 
# 
# new2$parameters=factor(new2$parameters,levels=new2$parameters[rev(order(new2$importance))])
# 
# new2=new2[!is.na(new2$importance),]
# 
# ggplot(head(new2,10),aes(y=importance,fill=variable,x=parameters))+
#   geom_bar(stat='identity')+theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10.5))+
#   ylab('gini importance')+
#   ggtitle('importance of 40 shared in any combo parameters')
# 
# 
