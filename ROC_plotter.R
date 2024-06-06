
library(ROCR)


outdir="/home/oscar/scripts/github/sheep_ML/outdir"

hits_on=10;hits_off=9;misses_on=30;misses_off=10;col=1;yaxt_var="n";xaxt_var="n";title_var="test"
ROC_curve_1=function(hits_on,misses_on,hits_off,misses_off,col=1,title_var="",xaxt_var="s",yaxt_var="s"){
    df=data.frame(predictions=c(rep(1,hits_on+misses_off),rep(0,misses_on+hits_off)),
        labels=c(rep(1,hits_on),rep(0,misses_off),rep(1,misses_on),rep(0,hits_off))
    )
    
    pred <- prediction(df$predictions, df$labels)
    perf <- performance(pred,"tpr","fpr")
    par(cex.main=1.6,cex.axis=1.4,cex.lab=1.4)
    plot(perf,col=col,lwd=2,
            main=title_var,
            yaxt=yaxt_var,
            xaxt=xaxt_var)
}
confusion_matrix=six_states
title="test"
run_ROC=function(confusion_matrix,title=""){

    for(j in 1:ncol(confusion_matrix)){
        if(j==ncol(confusion_matrix)){
               title_var=title
            xaxt_var="s"
            yaxt_var="s"
        }else{
         xaxt_var="n"
            title_var=""
            yaxt_var="n"
        }
        ROC_curve_1(confusion_matrix[j,j],
            sum(confusion_matrix[j,])-confusion_matrix[j,j],
            sum(confusion_matrix)-sum(confusion_matrix[j,])-sum(confusion_matrix[,j])+confusion_matrix[j,j],
            sum(confusion_matrix[,j])-confusion_matrix[j,j],
            col=j,
            title_var=title_var,
            yaxt_var=yaxt_var,
            xaxt_var=xaxt_var
        )
        par(new=T)
    }
    par(new=F)
    abline(a=0,b=1,lty=2)
    legend(x="bottomright",col=1:ncol(confusion_matrix),legend=colnames(confusion_matrix),bty="n",lwd=3,cex=2)

}


six_states=read.csv(paste0(outdir,"/six_states_performance_table_Nparams50.csv"),row.names=1)

four_states=read.csv(paste0(outdir,"/four_states_performance_table_Nparams17.csv"),row.names=1)

clinical=read.csv(paste0(outdir,"/clinical_RF_performance_table_100params.csv"),row.names=1)


png(paste0(outdir,"/six_states_performance_table_Nparams50_ROC.png"),width=600,height=600)
run_ROC(six_states,"Six states of infection")
dev.off()
pdf(paste0(outdir,"/six_states_performance_table_Nparams50_ROC.pdf"),width=9,height=9)
run_ROC(six_states,"Six states of infection")
dev.off()

png(paste0(outdir,"/four_states_performance_table_Nparams17_ROC.png"),width=600,height=600)
run_ROC(four_states,"Four states of infection")
dev.off()
pdf(paste0(outdir,"/sfour_states_performance_table_Nparams17_ROC.pdf"),width=9,height=9)
run_ROC(four_states,"Four states of infection")
dev.off()

colnames(clinical)=gsub("score.0.2","score 0-2",
                                gsub("score.6.8","score 6-8",
                                gsub("score.3.5","score 3-5",colnames(clinical))))
run_ROC(clinical,"Clinical score states of infection")


png(paste0(outdir,"/clinical_RF_performance_table_100params.png"),width=600,height=600)
run_ROC(clinical,"Clinical score states of infection")
dev.off()
pdf(paste0(outdir,"/clinical_RF_performance_table_100params.pdf"),width=9,height=9)
run_ROC(clinical,"Clinical score states of infection")
dev.off()
