library(scales)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(edgeR)
library(venn)
library(sva)
counts_in <- read.table("/home/oscar/Documents/sheep_megadata/RNA.seq.lab.nov.21/TruSeq-Nov2021/All_Count_test2.txt", header = T, row.names = 1)
colnames(counts_in) <- gsub("_20", "-20", gsub("B1", "BTV.1", gsub("B8", "BTV.8", colnames(counts_in))))
conditions <- sapply(colnames(counts_in), function(x) strsplit(x, "_")[[1]][1])
inf_states <- unique(grep("Mock", conditions, value = T, invert = T))
ISG_list_ov <- read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv", stringsAsFactors = F)
counts_in <- counts_in[!grepl("p6", colnames(counts_in))]
batches <- sapply(colnames(counts_in), function(x) {
  paste(strsplit(x, "_")[[1]][2],
    collapse = "_"
  )
})

gene_names <- read.csv("/home/oscar/Documents/BTV_RNA_seq/Ov_ISGlist.csv", stringsAsFactors = F)
all_gene_names <- read.csv("/home/oscar/Documents/orthologue_sets/sheep_and_cow_plus_manual_additions.csv")


outdir="/home/oscar/scripts/github/sheep_ML/outdir/RNA_seq/in_vitro_results/GLM.results.merge.tab.all.batch.combatseq"

times <- unique(sapply(colnames(counts_in), function(x) head(tail(strsplit(x, "_")[[1]], 2), 1)))[2:1]
time <- times[1]
inf <- "BTV.1-2006"
plots <- venns <- list()

counts <- counts_in

############################################################# main script

outdir=gsub("/$","",outdir)


# order columns
counts2 <- counts[, c(
  grep("Mock", colnames(counts)), grep(inf_states[1], colnames(counts)),
  grep(inf_states[2], colnames(counts)), grep(inf_states[3], colnames(counts))
)]
counts2 <- counts2[, c(grep("6h", colnames(counts)), grep("12h", colnames(counts)))]
# counts2=counts2[rowSums(counts2>0)>=(ncol(counts2)/2),]

conditions <- sapply(colnames(counts2), function(x) paste(strsplit(x, "_")[[1]][c(1, 3)], collapse = "_"))

conditions_notime <- sapply(colnames(counts2), function(x) paste(strsplit(x, "_")[[1]][c(1)], collapse = "_"))
plot(
  apply(counts2, 2, function(x) var(x / sum(x)))[c(
    grep("Mock", colnames(counts2)), grep(inf_states[1], colnames(counts2)),
    grep(inf_states[2], colnames(counts2)), grep(inf_states[3], colnames(counts2))
  )],
  col = as.numeric(factor(conditions_notime, levels = unique(conditions_notime))), pch = 19,
  ylab = "var(counts/sum(counts))", main = "6hr      \t\t\t\t          12hr",
  ylim = c(0, max(apply(counts2, 2, function(x) var(x / sum(x)))))
)
legend(x = "bottomright", pch = 19, col = 1:4, legend = unique(conditions_notime))
abline(v = 16.5)



treatments <- factor(conditions,
  levels = unique(conditions)
)
batches <- sapply(colnames(counts2), function(x) {
  paste(strsplit(x, "_")[[1]][2],
    collapse = "_"
  )
})
# batches=factor(batches,levels=unique(batches))

design <- model.matrix(~ 0 + treatments) #+batches)
rownames(design) <- colnames(counts)
conditions <- factor(conditions)
d <- DGEList(counts = counts2, group = treatments)

d$counts <- ComBat_seq(as.matrix(counts2), batch = batches, group = as.character(treatments))
d <- calcNormFactors(d)
d <- estimateDisp(d, design, robust = T)
#   d =estimateCommonDisp(d,design)
#    d =estimateTagwiseDisp(d,design)

d_all <- glmQLFit(d, design)


plot_dat_in <- matrix(nrow = 0, ncol = 5)
inf_states <- grep("Mock", unique(conditions), invert = T, value = T)
inf <- inf_states[1]
for (inf in inf_states) {
  cont <- c((unique(conditions) %in% inf))

  if (grepl("6h", inf)) {
    cont[1] <- -1
  } else {
    cont[5] <- -1
  }

  d_out <- glmQLFTest(d_all, contrast = cont)


  TT <- topTags(d_out, nrow(counts))
  TT_save <- cbind.data.frame(TT$table, all_gene_names$gene_name[match(rownames(TT$table), all_gene_names$sheep_ID)])

  write.csv(TT_save, paste(outdir,"/", inf, ".csv", sep = ""))
  print(c(
    inf, length(which(TT$table$FDR < 0.05)),
    length(which(TT$table$FDR < 0.05 & TT$table$logFC > 0)),
    length(which(TT$table$FDR < 0.05 & TT$table$logFC > 0 & rownames(TT$table) %in% ISG_list_ov$Gene.ID)),
    length(which(TT$table$logFC > 2)),
    length(which(TT$table$logFC > 2 & rownames(TT$table) %in% ISG_list_ov$Gene.ID))
  ))
  TT <- TT[rownames(TT$table) %in% ISG_list_ov$Gene.ID, ]
  plot_dat_in <- rbind(plot_dat_in, as.matrix(
    cbind(
      TT$table$logFC, TT$table$FDR < 0.05,
      rep(paste(inf, ".", sep = ""), nrow(TT)),
      rep(
        paste(inf, ".", "\nsig&logFC>2 = ",
          length(which(TT$table$FDR < 0.05 & abs(TT$table$logFC) > 2)), "",
          sep = ""
        ),
        nrow(TT)
      ),
      rownames(TT)
    )
  ))
}
colnames(plot_dat_in) <- c("logFC", "significant", "condition", "condition_lab", "gene")

unique(plot_dat_in[])


plots=list()
plot_i <- 0
plot_dat_save <- as.data.frame(plot_dat_in)

for (time in c("6h", "12h")) {
  plot_i <- plot_i + 1
  plot_dat <- plot_dat_save[grepl(time, plot_dat_save$condition), ]
  plot_dat$colour <- plot_dat$condition
  plot_dat$colour[!as.logical(plot_dat$significant)] <- "nonsig"
  plot_dat$colour <- factor(plot_dat$colour, levels = c("nonsig", unique(plot_dat$condition)))
  colpal <- c(scales::alpha("#AAAAAA", .2), scales::alpha(RColorBrewer::brewer.pal(length(unique(plot_dat$condition)), "Dark2"), .5))
  plot_dat$condition <- factor(plot_dat$condition, levels = unique(plot_dat$condition))

  order=c(grep("BTV.8-2017",unique(plot_dat$condition_lab),value=T),
      grep("BTV.1-2013",unique(plot_dat$condition_lab),value=T),
      grep("BTV.1-2006",unique(plot_dat$condition_lab),value=T))

  plot_dat$condition_lab <- factor(plot_dat$condition_lab, levels = order)
  plot_dat$logFC <- as.numeric(plot_dat$logFC)
  plot_dat <- plot_dat[order(plot_dat$colour), ]
  plots[[plot_i]] <- ggplot(plot_dat, aes(x = condition_lab, y = logFC, color = colour)) +
    geom_hline(yintercept = c(-2, 2)) +
    geom_jitter(size=1.7) +
    scale_colour_manual(values = colpal[sort(unique(as.numeric(plot_dat$colour)))]) +
    theme_bw() +
    theme(axis.text.x=element_text(size=22),axis.text.y=element_text(size=20),
       plot.title=element_text(size=24),
       axis.title.x=element_text(size=1),
       axis.title.y=element_text(size=22),
       legend.text=element_text(size=17))+
    ylim(-3, 8) +
    xlab("") +
    ggtitle(paste("ISG expression", time))+
    #drop legend
    guides(color=FALSE)

  plot_dat_sig <- plot_dat[plot_dat$significant == T, ]
  conds <- unique(plot_dat$condition)
  genes <- unique(plot_dat_sig$gene)
  pres_abs <- matrix(
    ncol = length(conds),
    nrow = length(genes), dat = 0
  )
  for (i in 1:nrow(pres_abs)) {
    for (j in 1:ncol(pres_abs)) {
      if (any(plot_dat_sig$gene[plot_dat_sig$condition == conds[j]] == genes[i])) {
        pres_abs[i, j] <- 1
      }
    }
  }

  rownames(pres_abs) <- genes
  pres_abs[rownames(pres_abs) %in% ISG_list_ov$Gene.ID, ]
  colnames(pres_abs) <- paste(conds, "\nN=", sapply(1:length(conds), function(x) {
    sum(pres_abs[, x])
  }), sep = "")
  venns[[time]] <- venn((as.data.frame((pres_abs),
    row.names = NULL,
    col.names = NULL,
  )), zcolor = "style", ggplot = T,ilcs = 1.7, sncs = 1.5) +
    geom_text(size=8,aes(size = 1000, x = 500, y = 1020, label = paste("DE ISGs on dpi", (time))))
  print(venns[[time]])

  pres_abs2 <- as.data.frame(pres_abs[order(rowSums(pres_abs), decreasing = T), ])
  pres_abs2$gene_names <- gene_names$Gene_name[match(rownames(pres_abs), gene_names$Gene.ID)]
  write.csv(pres_abs2, file = paste(outdir,"/combat_batch.all.tab",
    time, ".csv",
    sep = ""
  ))
}



png(paste(outdir,"/combat_batch.all.tab",
    ".venns.png",
    sep = ""), width=850,height=1050)
  grid.arrange(venns[["6h"]], venns[["12h"]], ncol = 1)
dev.off()

pdf(paste(outdir,"/combat_batch.all.tab",
    ".venns.pdf",
    sep = ""), width=8.50,height=10.50)
  grid.arrange(venns[["6h"]], venns[["12h"]], ncol = 1)
dev.off()

grid.arrange(plots[[1]], plots[[2]], ncol = 1)


png(paste(outdir,"/combat_batch.all.tab",
    ".dotplot.png",
    sep = ""), width=750,height=1000)
  grid.arrange(plots[[1]], plots[[2]], ncol = 1)
dev.off()

pdf(paste(outdir,"/combat_batch.all.tab",
    ".dotplot.pdf",
    sep = ""), width=7.50,height=10.00)
  grid.arrange(plots[[1]], plots[[2]], ncol = 1)
dev.off()


###############################################
###############################################
###############################################
###############################################
###############################################



inf_plot <- "BTV.8-2017"
inf_plot <- "BTV.1-2006"
inf_plot <- "BTV.1-2013"
par(mfrow = c(1, 2))
dat6 <- read.csv(paste("~/Documents/sheep_megadata/RNA.seq.lab.nov.21/GLM.results.merge.tab.all/", inf_plot, "_6h.csv", sep = ""))
dat12 <- read.csv(paste("~/Documents/sheep_megadata/RNA.seq.lab.nov.21/GLM.results.merge.tab.all/", inf_plot, "_12h.csv", sep = ""))
dat6 <- dat6[dat6$X %in% dat12$X, ]
dat12 <- dat12[dat12$X %in% dat6$X, ]
dat6 <- dat6[match(dat12$X, dat6$X), ]
col <- (dat12$FDR < 0.05) + 1
col[dat6$FDR < 0.05] <- 3
col[dat6$FDR < 0.05 & dat12$FDR < 0.05] <- 4
col2 <- alpha(col, (col / 4))
plot(dat6$logFC, dat12$logFC,
  col = col2, pch = 19, main = paste(inf_plot, "one combined 6+12 hr table"),
  xlab = "6hrs logFC", ylab = "12hrs logFC"
)
legend(x = "bottomright", col = 1:4, legend = paste(
  c("DE neither ;N=", "DE only at 12hr; N=", "DE only at 6hr;N=", "DE both ; N="),
  c(
    length(which(col == 1)), length(which(col == 2)),
    length(which(col == 3)), length(which(col == 4))
  )
), pch = 19, bty = "n")


dat6 <- read.csv(paste("~/Documents/sheep_megadata/RNA.seq.lab.nov.21/GLM.results.merge.tab/6h", inf_plot, ".csv", sep = ""))
dat12 <- read.csv(paste("~/Documents/sheep_megadata/RNA.seq.lab.nov.21/GLM.results.merge.tab/12h", inf_plot, ".csv", sep = ""))
dat6 <- dat6[dat6$X %in% dat12$X, ]
dat12 <- dat12[dat12$X %in% dat6$X, ]
dat6 <- dat6[match(dat12$X, dat6$X), ]
col <- (dat12$FDR < 0.05) + 1
col[dat6$FDR < 0.05] <- 3
col[dat6$FDR < 0.05 & dat12$FDR < 0.05] <- 4
col2 <- alpha(col, (col / 4))
plot(dat6$logFC, dat12$logFC,
  col = col2, pch = 19, main = paste(inf_plot, "two tables"),
  xlab = "6hrs logFC", ylab = "12hrs logFC"
)
legend(x = "bottomright", col = 1:4, legend = paste(
  c("DE neither ;N=", "DE only at 12hr; N=", "DE only at 6hr;N=", "DE both ; N="),
  c(
    length(which(col == 1)), length(which(col == 2)),
    length(which(col == 3)), length(which(col == 4))
  )
), pch = 19, bty = "n")
###########################################################################

#
# par(mar=c(4,4,1,1))
#
# dat6=read.csv("/home/oscar/Documents/sheep_megadata/RNA.seq.lab.nov.21/TruSeq-Nov2021/BTV2006_6h-VS-mock_6h_siggenes_edgeR_final.csv")
# dat12=read.csv("/home/oscar/Documents/sheep_megadata/RNA.seq.lab.nov.21/TruSeq-Nov2021/BTV2006_12h-VS-mock_12h_siggenes_edgeR_final.csv")
# dat6=dat6[dat6$ID%in%dat12$ID,]
# dat12=dat12[dat12$ID%in%dat6$ID,]
# dat6=dat6[match(dat12$ID,dat6$ID),]
# col=(dat12$FDR<0.05)+1
# col=alpha(col,(dat12$FDR<0.05)/3+.2)
# plot(dat6$logFC,dat12$logFC,col=col,pch=19)
