options(stringsAsFactors = FALSE)


###
# assuming data has proper column names and row names
###


my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
normalizeDESeq2 = function ( rawCounts, benignId, malignantId, reportPath, groupNameA, groupNameB, saveReport = T)
{
  library(DESeq2)
  
  d = rawCounts
  sampleIds = colnames(rawCounts)
  condition <- factor(rep("A",dim(d)[2]))
  
  dds <- DESeqDataSetFromMatrix(countData = d, DataFrame(condition), design = ~ 1)
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  
  logs_norm_arr_ord = assay(vsd)
  
  logsBenign = logs_norm_arr_ord[,benignId]
  logsMalignant = logs_norm_arr_ord[,malignantId]
  
  write.table(logsBenign, file = paste(reportPath, groupNameA,"_", Sys.Date(), ".tsv", sep = ""), sep = "\t")
  write.table(logsMalignant, file = paste(reportPath, groupNameB, "_", Sys.Date(), ".tsv", sep = ""), sep = "\t")
  
  logsBenign = t(logsBenign)
  logsMalignant = t(logsMalignant)
  meanBenign <- sapply(1:dim(t(logs_norm_arr_ord))[2], function(i) mean(logsBenign[,i]))
  meanMalignant <- sapply(1:dim(t(logs_norm_arr_ord))[2], function(i) mean(logsMalignant[,i])) 
  
  medianBenign <- sapply(1:dim(t(logs_norm_arr_ord))[2], function(i) median(logsBenign[,i]))
  medianMalignant <- sapply(1:dim(t(logs_norm_arr_ord))[2], function(i) median(logsMalignant[,i])) 
  
  folds_mean_Benign_Malignant = 2^(meanBenign - meanMalignant)  
  ttest_Benign_Malignant_pvalue = sapply(1:dim(t(logs_norm_arr_ord))[2], function(i)
    my.t.test.p.value(logsBenign[,i], logsMalignant[,i]))
  
  fdrt_Benign_Malignant_qvalue = p.adjust(ttest_Benign_Malignant_pvalue, method = "fdr", n = length(ttest_Benign_Malignant_pvalue))
  summary_table_Benign_Malignant = cbind.data.frame(rownames(logs_norm_arr_ord), meanBenign, meanMalignant,
                                                    medianBenign, medianMalignant, folds_mean_Benign_Malignant,
                                                    ttest_Benign_Malignant_pvalue, fdrt_Benign_Malignant_qvalue )
  
  
  
  tableCols_Benign_Malignant=  cbind( "Gene",paste("Mean ", groupNameA, sep = ""), paste("Mean ", groupNameB, sep = ""),
                                      paste("Median ", groupNameA, sep = ""), paste("Median ", groupNameB, sep = ""),
                                      paste("Mean fold change ", groupNameA, "-", groupNameB, sep = ""),
                                      "T-test p-value ",
                                      "FDR q-value " ) 
  colnames(summary_table_Benign_Malignant) = tableCols_Benign_Malignant
  if (saveReport == T)
    write.table(summary_table_Benign_Malignant,file = paste(reportPath,  "Table_I_", groupNameA, "_vs_",
                                                            groupNameB, "_", Sys.Date(), ".tsv", sep = ""), 
                sep = "\t", row.names = TRUE)
  return(logs_norm_arr_ord[, c(benignId, malignantId)])
}


normalizeDESeq2NoReport = function ( rawCounts )
{
  library(DESeq2)
  
  d = rawCounts
  sampleIds = colnames(rawCounts)
  condition <- factor(rep("A",dim(d)[2]))
  
  dds <- DESeqDataSetFromMatrix(countData = d, DataFrame(condition), design = ~ 1)
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  
  logs_norm_arr_ord = assay(vsd)
  return(logs_norm_arr_ord) 
}


plotUmap = function(x, labels,
                    main="A UMAP visualization of the TEPs dataset",
                    colors=c("#ff7f00", "#e377c2", "#17becf"),
                    pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
                    cex.main=1, cex.legend=1, 
                    legend.pos = "topright") {
  library(umap)
  layout = x
  if (class(x)=="umap") {
    layout = x$layout
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}


RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
visualizeAnalysisData = function ( data, benignId, malignantId, reportPath, groupNameA, groupNameB, saveReport = T)
{
  
  # Filter 2 - remove flat genes with variation < 1st quartile 
  data = data[RowVar(data) > summary(RowVar(data))[2],]
  
  pdf(paste(reportPath, "Figure_I_histogram_", groupNameA,"_", Sys.Date(), ".pdf", sep = ""))
  dend_ec_oc <- hclust(dist(data))
  par(cex = 0.7)
  plot(dend_ec_oc, hang = -1, xlab =paste(groupNameA, " vs ", groupNameB, sep = "")  , sub = "")
  dev.off()
  
  pdf(paste(reportPath, "Figure_II_dendrogram_","_", Sys.Date(), ".pdf", sep = ""))#, width = 1024, height = 1024)
  hist(data)
  dev.off()
  
  pdf(paste(reportPath, "Figure_III_pca_","_", Sys.Date(), ".pdf", sep = ""))#, width = 1024, height = 1024)
  log.pca <- prcomp(t(data), center = TRUE, scale. = TRUE) 
  
  plot(log.pca$x[,c(1,2)], xlim = c(1.1*min(log.pca$x[,1])
                                    , 1.1*max(log.pca$x[,1])), ylim =c(1.1*min(log.pca$x[,1])
                                                                       , 1.1*max(log.pca$x[,1])))
  dev.off()
}