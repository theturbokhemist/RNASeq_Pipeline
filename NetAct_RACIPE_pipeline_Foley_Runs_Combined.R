# read counts data
setwd("/Users/gordid/Desktop/Foley_Runs_Combined")

SkipCalcs = TRUE


RNAseqDegs_limma_dge <- function(dge, design, contr.matrix, pval = 0.05){
  require(edgeR)
  require(limma)
  # celltype = phenodata$celltype
  # levs = levels(celltype)
  # # if(length(levs) > 2) print("Warning ... more than two types!")
  # compList = paste(levs, collapse = "-")
  

  # if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
  #   if_batch = FALSE
  #   design = model.matrix(~ 0 + celltype)
  #   colnames(design) = levs
  # }else{
  #   if_batch = TRUE
  #   batch = phenodata$batch
  #   levs_batch = levels(batch)
  #   design = model.matrix(~ 0 + celltype + batch)
  #   colnames(design) = c(levs, levs_batch[2:length(levs_batch)])
  # }
  # contr.matrix = makeContrasts(contrasts = compList, levels = design)
  v <- voom(dge, design, plot=TRUE)
  vfit = lmFit(v, design)
  vfit = contrasts.fit(vfit, contrasts=contr.matrix)
  efit = eBayes(vfit)
  results = decideTests(efit, p.value = pval)
  print(summary(results))
  rslt = topTable(efit, number = Inf, sort.by = "P")
  rslt$padj = rslt$adj.P.Val
  rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
  degs = rownames(rslt[rslt$padj< qval, ])
  
  # if (if_batch) {
  #   design1 = model.matrix(~ 0 + celltype)
  #   data.rm.batch <- removeBatchEffect(v$E, batch, design=design1)  # limma remove batch
  # }else{
     data.rm.batch = NULL
  # }
  
  return(list(table = rslt, rank_vector = rank_vector, degs = degs, 
              e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch)))
}

if (TwoGroups == FALSE) {
phenoData_Foley_Runs_Combined = new("AnnotatedDataFrame", data = data.frame(celltype_Foley_Runs_Combined = group_name_Foley_Runs_Combined))
rownames(phenoData_Foley_Runs_Combined) = sample_name_Foley_Runs_Combined
phenoData_Foley_Runs_Combined
}

# DE limma + voom
library(NetAct)
load("mgs.rdata")
load("hdb_v2.rdata")

if (Mouse == TRUE) {
  specimen <- mgs
} else {
  specimen <- hDB_v2  
}

compare_list <- c()
for (i in 1:(length(levels(group_name_Foley_Runs_Combined))-1)) {
  for (j in (i+1):length(levels(group_name_Foley_Runs_Combined))) {
    compare_list <- append(compare_list, paste0(levels(group_name_Foley_Runs_Combined)[i], "-", levels(group_name_Foley_Runs_Combined)[j]))
    
  }
}
compare_list

if (TwoGroups == TRUE) {
  DErslt_test <- RNAseqDegs_limma_dge(combined_data_Foley_Runs_Combined, design_matrix_Foley_Runs_Combined_no.int, contr.matrix_no.int)
  str(DErslt_test)
  # DErslt = RNAseqDegs_limma(Combined_data_Foley_Runs_Combined_NetAct$counts, phenoData_Foley_Runs_Combined)
  # str(DErslt)
  
  e = DErslt$e
  
  neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_Foley_Runs_Combined)
  
  save(neweset, DErslt, file = paste(Directory_CSVs_NetAct_RACIPE, "/macro_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".RData", sep = ""))
  
  write.csv(DErslt$table, file = paste(Directory_CSVs_NetAct_RACIPE, "/DEGenes_Foley_Runs_Combined_RACIPE", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                       filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))
  
  if (SkipCalcs == FALSE) {
    gsearslt_Foley_Runs_Combined <- TF_GSEA(specimen, DErslt, minSize = 5, nperm = 1000, qval = T)
    write.csv(gsearslt_Foley_Runs_Combined, file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                           filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),".csv", sep = ""))
  }
  setwd(Directory_CSVs_NetAct_RACIPE)
  gsearslt_Foley_Runs_Combined <- read.csv(file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                          filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))
  qval <- 0.7
  Tfs_Foley_Runs_Combined <- as.character(gsearslt_Foley_Runs_Combined$tf[gsearslt_Foley_Runs_Combined$qvals < qval])
  tfs <- sort(unique(as.character(c(Tfs_Foley_Runs_Combined))))
  acts_mat = TF_Activity(tfs, specimen, neweset, DErslt, useDatabaseSign = FALSE)$all_activities
  test <- TF_Activity(tfs, specimen, neweset, DErslt, useDatabaseSign = FALSE)
  test
  
} else {
  DErslt = MultiRNAseqDegs_limma(Combined_data_Foley_Runs_Combined_NetAct$counts, phenoData_Foley_Runs_Combined, compare_list)
  str(DErslt)
  e = DErslt$Overall$e
  
  neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_Foley_Runs_Combined)
  
  save(neweset, DErslt, file = paste(Directory_CSVs_NetAct_RACIPE, "/macro_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".RData", sep = ""))
  
  NetAct_gene_list_Foley_Runs_Combined <- list()
  gsearslt_list <- list()
  tfs_list <- list()
  TFs_vector <- c()
  for (i in 1:length(compare_list)) {
    NetAct_gene_list_Foley_Runs_Combined[[i]] <- DErslt[[i]]$table
    write.csv(NetAct_gene_list_Foley_Runs_Combined[[i]], file = paste(Directory_CSVs_NetAct_RACIPE, "/DEGenes_Foley_Runs_Combined_RACIPE", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                        filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))
    
    if (SkipCalcs == FALSE) {
      gsearslt_Foley_Runs_Combined_multi <- TF_GSEA(specimen, DErslt[[i]], minSize = 5, nperm = 10000, qval = T)
      write.csv(gsearslt_Foley_Runs_Combined_multi, file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs_Combined", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                   filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),".csv", sep = ""))
    }
    setwd(Directory_CSVs_NetAct_RACIPE)
    gsearslt_list[[i]] <- read.csv(paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs_Combined", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                         filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))
    
    qval <- 0.05
    tfs_list[[i]] <- as.character(gsearslt_list[[i]]$tf[gsearslt_list[[i]]$qvals < qval])
    TFs_vector <- append(TFs_vector, tfs_list[[i]])
  }
  names(NetAct_gene_list_Foley_Runs_Combined) <- Comparisons_vector
  names(gsearslt_list) <- Comparisons_vector
  names(tfs_list) <- Comparisons_vector
  
  tfs <- sort(unique(as.character(TFs_vector)))
  lengths(tfs_list)
  length(tfs)
  tfs
  acts_mat = TF_Activity(tfs, specimen, neweset, DErslt$Overall)$all_activities
}
setwd("/Users/gordid/Desktop/Foley_Runs_Combined")

if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_NetAct_RACIPE, "/TFsHeatMap_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), "_qval", 
                                         gsub(".", "", as.character(qval), fixed = TRUE), filter_type, gsub(".*filter", "", Directory_Plots_NetAct_RACIPE),
                                         ".pdf", sep = ""), width = 15, height = 15)}

Combine_heatmap(acts_mat, neweset) #Heatmap

if (SavePlots == TRUE) {dev.off()}

tf_links <- TF_Filter(acts_mat, specimen, miTh = 0, nbins = 3)
tf_links
miTh = 0.0
nbins = 3
dim(tf_links)

####Generating Network####
if (SavePlots == TRUE) {
  setwd(Directory_Plots_NetAct_RACIPE)
  plot_network_v(tf_links)
  # webshot::webshot(paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs_Combined_html", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
  #                        "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
  #                        filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),".html", sep = ""),
  #                  file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
  #                               "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
  #                               filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".png", sep = ""), delay = 4,  vwidth = 1500, vheight = 1500)
  # setwd("/Users/gordid/Desktop/Foley_Runs_Combined")
  # plot_network_v(tf_links)
} else {
  plot_network_v(tf_links)
}

x <- union(tf_links$Source,tf_links$Target)



#if (SavePlots == TRUE) {
# orca(plot_network_v(tf_links), file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
#    "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
#                              if (outliers == TRUE) {print("_outliers-removed")}, ".png", sep = "")) 
#} else {
#  plot_network_v(tf_links)
#}

write.csv(tf_links, file = paste(Directory_CSVs_NetAct_RACIPE, "/TFNetwork_InteractionsTable_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                 "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                 filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))

#write.table(tf_links, "net_emt.tpo", quote = FALSE, row.names = FALSE)
write.table(tf_links, file = paste(Directory_Plots_NetAct_RACIPE, "/net_emt_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                   "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                   filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".tpo", sep = ""), quote = FALSE, row.names = FALSE)

####PERTURBATION####
if (SavePlots == TRUE) {
  setwd(Directory_Plots_NetAct_RACIPE)
  library(sRACIPE)
  
  racipe <- sRACIPE::simulateGRC(paste(Directory_Plots_NetAct_RACIPE, "/net_emt_Foley_Runs_Combined", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                       "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                       filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".tpo", sep = ""), quote = FALSE, row.names = FALSE, plots = TRUE, plotToFile = TRUE)
  racipePer <- sRACIPE::knockdownAnalysis(racipe,reduceProduction = 10)
  setwd("/Users/gordid/Desktop/Foley_Runs_Combined")
} else {
  print("No Perturbation")
}
#gsearslt_1 = read.csv(paste("Foley_Runs_Combined_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = read.csv(paste("Foley_Runs_Combined_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = read.csv(paste("Foley_Runs_Combined_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#tfs_1 = as.character(gsearslt_1$tf[gsearslt_1$qvals <0.2])
#tfs_2 = as.character(gsearslt_2$tf[gsearslt_2$qvals <0.2])
#tfs_3 = as.character(gsearslt_3$tf[gsearslt_3$qvals <0.2])
#tfs = sort(unique(as.character(c(tfs_1, tfs_2, tfs_3))))
#lengths(list(tfs_1,tfs_2,tfs_3, tfs))
#tfs
# for normal use ####
#e = DErslt$Overall$e
#neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_Foley_Runs_Combined)
#save(neweset, DErslt, file = "macro.RData")

#Lupus_Foley_Runs_Combined_BvC_genes <- (DErslt$`B-C`$table)
#head(Lupus_Foley_Runs_Combined_BvC_genes)

#write.csv(Lupus_Foley_Runs_Combined_BvC_genes, file = paste("Foley_Runs_Combined_BvC_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#Lupus_Foley_Runs_Combined_AvC_genes <- (DErslt$`A-C`$table)
#write.csv(Lupus_Foley_Runs_Combined_BvC_genes, file = paste("Foley_Runs_Combined_AvC_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#Lupus_Foley_Runs_Combined_AvB_genes <- (DErslt$`A-B`$table)
#write.csv(Lupus_Foley_Runs_Combined_BvC_genes, file = paste("Foley_Runs_Combined_AvB_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))


#### NetAct####

#gsearslt_1 = TF_GSEA(mgs, DErslt$`B-C`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_1, file = paste("Foley_Runs_Combined_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = TF_GSEA(mgs, DErslt$`A-B`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_2, file = paste("Foley_Runs_Combined_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = TF_GSEA(mgs, DErslt$`A-C`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_3, file = paste("Foley_Runs_Combined_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#gsearslt_1 = read.csv(paste("Foley_Runs_Combined_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = read.csv(paste("Foley_Runs_Combined_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = read.csv(paste("Foley_Runs_Combined_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#tfs_1 = as.character(gsearslt_1$tf[gsearslt_1$qvals <0.2])
#tfs_2 = as.character(gsearslt_2$tf[gsearslt_2$qvals <0.2])
#tfs_3 = as.character(gsearslt_3$tf[gsearslt_3$qvals <0.2])
#tfs = sort(unique(as.character(c(tfs_1, tfs_2, tfs_3))))
#lengths(list(tfs_1,tfs_2,tfs_3, tfs))
#tfs


#acts_mat = TF_Activity(tfs, mgs, neweset, DErslt$Overall)$all_activities #For all three comparisons
#acts_mat_test = TF_Activity(tfs, mgs, neweset, DErslt$`B-C`)$all_activities #Only for BvC comparison




#Combine_heatmap(acts_mat, neweset) #Heatmap

#tf_links = TF_Filter(acts_mat, mgs, miTh = 0.2, nbins = 8, method = "spearman")

#x = union(tf_links$Source,tf_links$Target)

#dim(tf_links)
#plot_network_v(tf_links)

#write.table(tf_links, "net_emt.tpo", quote = FALSE, row.names = FALSE)
#write.table(tf_links, paste("Foley_Runs_Combined_net_emt_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
#           quote = FALSE, row.names = FALSE)

####PERTURBATION
#library(sRACIPE)
#racipe <- sRACIPE::simulateGRC(paste("Foley_Runs_Combined_net_emt_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
#                              quote = FALSE, row.names = FALSE, plots = TRUE, plotToFile = TRUE)

#racipePer <- sRACIPE::knockdownAnalysis(racipe,reduceProduction = 10)

