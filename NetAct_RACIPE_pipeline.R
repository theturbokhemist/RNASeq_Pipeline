# read counts data
setwd("/Users/gordid/Desktop/Lucas")

SkipCalcs = TRUE

phenoData_lucas = new("AnnotatedDataFrame", data = data.frame(celltype_lucas = group_name_lucas))
rownames(phenoData_lucas) = sample_names_lucas_final
head(Combined_data_lucas_NetAct$counts)
phenoData_lucas
colnames(Combined_data_lucas_NetAct$counts) <- sample_names_lucas_final
head(Combined_data_lucas_NetAct$counts)

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
for (i in 1:(length(levels(group_name_lucas))-1)) {
  for (j in (i+1):length(levels(group_name_lucas))) {
    compare_list <- append(compare_list, paste0(levels(group_name_lucas)[i], "-", levels(group_name_lucas)[j]))
    
  }
}
compare_list

if (TwoGroups == TRUE) {
  DErslt = RNAseqDegs_limma(Combined_data_lucas_NetAct$counts, phenoData_lucas, compare_list)
  str(DErslt)
  
  e = DErslt$e
  
  neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_lucas)
  
  save(neweset, DErslt, file = paste(Directory_CSVs_NetAct_RACIPE, "/macro_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".RData", sep = ""))
  
  write.csv(DErslt$table, file = paste(Directory_CSVs_NetAct_RACIPE, "/DEGenes_Lucas_RACIPE", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                       filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
  
  if (SkipCalcs == FALSE) {
    gsearslt_lucas <- TF_GSEA(specimen, DErslt, minSize = 5, nperm = 10000, qval = T)
    write.csv(gsearslt_lucas, file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
  }
  setwd(Directory_CSVs_NetAct_RACIPE)
  gsearslt_lucas <- read.csv(file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                          filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
  qval <- 0.2
  Tfs_lucas <- as.character(gsearslt_lucas$tf[gsearslt_lucas$qvals < qval])
  tfs <- sort(unique(as.character(c(Tfs_lucas))))
  acts_mat = TF_Activity(tfs, specimen, neweset, DErslt)$all_activities
  
  
  } else {
  DErslt = MultiRNAseqDegs_limma(Combined_data_lucas_NetAct$counts, phenoData_lucas, compare_list)
  str(DErslt)
  e = DErslt$Overall$e
  
  neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_lucas)
  
  save(neweset, DErslt, file = paste(Directory_CSVs_NetAct_RACIPE, "/macro_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".RData", sep = ""))
  
  NetAct_gene_list_lucas <- list()
  gsearslt_list <- list()
  tfs_list <- list()
  TFs_vector <- c()
  for (i in 1:length(compare_list)) {
  NetAct_gene_list_lucas[[i]] <- DErslt[[i]]$table
  write.csv(NetAct_gene_list_lucas[[i]], file = paste(Directory_CSVs_NetAct_RACIPE, "/DEGenes_Lucas_RACIPE", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                      filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
  
  if (SkipCalcs == FALSE) {
  gsearslt_lucas_multi <- TF_GSEA(specimen, DErslt[[i]], minSize = 5, nperm = 10000, qval = T)
  write.csv(gsearslt_lucas_multi, file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Lucas", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
  }
  setwd(Directory_CSVs_NetAct_RACIPE)
  gsearslt_list[[i]] <- read.csv(paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Lucas", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                       filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

  qval <- 0.2
  tfs_list[[i]] <- as.character(gsearslt_list[[i]]$tf[gsearslt_list[[i]]$qvals < qval])
  TFs_vector <- append(TFs_vector, tfs_list[[i]])
}
names(NetAct_gene_list_lucas) <- Comparisons_vector
names(gsearslt_list) <- Comparisons_vector
names(tfs_list) <- Comparisons_vector

tfs <- sort(unique(as.character(TFs_vector)))
lengths(tfs_list)
length(tfs)
tfs
acts_mat = TF_Activity(tfs, specimen, neweset, DErslt$Overall)$all_activities
}
setwd("/Users/gordid/Desktop/Lucas")

if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_NetAct_RACIPE, "/TFsHeatMap_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), "_qval", 
                       gsub(".", "", as.character(qval), fixed = TRUE), filter_type, gsub(".*filter", "", Directory_Plots_NetAct_RACIPE),
                       if (outliers == TRUE) {print("_outliers-removed")}, ".pdf", sep = ""))}

Combine_heatmap(acts_mat, neweset) #Heatmap

if (SavePlots == TRUE) {dev.off()}

tf_links <- TF_Filter(acts_mat, specimen, miTh = 0.05, nbins = 8, method = "spearman")
miTh = 0.05
dim(tf_links)

####Generating Network####
if (SavePlots == TRUE) {
setwd(Directory_Plots_NetAct_RACIPE)
  plot_network_v_html(tf_links)
webshot::webshot(paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Lucas_html", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                       "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
                       if (outliers == TRUE) {print("_outliers-removed")}, ".html", sep = ""),
                 file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                       "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
                              if (outliers == TRUE) {print("_outliers-removed")}, ".png", sep = ""), delay = 2)
setwd("/Users/gordid/Desktop/Lucas")
plot_network_v(tf_links)
} else {
  plot_network_v(tf_links)
}

x <- union(tf_links$Source,tf_links$Target)



#if (SavePlots == TRUE) {
 # orca(plot_network_v(tf_links), file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
  #    "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
   #                              if (outliers == TRUE) {print("_outliers-removed")}, ".png", sep = "")) 
#} else {
#  plot_network_v(tf_links)
#}

write.csv(tf_links, file = paste(Directory_CSVs_NetAct_RACIPE, "/TFNetwork_InteractionsTable_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                          "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
                          if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#write.table(tf_links, "net_emt.tpo", quote = FALSE, row.names = FALSE)
write.table(tf_links, file = paste(Directory_Plots_NetAct_RACIPE, "/net_emt_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                      "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
                                   if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""), quote = FALSE, row.names = FALSE)

####PERTURBATION####
if (SavePlots == TRUE) {
setwd(Directory_Plots_NetAct_RACIPE)
library(sRACIPE)

racipe <- sRACIPE::simulateGRC(paste(Directory_Plots_NetAct_RACIPE, "/net_emt_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                  "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
                                     if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
                                     quote = FALSE, row.names = FALSE, plots = TRUE, plotToFile = TRUE)
racipePer <- sRACIPE::knockdownAnalysis(racipe,reduceProduction = 10)
setwd("/Users/gordid/Desktop/Lucas")
} else {
  print("No Perturbation")
}
#gsearslt_1 = read.csv(paste("Lucas_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = read.csv(paste("Lucas_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = read.csv(paste("Lucas_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#tfs_1 = as.character(gsearslt_1$tf[gsearslt_1$qvals <0.2])
#tfs_2 = as.character(gsearslt_2$tf[gsearslt_2$qvals <0.2])
#tfs_3 = as.character(gsearslt_3$tf[gsearslt_3$qvals <0.2])
#tfs = sort(unique(as.character(c(tfs_1, tfs_2, tfs_3))))
#lengths(list(tfs_1,tfs_2,tfs_3, tfs))
#tfs
# for normal use ####
#e = DErslt$Overall$e
#neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_lucas)
#save(neweset, DErslt, file = "macro.RData")

#Lupus_Lucas_BvC_genes <- (DErslt$`B-C`$table)
#head(Lupus_Lucas_BvC_genes)

#write.csv(Lupus_Lucas_BvC_genes, file = paste("Lucas_BvC_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#Lupus_Lucas_AvC_genes <- (DErslt$`A-C`$table)
#write.csv(Lupus_Lucas_BvC_genes, file = paste("Lucas_AvC_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#Lupus_Lucas_AvB_genes <- (DErslt$`A-B`$table)
#write.csv(Lupus_Lucas_BvC_genes, file = paste("Lucas_AvB_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))


#### NetAct####

#gsearslt_1 = TF_GSEA(mgs, DErslt$`B-C`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_1, file = paste("Lucas_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = TF_GSEA(mgs, DErslt$`A-B`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_2, file = paste("Lucas_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = TF_GSEA(mgs, DErslt$`A-C`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_3, file = paste("Lucas_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#gsearslt_1 = read.csv(paste("Lucas_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = read.csv(paste("Lucas_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = read.csv(paste("Lucas_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
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
#write.table(tf_links, paste("Lucas_net_emt_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
 #           quote = FALSE, row.names = FALSE)

####PERTURBATION
#library(sRACIPE)
#racipe <- sRACIPE::simulateGRC(paste("Lucas_net_emt_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
 #                              quote = FALSE, row.names = FALSE, plots = TRUE, plotToFile = TRUE)

#racipePer <- sRACIPE::knockdownAnalysis(racipe,reduceProduction = 10)

