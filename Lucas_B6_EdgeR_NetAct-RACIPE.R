#read in data file on experiment
setwd("/Users/gordid/Desktop/Lucas_B6/")

#Rename files
Rename_file_extensions <- function(directory) {
  files_list <<- list.files(path = directory, pattern = "*.results")
  
  if (length(which(grepl(glob2rx("*.results"), list.files(directory)))) > 0) {
    files_list_renamed <<- gsub(".results", ".txt", files_list)
    setwd(directory)
    for (i in 1:length(files_list)) {
      file.rename(from = file.path(directory, files_list[i]), to = file.path(directory ,files_list_renamed[i]))
    }
  } else {
    print("No '.results' file found")                                      
  }
}
Rename_file_extensions("/Users/gordid/Desktop/Lucas_B6/genes_results/")


#reading in design CSV
experiment_Lucas_B6 <- read.csv(paste0("/Users/gordid/Desktop/Lucas_B6/", "Lucas_B6_experimental_design", ".csv"))
experiment_Lucas_B6

#creating vector of data files
files_Lucas_B6_unsorted <- list.files("/Users/gordid/Desktop/Lucas_B6/genes_results/")
files_Lucas_B6_unsorted

#Putting data files character vector into the order of the sample IDS in the design experiment CSV
files_Lucas_B6 <- c()

for (i in 1:length(files_Lucas_B6_unsorted)) {
  for (j in 1:length(experiment_Lucas_B6$Sample.Name)) {
    if (grepl(as.character(experiment_Lucas_B6$Sample.Name)[j], files_Lucas_B6_unsorted[i]) == TRUE) {
      files_Lucas_B6[j] <- files_Lucas_B6_unsorted[i]
    }
  }
}
files_Lucas_B6


####Save plots?####
SavePlots = TRUE
SaveCSVs = TRUE
EnrichR = TRUE
FGSEA = TRUE

####What Organism?####
Mouse = TRUE 

####Comparing just 2 groups?####
TwoGroups = TRUE

experiment_Lucas_B6$Group.name

if (TwoGroups == TRUE) {
  
  Group1 = "T4wS"
  Group2 = "T96S"
  

  WhichGroups = c(which(experiment_Lucas_B6$Group.Name == Group1),which(experiment_Lucas_B6$Group.Name == Group2))
  
  experiment_Lucas_B6 <- experiment_Lucas_B6[WhichGroups,]
  
  files_Lucas_B6 <- files_Lucas_B6[WhichGroups]
}


####Are there Outliers?####
outliers = TRUE

if (outliers == TRUE) {
  Outlier1 = "A091"
  Outlier2 = ""
  
  outliers_index <- c(which(experiment_Lucas_B6$Sample.Name == Outlier1),which(experiment_Lucas_B6$Sample.Name == Outlier2))
  
  outliers_name <- as.character(experiment_Lucas_B6$Sample.Name[outliers_index])
  
  experiment_Lucas_B6 <- experiment_Lucas_B6[-outliers_index,]
  
  files_Lucas_B6 <- files_Lucas_B6[-outliers_index]
  
}
#Better if you could put in comparisons you are interested in

####Are you using CPM or counts to filter?####
CPMFilter = TRUE

if (CPMFilter == TRUE) {
  filter_type = "_cpm-filter"
} else {
  filter_type = "_counts-filter"
}

#reading in one of the data files
setwd("/Users/gordid/Desktop/Lucas_B6/genes_results/")
file_1_Lucas_B6 <- read.delim(files_Lucas_B6[1], header=TRUE)
head(file_1_Lucas_B6)
str(file_1_Lucas_B6)
dim(file_1_Lucas_B6)
setwd("/Users/gordid/Desktop/Lucas_B6")

####Renaming columns and creating groups####

group_name_Lucas_B6 <- as.character(experiment_Lucas_B6$Group.Name)
group_name_Lucas_B6 <- factor(as.factor(group_name_Lucas_B6), levels = unique(as.factor(group_name_Lucas_B6)) ) #preserves order of factors
group_name_Lucas_B6
levels(group_name_Lucas_B6)


sample_name_Lucas_B6 <- as.character(experiment_Lucas_B6$Sample.Name)
sample_name_Lucas_B6

Age.at.Harvest_Lucas_B6 <- factor(as.factor(as.character(experiment_Lucas_B6$Age.at.Harvest)), levels = unique(as.factor(as.character(experiment_Lucas_B6$Age.at.Harvest))) )
Age.at.Harvest_Lucas_B6

Time.on.Treatment_Lucas_B6 <- factor(as.factor(as.character(experiment_Lucas_B6$Time.on.Treatment)), levels = unique(as.factor(as.character(experiment_Lucas_B6$Time.on.Treatment))) )
Time.on.Treatment_Lucas_B6

Lane_Lucas_B6 <- factor(as.factor(as.character(experiment_Lucas_B6$Lane)), levels = unique(as.factor(as.character(experiment_Lucas_B6$Lane))) )
Lane_Lucas_B6

Remarks_Lucas_B6 <- factor(as.factor(as.character(experiment_Lucas_B6$Remarks)), levels = unique(as.factor(as.character(experiment_Lucas_B6$Remarks))) )
Remarks_Lucas_B6

####creating vector of gene names and vector of possible sample comparisons####

Comparisons_vector <- c()
for (i in 1:(length(levels(group_name_Lucas_B6))-1)) {
  for (j in (i+1):length(levels(group_name_Lucas_B6))) {
    
    Comparisons_vector <- append(Comparisons_vector, paste0(levels(group_name_Lucas_B6)[i], "vs", levels(group_name_Lucas_B6)[j]))
    
  }
}

Comparisons_vector
length(Comparisons_vector)

#Creating new directorys####

#plots/RNASeq folder
Directory_Plots_RNASeq <- paste0("/Users/gordid/Desktop/Lucas_B6/All_results/Plots/RNASeq/",
                                 (if(TwoGroups == TRUE) { 
                                   paste0(levels(group_name_Lucas_B6)[1], "vs", levels(group_name_Lucas_B6)[2],
                                          if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                 } else {
                                   paste0("all-data", if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                 }))

if (dir.exists(Directory_Plots_RNASeq) == FALSE) {
  dir.create(Directory_Plots_RNASeq, recursive = TRUE)
} else {
  "Already made this directory!"
}

#plots/NetAct_RACIPE folder

Directory_Plots_NetAct_RACIPE <- paste0("/Users/gordid/Desktop/Lucas_B6/All_results/Plots/NetAct_RACIPE/",
                                        (if(TwoGroups == TRUE) { 
                                          paste0(levels(group_name_Lucas_B6)[1], "vs", levels(group_name_Lucas_B6)[2],
                                                 if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                        } else {
                                          paste0("all-data", if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                        }))

if (dir.exists(Directory_Plots_NetAct_RACIPE) == FALSE) {
  dir.create(Directory_Plots_NetAct_RACIPE, recursive = TRUE)
} else {
  "Already made this directory!"
}


#CSVs/RNASeq folder
Directory_CSVs_DEGenes <- paste0("/Users/gordid/Desktop/Lucas_B6/All_results/CSVs/DEGenes/", 
                                 (if(TwoGroups == TRUE) {
                                   paste0(levels(group_name_Lucas_B6)[1], "vs", levels(group_name_Lucas_B6)[2],
                                          if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                 } else {
                                   paste0("all-data", if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                 }))

if (dir.exists(Directory_CSVs_DEGenes) == FALSE) {
  dir.create(Directory_CSVs_DEGenes, recursive = TRUE)
} else {
  "Already made this directory!"
}

#CSVs/NetAct_RACIPE folder
Directory_CSVs_NetAct_RACIPE <- paste0("/Users/gordid/Desktop/Lucas_B6/All_results/CSVs/NetAct_RACIPE/",
                                       (if(TwoGroups == TRUE) { 
                                         paste0(levels(group_name_Lucas_B6)[1], "vs", levels(group_name_Lucas_B6)[2],
                                                if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                       } else {
                                         paste0("all-data", if(outliers == TRUE) {paste0("_", "outlier-", paste(outliers_name, collapse = "_"))}, filter_type)
                                       }))


if (dir.exists(Directory_CSVs_NetAct_RACIPE) == FALSE) {
  dir.create(Directory_CSVs_NetAct_RACIPE, recursive = TRUE)
} else {
  "Already made this directory!"
}

####CREATE DGEList OBJECT####
setwd("/Users/gordid/Desktop/Lucas_B6/genes_results/")
combined_data_Lucas_B6 <- readDGE(files_Lucas_B6, columns = c(1, 5))
dim(combined_data_Lucas_B6) #showing dimensions of "counts" matrix
class(combined_data_Lucas_B6) #showing class of the whole object
summary(combined_data_Lucas_B6) #showing names, size, class, and mode of each element of the object, not sure why class is "none" for counts
str(combined_data_Lucas_B6) #showing the structure of the object
head(combined_data_Lucas_B6$counts)
colnames(combined_data_Lucas_B6) <- sample_name_Lucas_B6
setwd("/Users/gordid/Desktop/Lucas_B6")

####Assigning experimental info to categories in the DGEList obect####
combined_data_Lucas_B6$samples$group <- group_name_Lucas_B6
combined_data_Lucas_B6$samples

#combined_data_Lucas_B6$samples$remarks <- remarks_Lucas_B6

####Create new genes matrix by retrieving gene information from the mus library. The gene information of their SYMBOL and TXCHROM#### 
#number of gene IDs
gene_id_Lucas_B6 <- rownames(combined_data_Lucas_B6)
head(gene_id_Lucas_B6)
length(gene_id_Lucas_B6)
#number of duplicated ENSEMBL IDS PRE select
dup_genes_Lucas_B6_PRESelect <- which(duplicated(gene_id_Lucas_B6))
length(dup_genes_Lucas_B6_PRESelect)

columns(Mus.musculus)
if (Mouse == TRUE) {
  genes_Lucas_B6 <- biomaRt::select(Mus.musculus, keys = gene_id_Lucas_B6, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
} else {
  genes_Lucas_B6 <- biomaRt::select(Homo.sapiens, keys = gene_id_Lucas_B6, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
  
}
dim(genes_Lucas_B6)
head(genes_Lucas_B6)
length(genes_Lucas_B6$ENSEMBL) - length(gene_id_Lucas_B6)

####symbols are duplicated too often and so I dont remove genes based on duplicated symbols (cause of NA mostly, will remove NAs later)####
head((which(duplicated(genes_Lucas_B6$SYMBOL))))
length(which(duplicated(genes_Lucas_B6$SYMBOL)))
genes_Lucas_B6[head((which(duplicated(genes_Lucas_B6$SYMBOL)))),]

#number of duplicated ENSEMBL IDS POST select
dup_genes_Lucas_B6 <- which(duplicated(genes_Lucas_B6$ENSEMBL))
length(dup_genes_Lucas_B6)
head(dup_genes_Lucas_B6)
genes_Lucas_B6[c(head(which(duplicated(genes_Lucas_B6$ENSEMBL))), head(which(duplicated(genes_Lucas_B6$ENSEMBL))) - 1),]

####remove duplicated genes based on duplicated ENSEMBL ID####
genes_Lucas_B6_final <- genes_Lucas_B6[-dup_genes_Lucas_B6, ]
dim(genes_Lucas_B6_final)
dim(combined_data_Lucas_B6$counts)

#Number of removed Symbol IDS must be equal to or less than the number of removed ENSEMBL IDs because during the select process, some ENSEMBL ID's may match to two different Symbol IDs.
length(which(duplicated(genes_Lucas_B6$SYMBOL))) - length(which(duplicated(genes_Lucas_B6_final$SYMBOL)))
length(dup_genes_Lucas_B6)

#assign data frame of gene annotations to the combined_data object
combined_data_Lucas_B6$genes <- genes_Lucas_B6_final 
summary(combined_data_Lucas_B6)

#RAW COUNTS DATA
if (SaveCSVs == TRUE) {
  
  write.csv(cbind2(combined_data_Lucas_B6$counts, combined_data_Lucas_B6$genes), file = paste0(Directory_CSVs_DEGenes,
                                                                                                   "/Raw_Counts_Data_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), ".csv"))
  
}

####CPM normalization on count data####
cpm_Lucas_B6 <- cpm(combined_data_Lucas_B6)
lcpm_Lucas_B6 <-cpm(combined_data_Lucas_B6, log = TRUE)

L_Lucas_B6 <- mean(combined_data_Lucas_B6$samples$lib.size) * 1e-6 #calculation of L for the prior count.
M_Lucas_B6 <- median(combined_data_Lucas_B6$samples$lib.size) * 1e-6
c(L_Lucas_B6, M_Lucas_B6)
head(lcpm_Lucas_B6)

table(rowSums(combined_data_Lucas_B6$counts==0)==length(colnames(lcpm_Lucas_B6)))

###FILTER THE GENES####
if (CPMFilter == TRUE) {
  kept_expression_Lucas_B6 <- filterByExpr(combined_data_Lucas_B6, group = combined_data_Lucas_B6$samples$group) #filterByExpr function
  # keeps genes with 10 read counts or more in a minimum number of samples, where the minimum number of samples is chosen according 
  #to the minimum group sample size. 
}  else {
  #filter by counts, rather than by cpm
  kept_expression_Lucas_B6 <- (rowSums(combined_data_Lucas_B6$counts > 10) >= 1)
}
length(kept_expression_Lucas_B6)
head(kept_expression_Lucas_B6)

#removing all the lowly expressed genes from the combined_data object. When you subset a DGEList and specify keep.lib.sizes=FALSE,
#the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
filtered_combined_data_Lucas_B6 <- combined_data_Lucas_B6[kept_expression_Lucas_B6, , keep.lib.sizes = FALSE] 
dim(filtered_combined_data_Lucas_B6)
dim(combined_data_Lucas_B6)

####Removing duplicated gene symbols for the NetAct analysis####
Combined_data_Lucas_B6_NetAct <- filtered_combined_data_Lucas_B6
rownames(Combined_data_Lucas_B6_NetAct$counts) <- Combined_data_Lucas_B6_NetAct$genes$SYMBOL
dim(Combined_data_Lucas_B6_NetAct)
#Number of Symbols that are NA
dup_genes_Lucas_B6_NA <- which(is.na(filtered_combined_data_Lucas_B6$genes$SYMBOL))
length(dup_genes_Lucas_B6_NA)
#Number of duplicated Symbols that are not NA
length(which(duplicated(na.omit(filtered_combined_data_Lucas_B6$genes$SYMBOL))))
#Number of duplicated symbols + all NAs
length(dup_genes_Lucas_B6_NA) + length(which(duplicated(na.omit(filtered_combined_data_Lucas_B6$genes$SYMBOL))))
# of genes theoretically left over at removing of all duplicated
length(filtered_combined_data_Lucas_B6$genes$SYMBOL) -(length(dup_genes_Lucas_B6_NA) + length(which(duplicated(na.omit(filtered_combined_data_Lucas_B6$genes$SYMBOL)))))

#removing all NAs
Combined_data_Lucas_B6_NetAct <-Combined_data_Lucas_B6_NetAct[-dup_genes_Lucas_B6_NA, ]
dim(Combined_data_Lucas_B6_NetAct)

#removing all duplicated symbols
dup_genes_Lucas_B6_SYMBOL <- which(duplicated(Combined_data_Lucas_B6_NetAct$genes$SYMBOL))
length(dup_genes_Lucas_B6_SYMBOL)
Combined_data_Lucas_B6_NetAct <- Combined_data_Lucas_B6_NetAct[-dup_genes_Lucas_B6_SYMBOL, ]
dim(Combined_data_Lucas_B6_NetAct)

####Plotting filtered Data####
lcpm.cutoff_Lucas_B6 <- log2(10/M_Lucas_B6 + 2/L_Lucas_B6)
lcpm.cutoff_Lucas_B6
nsamples <- ncol(combined_data_Lucas_B6$counts)
col <- brewer.pal(nsamples, "Paired")

if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/FilteredData_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                  ".pdf", sep = ""))}
par(mfrow = c(1,2))
plot(density(lcpm_Lucas_B6[, 1]), col = col[1], lwd=2, ylim=c(0,0.46), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_Lucas_B6, lty=3)
title(main="A. Raw data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_Lucas_B6[, i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

combined_data_Lucas_B6 <- filtered_combined_data_Lucas_B6
lcpm_Lucas_B6 <- cpm(combined_data_Lucas_B6, log = TRUE)

plot(density(lcpm_Lucas_B6[, 1]), col = col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_Lucas_B6, lty=3)
title(main="B. Filtered data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_Lucas_B6[, i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
if (SavePlots == TRUE) {dev.off() }


####Calculate normalization factors and creates new DGEList object with these new normalization factors inserted. THEN GENERATING BOX PLOT####
colnames(lcpm_Lucas_B6) <- sample_name_Lucas_B6
combined_data_Lucas_B6$samples$norm.factors

if (SavePlots == TRUE) {pdf(paste0(Directory_Plots_RNASeq, "/TMM_Normalized_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), ".pdf"))}


par(mfrow = c(1,2))
boxplot(lcpm_Lucas_B6, las = 2, col = col, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")

combined_data_Lucas_B6 = calcNormFactors(combined_data_Lucas_B6, method = "TMM") #calculates normalization factors and creates new DGEList object
#these normalization factors are important because if there is a sample where expression level is different accross the board, that may have been caused
#by an eternal factor that is not of biological interest. TMM helps normalize for that, using the assumption that should have a similar range
#and distribution of expression values.
lcpm_Lucas_B6 <- cpm(combined_data_Lucas_B6, log = TRUE) 
colnames(lcpm_Lucas_B6) <- sample_name_Lucas_B6
head(lcpm_Lucas_B6)

combined_data_Lucas_B6$samples$norm.factors
boxplot(lcpm_Lucas_B6, las = 2, col = col, main="")
title(main="B.Normalised with TMM data", ylab="Log-cpm")
if (SavePlots == TRUE) {dev.off() }


combined_data_Lucas_B6$samples$names <- sample_name_Lucas_B6
combined_data_Lucas_B6$samples

combined_data_Lucas_B6$samples$lanes <- Lane_Lucas_B6

####PCA and PCoA Plots####
if (SavePlots == TRUE) {pdf(paste0(Directory_Plots_RNASeq, "/PCA_bygroup_1vs2_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                   ".pdf"))}
color_groups_Lucas_B6 <- combined_data_Lucas_B6$samples$group
levels(color_groups_Lucas_B6) <- brewer.pal(nlevels(color_groups_Lucas_B6), "Set1")
color_groups_Lucas_B6 <- as.character(color_groups_Lucas_B6)
color_groups_Lucas_B6

color_lanes_Lucas_B6 <- combined_data_Lucas_B6$samples$lanes; levels(color_lanes_Lucas_B6) <- brewer.pal(nlevels(color_lanes_Lucas_B6), "Set1"); color_lanes_Lucas_B6 <- as.character(color_lanes_Lucas_B6)

par(mfrow = c(1,1))
#plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$names, col = color_groups_Lucas_B6, main = "PCoA")
plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$group, col = color_groups_Lucas_B6, gene.selection="common", 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))
if (SavePlots == TRUE) {dev.off() }

par(mfrow = c(1,2))
plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$names, col = color_groups_Lucas_B6, dim = c(2,3), 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCoA")) 
plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$names, col = color_groups_Lucas_B6, gene.selection="common", dim = c(2,3),
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))

par(mfrow = c(1,1))
plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$lanes, col = color_lanes_Lucas_B6, gene.selection="common", 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))


glMDSPlot(lcpm_Lucas_B6, labels=paste(combined_data_Lucas_B6$samples$names, combined_data_Lucas_B6$samples$group, sep="_"), 
          groups = combined_data_Lucas_B6$samples[,c(2,5)], 
          gene.selection="common", main = paste0(sub("_[^_]+$*" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"), launch=TRUE, path = Directory_Plots_RNASeq,
          html = paste0("InteractivePCA_ScreePlot_1vs2_Lucas_B6", "_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")))

#scree.plot(plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$names, col = color_groups_Lucas_B6, gene.selection="common", main = "PCA"))
#a <- (plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$names, col = color_groups_Lucas_B6, main = "PCoA"))
#screeplot(a, xlim = c(-8,13), ylim=c(0,0.46))
####PCA Plots for other phenotypic data and meta-data####
#color_age_Lucas_B6 <- combined_data_Lucas_B6$samples$age
#levels(color_age_Lucas_B6) <- brewer.pal(nlevels(color_age_Lucas_B6), "Set2")
#color_age_Lucas_B6 <- as.character(color_age_Lucas_B6)
#color_age_Lucas_B6
#plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$age, col = color_groups_Lucas_B6)
#plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$age, col = color_groups_Lucas_B6, gene.selection="common")

#color_x2dg_Lucas_B6 <- combined_data_Lucas_B6$samples$x2dg
#levels(color_x2dg_Lucas_B6) <- brewer.pal(nlevels(color_x2dg_Lucas_B6), "Set3")
#color_x2dg_Lucas_B6 <- as.character(color_x2dg_Lucas_B6)
#color_x2dg_Lucas_B6
#plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$x2dg, col = color_groups_Lucas_B6)
#plotMDS(lcpm_Lucas_B6, labels = combined_data_Lucas_B6$samples$x2dg, col = color_groups_Lucas_B6, gene.selection="common")

####Building Design Matrix####

design_matrix_Lucas_B6_no.int <- model.matrix(~0+combined_data_Lucas_B6$samples$group) #no intercept
colnames(design_matrix_Lucas_B6_no.int) = gsub("combined_data_Lucas_B6$samples$group", "", colnames(design_matrix_Lucas_B6_no.int), fixed = TRUE)
design_matrix_Lucas_B6_no.int

design_matrix_Lucas_B6_int <- model.matrix(~combined_data_Lucas_B6$samples$group) #with intercept
colnames(design_matrix_Lucas_B6_int) = gsub("combined_data_Lucas_B6$samples$group", "", colnames(design_matrix_Lucas_B6_int), fixed = TRUE)
design_matrix_Lucas_B6_int


if (TwoGroups == TRUE) {
  contr.matrix_no.int <- makeContrasts(
    paste(levels(group_name_Lucas_B6)[1]," - ", levels(group_name_Lucas_B6)[2], sep = ""),
    levels = colnames(design_matrix_Lucas_B6_no.int))
  colnames(contr.matrix_no.int) <- paste(levels(group_name_Lucas_B6)[1],"vs", levels(group_name_Lucas_B6)[2], sep = "")
  contr.matrix_no.int
} else {
  contrast_vector <- c()
  column_name_vector <- c()
  for (i in 1:(length(levels(group_name_Lucas_B6))-1)) {
    
    for (j in (i+1):length(levels(group_name_Lucas_B6))) {
      
      contrast_vector <- append(contrast_vector, paste(levels(group_name_Lucas_B6)[i]," - ", levels(group_name_Lucas_B6)[j], sep = ""))
      column_name_vector <- append(column_name_vector, paste(levels(group_name_Lucas_B6)[i], "vs", levels(group_name_Lucas_B6)[j], sep = ""))
    }
    
  }
  
  contr.matrix_no.int <- makeContrasts(
    contrasts = contrast_vector,
    levels = colnames(design_matrix_Lucas_B6_no.int)
  )
  colnames(contr.matrix_no.int) <- column_name_vector
  contr.matrix_no.int
}


#contr.matrix_no.int <- makeContrasts(
#AvsB = "A- B",
#AvsC = A - C,
#BvsC = B - C,
#levels = colnames(design_matrix_Lucas_B6_no.int))
#contr.matrix_no.int
#colnames(contr.matrix_no.int) <- c("a", "b", "c")

#contr.matrix_no.int <- makeContrasts(
#contrasts = contrast_vector,
#levels = colnames(design_matrix_Lucas_B6_no.int))
#contr.matrix_no.int

#contr.matrix_int <- makeContrasts( #Intercept column has replcaced "A", not sure how contrast matrix works now
# AvsB = A - B,
# AvsC = A - C,
# BvsC = B - C,
#levels = colnames(design_matrix_Lucas_B6_int))
#contr.matrix_int

####VOOM####
adjP = 0.05
if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Voom_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                  ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,2))
v_Lucas_B6 <- voom(combined_data_Lucas_B6, design_matrix_Lucas_B6_no.int, plot=TRUE)
str(v_Lucas_B6)
dim(v_Lucas_B6$genes)

vfit_Lucas_B6 <- lmFit(v_Lucas_B6, design_matrix_Lucas_B6_no.int)
summary(decideTests(vfit_Lucas_B6))
str(vfit_Lucas_B6)

vfit_Lucas_B6_contrasts <- contrasts.fit(vfit_Lucas_B6, contrasts=contr.matrix_no.int) 
summary(decideTests(vfit_Lucas_B6_contrasts))

efit_Lucas_B6 <- eBayes(vfit_Lucas_B6_contrasts)
summary(decideTests(efit_Lucas_B6))
de_Lucas_B6 <- decideTests(efit_Lucas_B6)

plotSA(efit_Lucas_B6)
title(main="voom: Mean-variance trend-Normalized")

if (SavePlots == TRUE) {dev.off()}


####Filter further with logfold change and Venn Diagram####
#Figure out the criteria data has to choose appropriate lfc
lfc = 0.1
tfit_Lucas_B6 <- treat(vfit_Lucas_B6_contrasts, lfc=lfc)
dt_Lucas_B6 <- decideTests(tfit_Lucas_B6)
summary(dt_Lucas_B6)


topTable_p <-  topTable(efit_Lucas_B6, n=Inf); head(topTable_p)
topTreat_p <- topTreat(efit_Lucas_B6, n=Inf); head(topTreat_p)
topTreat_lfc <- topTreat(tfit_Lucas_B6, n=Inf); head(topTreat_lfc)

dim(topTable_p)
dim(topTreat_p)
dim(topTreat_lfc)

sum(topTable_p$P.Val < 0.05); sum(topTreat_p$adj.P.Val < 0.05)
sum(topTreat_lfc$P.Value < 0.05); sum(topTreat_lfc$adj.P.Val < 0.05)

topTable_p_symbols_1 <- topTable_p[!is.na(topTable_p$SYMBOL), ]; dim(topTable_p_symbols_1)
topTable_p_symbols <- topTable_p_symbols_1[-which(duplicated(topTable_p_symbols_1$SYMBOL)),]; dim(topTable_p_symbols)
sum(topTable_p_symbols$adj.P.Val < 0.05)

#which genes are DE in multiple comparisons    
#Figure out how to add this step to the pipeline, ALSO why do I need to add print function in order for them to show up
if (TwoGroups == FALSE) {
  de.common_Lucas_B6_2_comparisons <- which(dt_Lucas_B6[,1]!=0 & dt_Lucas_B6[,2]!=0)
  print(length(de.common_Lucas_B6_2_comparisons))
  
  de.common_Lucas_B6_3_comparisons <- which(dt_Lucas_B6[,1]!=0 & dt_Lucas_B6[,2]!=0 & dt_Lucas_B6[,3]!=0)
  print(length(de.common_Lucas_B6_3_comparisons))
  
  if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Venn_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                    ".pdf", sep = ""), width=13, height=10)}
  par(mfrow = c(1,1))
  vennDiagram(dt_Lucas_B6[,1:3], circle.col=c("turquoise", "salmon", "orange"))
  if (SavePlots == TRUE) {dev.off()}
}

if (SaveCSVs == TRUE) {
  write.fit(tfit_Lucas_B6, dt_Lucas_B6, file = paste(Directory_CSVs_DEGenes, "/tfit_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                         "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".txt", sep = ""))
  
  tfit_Lucas_B6.txt = read.table(file = paste(Directory_CSVs_DEGenes, "/tfit_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                                gsub(".", "", as.character(lfc), fixed = TRUE),
                                                ".txt", sep = ""))
  
  write.csv(tfit_Lucas_B6.txt, file = paste(Directory_CSVs_DEGenes, "/tfit_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                              gsub(".", "", as.character(lfc), fixed = TRUE),
                                              ".csv", sep = ""))
}

####Mean-difference plots and arranging DE genes####
if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/Mean-Difference-Plots_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                         gsub(".", "", as.character(lfc), fixed = TRUE),
                                         ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,length(Comparisons_vector)))
for (i in 1:length(Comparisons_vector)) {
  plotMD(efit_Lucas_B6, column = i, status = de_Lucas_B6[,i], xlim = c(-8,13), main = paste0(Comparisons_vector[i], " Mean-Difference Plot"))
  
  glMDPlot(efit_Lucas_B6, coef=i, status=de_Lucas_B6, main=colnames(efit_Lucas_B6)[i], side.main="ENSEMBL", counts=lcpm_Lucas_B6,
           groups=combined_data_Lucas_B6$samples$group,launch=FALSE, 
           path = Directory_Plots_RNASeq, html = paste0("Interactive_MeanDifference_Plot_Lucas_B6", "_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), 
                                                        if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE)))
  
}
if (SavePlots == TRUE){ dev.off()}

#par(mfrow = c(1,length(levels(group_name_Lucas_B6))))
#Ranked_DEGenes_list <- list()
#for (i in 1:(length(levels(group_name_Lucas_B6))-1)) {
#  for (j in (i+1):length(levels(group_name_Lucas_B6))) {
#    print(paste0(levels(group_name_Lucas_B6)[i], "vs", levels(group_name_Lucas_B6)[j]))
#    assign(paste0(levels(group_name_Lucas_B6)[i], "vs", levels(group_name_Lucas_B6)[j]), (topTreat(tfit_Lucas_B6, coef = i, n = Inf)))
#    head(eval(as.name(paste0(levels(group_name_Lucas_B6)[i], "vs", levels(group_name_Lucas_B6)[j]))))
#  }
#  
#} 
####creating list of Toptreat-DEGenes####
# par(mfrow = c(1,length(levels(group_name_Lucas_B6))))
Ranked_DEGenes_list <- list()
Ranked_DEGenes_list_tfit <- list()
k = 0
for (i in 1:(length(levels(group_name_Lucas_B6))-1)) {
  for (j in (i+1):length(levels(group_name_Lucas_B6))) {
    k = k+1 
    Ranked_DEGenes_list[[k]] <- topTable(efit_Lucas_B6, coef = k, n = Inf)
    Ranked_DEGenes_list_tfit[[k]] <- topTreat(tfit_Lucas_B6, coef = k, n = Inf)
  }
} 
names(Ranked_DEGenes_list) <- Comparisons_vector
names(Ranked_DEGenes_list_tfit) <- Comparisons_vector
summary(Ranked_DEGenes_list)
dim(Ranked_DEGenes_list[[1]])

#remove all NAs
for (i in 1:length(Comparisons_vector)) {
  Ranked_DEGenes_list[[i]] <- Ranked_DEGenes_list[[i]][!is.na(Ranked_DEGenes_list[[i]]$SYMBOL),]
  Ranked_DEGenes_list_tfit[[i]] <- Ranked_DEGenes_list_tfit[[i]][!is.na(Ranked_DEGenes_list_tfit[[i]]$SYMBOL),]
}
str(Ranked_DEGenes_list)
dim(Ranked_DEGenes_list[[1]])
#list[1]: selects first element, list[[1]]: selects first element of first element, list[1][[1]]: selects first element of first element, list[1][[1]][[1]]: selects 1st element of 1st element of 1st elemnt

#removing duplicated Symbols#
for (i in 1:length(Comparisons_vector)) {
  Ranked_DEGenes_list[[i]] <- Ranked_DEGenes_list[[i]][-which(duplicated(Ranked_DEGenes_list[[i]]$SYMBOL)),]
  Ranked_DEGenes_list_tfit[[i]] <- Ranked_DEGenes_list_tfit[[i]][-which(duplicated(Ranked_DEGenes_list_tfit[[i]]$SYMBOL)),]
  
  if (SaveCSVs == TRUE) {
    write.csv(Ranked_DEGenes_list[[i]], file = paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                      if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
    write.csv(Ranked_DEGenes_list_tfit[[i]], file = paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                           if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".csv"))
  }
}

####GSEA####
Ranked_DEGenes_list_deframed_t <- list()
Ranked_DEGenes_list_deframed_t_upper <- list()

if (FGSEA == TRUE) {
  for (i in 1:length(Comparisons_vector)) {
    
    #list of deframed vectors
    Ranked_DEGenes_list_deframed_t[[i]] <- deframe(Ranked_DEGenes_list_tfit[[i]][,c(2,6)])
    Ranked_DEGenes_list_deframed_t_upper[[i]] <- setNames(Ranked_DEGenes_list_tfit[[i]]$t, toupper(Ranked_DEGenes_list_tfit[[i]]$SYMBOL))
    
    #create directory
    if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/FGSEA", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}))) 
    {dir.create(paste0(Directory_CSVs_DEGenes, "/FGSEA", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), recursive = TRUE) }
    
    #list of downloaded databases for FGSEA analysis
    Databases_GSEA <- list.files("/Users/gordid/Desktop/Databases_GSEA/")
    Databases_GSEA_names <- gsub(".txt|.gmt", "", Databases_GSEA)
    
    for (k in 1:length(Databases_GSEA)) {
      fgsea_results <- fgsea(pathways = gmtPathways(paste0("/Users/gordid/Desktop/Databases_GSEA/", Databases_GSEA[1])), stats = Ranked_DEGenes_list_deframed_t_upper[[1]], nperm = 1000)[order(pval)]
      
      fwrite(fgsea_results, file = paste0(paste0(Directory_CSVs_DEGenes, "/FGSEA", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), "/FGSEA_Lucas_B6_", 
                                          paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                          if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", Databases_GSEA_names[k],
                                          "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".csv"))
      print(Databases_GSEA_names[k])
    }
    
    
    
    if (Mouse == TRUE) {
      #create directory
      if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/FGSEA_GSKB", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}))) 
      {dir.create(paste0(Directory_CSVs_DEGenes, "/FGSEA_GSKB", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), recursive = TRUE) }
      
      
      #setting up gskb databases
      data(mm_GO)
      data(mm_metabolic)
      data(mm_pathway)
      data(mm_TF)
      Databases_GSKB <- list(mm_GO, mm_metabolic, mm_pathway, mm_TF) 
      names(Databases_GSKB) <- c("mm_GO", "mm_metabolic", "mm_pathway", "mm_TF")
      
      #FGSEA with gskb databases
      for (k in 1:length(Databases_GSKB)) {
        fgsea_gskb_results <- fgsea(pathways = Databases_GSKB[[k]], stats = Ranked_DEGenes_list_deframed_t_upper[[i]], nperm = 1000)[order(pval)]
        
        fwrite(fgsea_gskb_results, file = paste0(paste0(Directory_CSVs_DEGenes, "/FGSEA_GSKB", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), "/FGSEA_GSKB_Lucas_B6_", 
                                                 paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                 if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", names(Databases_GSKB)[k],
                                                 "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".csv"))
        
        print(names(Databases_GSKB)[k])
      }
    } 
  }
}

#Database for EnrichR
#install.packages("enrichR")
#library(enrichR)
#listEnrichrDbs()
#remove.packages("enrichR")

if (Mouse == TRUE) {
  enrichR_dbs <- c("WikiPathways_2019_Mouse", "KEGG_2019_Mouse", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "GWAS_Catalog_2019", "TRRUST_Transcription_Factors_2019")
} else {
  enrichR_dbs <- c("WikiPathways_2019_Human", "KEGG_2019_Human", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "GWAS_Catalog_2019", "TRRUST_Transcription_Factors_2019")
}

for (i in 1:length(Comparisons_vector)) {
  
  if (SaveCSVs == TRUE) {
    
    Lucas_B6_enrichR <- read.csv(paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                          if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
    
    if (EnrichR == TRUE) {
      
      Lucas_B6_enrichR_Down <- Lucas_B6_enrichR[which(Lucas_B6_enrichR$logFC < 0), ]
      Lucas_B6_enrichR_Down_0.05 <- Lucas_B6_enrichR_Down[which(Lucas_B6_enrichR_Down$adj.P.Val < 0.05), ]
      dim(Lucas_B6_enrichR_Down_0.05)
      
      if (length(which(Lucas_B6_enrichR$adj.P.Val < 0.05)) > 800) {
        
        Lucas_B6_enrichR_topDE <- enrichr(as.character(Lucas_B6_enrichR$SYMBOL[1:length(which(Lucas_B6_enrichR$adj.P.Val < 0.05))]), enrichR_dbs)
        
      } else {
        
        Lucas_B6_enrichR_topDE <- enrichr(as.character(Lucas_B6_enrichR$SYMBOL[1:800]), enrichR_dbs)
        
        
      }
      
      if (length(Lucas_B6_enrichR_Down_0.05$adj.P.Val) > 10) {
        Lucas_B6_enriched_Down <- enrichr(as.character(Lucas_B6_enrichR_Down_0.05$SYMBOL), enrichR_dbs)
      }
      
      Lucas_B6_enrichR_Up <- Lucas_B6_enrichR[which(Lucas_B6_enrichR$logFC > 0), ]
      Lucas_B6_enrichR_Up_0.05 <- Lucas_B6_enrichR_Up[which(Lucas_B6_enrichR_Up$adj.P.Val < 0.05), ]
      dim(Lucas_B6_enrichR_Up_0.05)
      
      if (length(Lucas_B6_enrichR_Up_0.05$adj.P.Val) > 10) {
        Lucas_B6_enriched_Up <- enrichr(as.character(Lucas_B6_enrichR_Up_0.05$SYMBOL), enrichR_dbs)
      }
      
      for (k in 1:length(enrichR_dbs)) {
        
        if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_topDE"))) { dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_topDE"), recursive = TRUE) }
        
        write.csv(Lucas_B6_enrichR_topDE[[k]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_topDE"), "/EnrichR_topDE_Lucas_B6_", 
                                                               paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                               if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                               "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        
        if (length(Lucas_B6_enrichR_Down_0.05$adj.P.Val) > 10) {
          if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"))) {dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"), recursive = TRUE) }
          
          write.csv(Lucas_B6_enriched_Down[[j]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"), "/EnrichR_Down_Lucas_B6_", 
                                                                 paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                                 if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                                 "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        }
        
        if (length(Lucas_B6_enrichR_Up_0.05$adj.P.Val) > 10) {
          if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"))) {dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"), recursive = TRUE)}
          
          write.csv(Lucas_B6_enriched_Up[[j]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"), "/EnrichR_Up_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                               if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                               "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        }
      }
      
      
      
      
      
    }
  }
} 




#par(mfrow = c(1,3))
#BvsC_Lucas_B6 <- topTreat(tfit_Lucas_B6, coef = 3, n = Inf)
#head(BvsC_Lucas_B6)
#class(BvsC_Lucas_B6)
#dim(BvsC_Lucas_B6)
#BvsC_Lucas_B6_NA_removed <- BvsC_Lucas_B6[!is.na(BvsC_Lucas_B6$SYMBOL),]
#head(BvsC_Lucas_B6_NA_removed$SYMBOL)
#dim(BvsC_Lucas_B6_NA_removed)
#par(mfrow = c(1,1))
#plotMD(tfit_Lucas_B6, column=3, status=dt_Lucas_B6[,3], xlim=c(-8,13), main = "BvsC mean-difference plot")
#glMDPlot(tfit_Lucas_B6, coef=3, status=dt_Lucas_B6, main=colnames(tfit_Lucas_B6)[3], side.main="ENSEMBL", 
#counts=lcpm_Lucas_B6, groups=combined_data_Lucas_B6$samples$group,launch=FALSE)
#main=colnames(tfit_Lucas_B6)[3]

#AvsB_Lucas_B6 <- topTreat(tfit_Lucas_B6, coef = 1, n = Inf)
#head(AvsB_Lucas_B6)
#dim(AvsB_Lucas_B6)
#AvsB_Lucas_B6_NA_removed <- AvsB_Lucas_B6[!is.na(AvsB_Lucas_B6$SYMBOL),]
#head(AvsB_Lucas_B6_NA_removed$SYMBOL)
#plotMD(tfit_Lucas_B6, column = 1, status = dt_Lucas_B6[,1], xlim = c(-8,13), main = "AvsB mean-difference plot")

#AvsC_Lucas_B6 <- topTreat(tfit_Lucas_B6, coef = 2, n = Inf)
#head(AvsC_Lucas_B6)
#AvsC_Lucas_B6_NA_removed <- AvsC_Lucas_B6[!is.na(AvsC_Lucas_B6$SYMBOL),]
#head(AvsC_Lucas_B6_NA_removed$SYMBOL)
#plotMD(tfit_Lucas_B6, column = 2, status = dt_Lucas_B6[,2], xlim = c(-8,13), main = "AvsC mean-difference plot")

####HEATMAPS####

par(mfrow = c(1,1))
for (i in 1:length(Comparisons_vector)) {
  if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/HeatMap_lcpm_Lucas_B6_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                           if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, ".pdf", sep = ""), width=16.25, height=12.5)}
  
  j <- which(v_Lucas_B6$genes$ENSEMBL %in% Ranked_DEGenes_list[[i]]$ENSEMBL[1:100])
  mycol <- colorpanel(1000, "blue", "white", "red")
  heatmap.2(lcpm_Lucas_B6[j,], scale = "row", labRow = v_Lucas_B6$genes$SYMBOL[j], 
            labCol = sample_name_Lucas_B6, col = mycol, trace = "none", Colv = FALSE, #Colv reorders the columns of heatmap
            density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'row', distfun = function(x) as.dist(1-cor(t(x), method = "s")),
            main = paste0(Comparisons_vector[i], " Heat Map"))
  if (SavePlots == TRUE) {dev.off()}
}          


#BvsC_topgenes <- BvsC_Lucas_B6_NA_removed$ENSEMBL[1:80]
#i <- which(v_Lucas_B6$genes$ENSEMBL %in% BvsC_topgenes)
#mycol <- colorpanel(1000, "blue", "white", "red")
#heatmap.2(lcpm_Lucas_B6[i,], scale = "row", labRow = v_Lucas_B6$genes$SYMBOL[i], labCol = sample_name_Lucas_B6, col = mycol, trace = "none",
#         density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'both', distfun = function(x) as.dist(1-cor(t(x), method = "s")), main = "BvsC HeatMap")


#}

#Bar plot of DEGenes (AdjP Value)
# experiment_name <- "Lucas_B6"
# bar_plot_names <- gsub(filter_type, "", list.files("/Users/gordid/Desktop/Lucas_B6/Ranked_DEGenes/")) %>% gsub(pattern = paste0("Ranked-DEGenes_", experiment_name), replacement = "") %>%
#   gsub(pattern = "_[^_]+$", replacement = "") %>% gsub(pattern = "^_", replacement = "")
# bar_plot_names
# 
# Ranked_DEGenes_data_Lucas_B6 <- list()
# Total_DEGenes_Lucas_B6 <- data.frame()
# setwd("/Users/gordid/Desktop/Lucas_B6/Ranked_DEGenes_lfc/")
# for (i in 1:length(list.files())) {
#   Ranked_DEGenes_data_Lucas_B6[[i]] <- read.csv(list.files()[i])
#   print(list.files()[i])
#   Total_DEGenes_Lucas_B6[i,1] <- bar_plot_names[i]
#   Total_DEGenes_Lucas_B6[i,2] <- length(which(Ranked_DEGenes_data_Lucas_B6[[i]]$adj.P.Val[which(Ranked_DEGenes_data_Lucas_B6[[i]]$logFC > 0)] < 0.05))
#   Total_DEGenes_Lucas_B6[i,3] <- length(which(Ranked_DEGenes_data_Lucas_B6[[i]]$adj.P.Val[which(Ranked_DEGenes_data_Lucas_B6[[i]]$logFC < 0)] < 0.05))
# }
# names(Ranked_DEGenes_data_Lucas_B6) <- bar_plot_names
# colnames(Total_DEGenes_Lucas_B6) <- c("Comparisons", "DEG_Up_Counts", "DEG_Down_Counts")
# Total_DEGenes_Lucas_B6 <- Total_DEGenes_Lucas_B6[order(rowSums(Total_DEGenes_Lucas_B6[,c(2,3)])),]
# Total_DEGenes_Lucas_B6$Comparisons <- factor(as.factor(Total_DEGenes_Lucas_B6$Comparisons), levels = unique(as.factor(Total_DEGenes_Lucas_B6$Comparisons)) )
# Total_DEGenes_Lucas_B6
# 
# Total_DEGenes_Lucas_B6 <- melt(Total_DEGenes_Lucas_B6, id.var="Comparisons")
# Total_DEGenes_Lucas_B6
# setwd("/Users/gordid/Desktop/Lucas_B6/")
# 
# if (SavePlots == TRUE) {pdf(file = paste0("/Users/gordid/Desktop/Lucas_B6/All_results", "/Total_DEGenes_AdjP05_lfc_perComparison_BarGraph_Lucas_B6",
#                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".pdf"), width=30, height=12.5)}
# 
# ggplot(Total_DEGenes_Lucas_B6, aes(x = Comparisons, y = value, fill = variable)) +
#   geom_bar(stat = "identity") + xlab("Sample Comparisons") + ylab("DEG Counts (AdjP < 0.05)") + coord_flip() + scale_fill_manual(values=c("red", "blue"))
# 
# if (SavePlots == TRUE) {dev.off()}
# 
# ##Bar plot of DEGenes (P Value)
# Ranked_DEGenes_data_Lucas_B6 <- list()
# Total_DEGenes_Lucas_B6 <- data.frame()
# setwd("/Users/gordid/Desktop/Lucas_B6/Ranked_DEGenes_lfc/")
# for (i in 1:length(list.files())) {
#   Ranked_DEGenes_data_Lucas_B6[[i]] <- read.csv(list.files()[i])
#   print(list.files()[i])
#   Total_DEGenes_Lucas_B6[i,1] <- bar_plot_names[i]
#   Total_DEGenes_Lucas_B6[i,2] <- length(which(Ranked_DEGenes_data_Lucas_B6[[i]]$P.Val[which(Ranked_DEGenes_data_Lucas_B6[[i]]$logFC > 0)] < 0.05))
#   Total_DEGenes_Lucas_B6[i,3] <- length(which(Ranked_DEGenes_data_Lucas_B6[[i]]$P.Val[which(Ranked_DEGenes_data_Lucas_B6[[i]]$logFC < 0)] < 0.05))
# }
# names(Ranked_DEGenes_data_Lucas_B6) <- bar_plot_names
# colnames(Total_DEGenes_Lucas_B6) <- c("Comparisons", "DEG_Up_Counts", "DEG_Down_Counts")
# Total_DEGenes_Lucas_B6 <- Total_DEGenes_Lucas_B6[order(rowSums(Total_DEGenes_Lucas_B6[,c(2,3)])),]
# Total_DEGenes_Lucas_B6$Comparisons <- factor(as.factor(Total_DEGenes_Lucas_B6$Comparisons), levels = unique(as.factor(Total_DEGenes_Lucas_B6$Comparisons)) )
# Total_DEGenes_Lucas_B6
# 
# Total_DEGenes_Lucas_B6 <- melt(Total_DEGenes_Lucas_B6, id.var="Comparisons")
# Total_DEGenes_Lucas_B6
# setwd("/Users/gordid/Desktop/Lucas_B6/")
# 
# if (SavePlots == TRUE) {pdf(file = paste0("/Users/gordid/Desktop/Lucas_B6/All_results", "/Total_DEGenes_P05_lfc_perComparison_BarGraph_Lucas_B6",
#                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".pdf"), width=30, height=12.5)}
# 
# ggplot(Total_DEGenes_Lucas_B6, aes(x = Comparisons, y = value, fill = variable)) +
#   geom_bar(stat = "identity") + xlab("Sample Comparisons") + ylab("DEG Counts (P < 0.05)") + coord_flip() + scale_fill_manual(values=c("red", "blue"))
# 
# if (SavePlots == TRUE) {dev.off()}
