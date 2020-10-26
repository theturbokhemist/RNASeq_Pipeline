#read in data file on experiment
setwd("/Users/gordid/Desktop/Foley_Runs/")

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
Rename_file_extensions("/Users/gordid/Desktop/Foley_Runs/genes_results/")


#reading in design CSV
experiment_Foley_Runs <- read.csv(paste0("/Users/gordid/Desktop/Foley_Runs/", "Foley_Runs_experimental_design", ".csv"))
experiment_Foley_Runs

#creating vector of data files
files_Foley_Runs_unsorted <- list.files("/Users/gordid/Desktop/Foley_Runs/genes_results/")
files_Foley_Runs_unsorted

#Putting data files character vector into the order of the sample IDS in the design experiment CSV
files_Foley_Runs <- c()

for (i in 1:length(files_Foley_Runs_unsorted)) {
  for (j in 1:length(experiment_Foley_Runs$Sample.Name)) {
    if (grepl(as.character(experiment_Foley_Runs$Sample.Name)[j], files_Foley_Runs_unsorted[i]) == TRUE) {
      files_Foley_Runs[j] <- files_Foley_Runs_unsorted[i]
    }
  }
}
files_Foley_Runs
#Files need to be in the same order as the samples in the excel file. Find a way to call by name instead of index.

####Save plots?####
SavePlots = TRUE
SaveCSVs = TRUE
EnrichR = FALSE
FGSEA = FALSE
NetAct = FALSE
Calcs = TRUE
RACIPE = FALSE
####What Organism?####
Mouse = TRUE 

####Comparing just 2 groups?####
TwoGroups = TRUE

experiment_Foley_Runs$Group.Name

if (TwoGroups == TRUE) {
  
  Group1 = "CD4weeksSedC"
  Group2 = "WD4weeksSedC"

  WhichGroups = c(which(experiment_Foley_Runs$Group.Name == Group1),which(experiment_Foley_Runs$Group.Name == Group2))
  
  experiment_Foley_Runs <- experiment_Foley_Runs[WhichGroups,]
  
  files_Foley_Runs <- files_Foley_Runs[WhichGroups]
}


####Are there Outliers?####
outliers = FALSE

if (outliers == TRUE) {
  Outlier1 = "17979H"
  Outlier2 = ""
  
  outliers_index <- c(which(experiment_Foley_Runs$Sample.Name == Outlier1),which(experiment_Foley_Runs$Sample.Name == Outlier2))
  
  outliers_name <- as.character(experiment_Foley_Runs$Sample.Name[outliers_index])
  
  experiment_Foley_Runs <- experiment_Foley_Runs[-outliers_index,]
  
  files_Foley_Runs <- files_Foley_Runs[-outliers_index]
  
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
setwd("/Users/gordid/Desktop/Foley_Runs/genes_results/")
list.files()
file_1_Foley_Runs <- read.delim(files_Foley_Runs[1], header=TRUE)
head(file_1_Foley_Runs)
str(file_1_Foley_Runs)
dim(file_1_Foley_Runs)
setwd("/Users/gordid/Desktop/Foley_Runs")

####Renaming columns and creating groups####

group_name_Foley_Runs <- as.character(experiment_Foley_Runs$Group.Name)
group_name_Foley_Runs <- factor(as.factor(group_name_Foley_Runs), levels = unique(as.factor(group_name_Foley_Runs)) ) #preserves order of factors
group_name_Foley_Runs
levels(group_name_Foley_Runs)


sample_name_Foley_Runs <- as.character(experiment_Foley_Runs$Sample.Name)
sample_name_Foley_Runs

Tissue_Foley_Runs <- factor(as.factor(as.character(experiment_Foley_Runs$Tissue)), levels = unique(as.factor(as.character(experiment_Foley_Runs$Tissue))) )
Tissue_Foley_Runs
length(levels(Tissue_Foley_Runs))

Diet.Duration_Foley_Runs <- factor(as.factor(as.character(experiment_Foley_Runs$Diet.Duration)), levels = unique(as.factor(as.character(experiment_Foley_Runs$Diet.Duration))) )
Diet.Duration_Foley_Runs

Remarks_Foley_Runs <- factor(as.factor(as.character(experiment_Foley_Runs$Remarks)), levels = unique(as.factor(as.character(experiment_Foley_Runs$Remarks))) )
Remarks_Foley_Runs

####creating vector of gene names and vector of possible sample comparisons####

Comparisons_vector <- c()
for (i in 1:(length(levels(group_name_Foley_Runs))-1)) {
  for (j in (i+1):length(levels(group_name_Foley_Runs))) {
    
    Comparisons_vector <- append(Comparisons_vector, paste0(levels(group_name_Foley_Runs)[i], "vs", levels(group_name_Foley_Runs)[j]))
    
  }
}

Comparisons_vector
length(Comparisons_vector)

#Creating new directorys####

#plots/RNASeq folder
Directory_Plots_RNASeq <- paste0("/Users/gordid/Desktop/Foley_Runs/All_results/Plots/RNASeq/",
                                 (if(TwoGroups == TRUE) { 
                                   paste0(levels(group_name_Foley_Runs)[1], "vs", levels(group_name_Foley_Runs)[2],
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

Directory_Plots_NetAct_RACIPE <- paste0("/Users/gordid/Desktop/Foley_Runs/All_results/Plots/NetAct_RACIPE/",
                                        (if(TwoGroups == TRUE) { 
                                          paste0(levels(group_name_Foley_Runs)[1], "vs", levels(group_name_Foley_Runs)[2],
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
Directory_CSVs_DEGenes <- paste0("/Users/gordid/Desktop/Foley_Runs/All_results/CSVs/DEGenes/", 
                                 (if(TwoGroups == TRUE) {
                                   paste0(levels(group_name_Foley_Runs)[1], "vs", levels(group_name_Foley_Runs)[2],
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
Directory_CSVs_NetAct_RACIPE <- paste0("/Users/gordid/Desktop/Foley_Runs/All_results/CSVs/NetAct_RACIPE/",
                                       (if(TwoGroups == TRUE) { 
                                         paste0(levels(group_name_Foley_Runs)[1], "vs", levels(group_name_Foley_Runs)[2],
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
setwd("/Users/gordid/Desktop/Foley_Runs/genes_results/")
combined_data_Foley_Runs <- readDGE(files_Foley_Runs, columns = c(1, 5))
dim(combined_data_Foley_Runs) #showing dimensions of "counts" matrix
class(combined_data_Foley_Runs) #showing class of the whole object
summary(combined_data_Foley_Runs) #showing names, size, class, and mode of each element of the object, not sure why class is "none" for counts
str(combined_data_Foley_Runs) #showing the structure of the object
head(combined_data_Foley_Runs$counts)
colnames(combined_data_Foley_Runs) <- sample_name_Foley_Runs
setwd("/Users/gordid/Desktop/Foley_Runs")

####Assigning experimental info to categories in the DGEList obect####
combined_data_Foley_Runs$samples$group <- group_name_Foley_Runs
combined_data_Foley_Runs$samples

#combined_data_Foley_Runs$samples$remarks <- remarks_Foley_Runs

####Create new genes matrix by retrieving gene information from the mus library. The gene information of their SYMBOL and TXCHROM#### 
#number of gene IDs
gene_id_Foley_Runs <- rownames(combined_data_Foley_Runs)
head(gene_id_Foley_Runs)
length(gene_id_Foley_Runs)
#number of duplicated ENSEMBL IDS PRE select
dup_genes_Foley_Runs_PRESelect <- which(duplicated(gene_id_Foley_Runs))
length(dup_genes_Foley_Runs_PRESelect)

columns(Mus.musculus)
if (Mouse == TRUE) {
  genes_Foley_Runs <- biomaRt::select(Mus.musculus, keys = gene_id_Foley_Runs, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
} else {
  genes_Foley_Runs <- biomaRt::select(Homo.sapiens, keys = gene_id_Foley_Runs, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
  
}
dim(genes_Foley_Runs)
head(genes_Foley_Runs)
length(genes_Foley_Runs$ENSEMBL) - length(gene_id_Foley_Runs)

####symbols are duplicated too often and so I dont remove genes based on duplicated symbols (cause of NA mostly, will remove NAs later)####
head((which(duplicated(genes_Foley_Runs$SYMBOL))))
length(which(duplicated(genes_Foley_Runs$SYMBOL)))
genes_Foley_Runs[head((which(duplicated(genes_Foley_Runs$SYMBOL)))),]

#number of duplicated ENSEMBL IDS POST select
dup_genes_Foley_Runs <- which(duplicated(genes_Foley_Runs$ENSEMBL))
length(dup_genes_Foley_Runs)
head(dup_genes_Foley_Runs)
genes_Foley_Runs[c(head(which(duplicated(genes_Foley_Runs$ENSEMBL))), head(which(duplicated(genes_Foley_Runs$ENSEMBL))) - 1),]

####remove duplicated genes based on duplicated ENSEMBL ID####
genes_Foley_Runs_final <- genes_Foley_Runs[-dup_genes_Foley_Runs, ]
dim(genes_Foley_Runs_final)
dim(combined_data_Foley_Runs$counts)

#Number of removed Symbol IDS must be equal to or less than the number of removed ENSEMBL IDs because during the select process, some ENSEMBL ID's may match to two different Symbol IDs.
length(which(duplicated(genes_Foley_Runs$SYMBOL))) - length(which(duplicated(genes_Foley_Runs_final$SYMBOL)))
length(dup_genes_Foley_Runs)

#assign data frame of gene annotations to the combined_data object
combined_data_Foley_Runs$genes <- genes_Foley_Runs_final 
summary(combined_data_Foley_Runs)

#RAW COUNTS DATA
if (SaveCSVs == TRUE) {
  
  write.csv(cbind2(combined_data_Foley_Runs$counts, combined_data_Foley_Runs$genes), file = paste0(Directory_CSVs_DEGenes,
                                                                                                   "/Raw_Counts_Data_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), ".csv"))
  
}

####CPM normalization on count data####
cpm_Foley_Runs <- cpm(combined_data_Foley_Runs)
lcpm_Foley_Runs <-cpm(combined_data_Foley_Runs, log = TRUE)

L_Foley_Runs <- mean(combined_data_Foley_Runs$samples$lib.size) * 1e-6 #calculation of L for the prior count.
M_Foley_Runs <- median(combined_data_Foley_Runs$samples$lib.size) * 1e-6
c(L_Foley_Runs, M_Foley_Runs)
head(lcpm_Foley_Runs)

table(rowSums(combined_data_Foley_Runs$counts==0)==length(colnames(lcpm_Foley_Runs)))

###FILTER THE GENES####
if (CPMFilter == TRUE) {
  kept_expression_Foley_Runs <- filterByExpr(combined_data_Foley_Runs, group = combined_data_Foley_Runs$samples$group) #filterByExpr function
  # keeps genes with 10 read counts or more in a minimum number of samples, where the minimum number of samples is chosen according 
  #to the minimum group sample size. 
}  else {
  #filter by counts, rather than by cpm
  kept_expression_Foley_Runs <- (rowSums(combined_data_Foley_Runs$counts > 10) >= 1)
}
length(kept_expression_Foley_Runs)
head(kept_expression_Foley_Runs)

#removing all the lowly expressed genes from the combined_data object. When you subset a DGEList and specify keep.lib.sizes=FALSE,
#the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
filtered_combined_data_Foley_Runs <- combined_data_Foley_Runs[kept_expression_Foley_Runs, , keep.lib.sizes = FALSE] 
dim(filtered_combined_data_Foley_Runs)
dim(combined_data_Foley_Runs)

####Removing duplicated gene symbols for the NetAct analysis####
Combined_data_Foley_Runs_NetAct <- filtered_combined_data_Foley_Runs
rownames(Combined_data_Foley_Runs_NetAct$counts) <- Combined_data_Foley_Runs_NetAct$genes$SYMBOL
dim(Combined_data_Foley_Runs_NetAct)
#Number of Symbols that are NA
dup_genes_Foley_Runs_NA <- which(is.na(filtered_combined_data_Foley_Runs$genes$SYMBOL))
length(dup_genes_Foley_Runs_NA)
#Number of duplicated Symbols that are not NA
length(which(duplicated(na.omit(filtered_combined_data_Foley_Runs$genes$SYMBOL))))
#Number of duplicated symbols + all NAs
length(dup_genes_Foley_Runs_NA) + length(which(duplicated(na.omit(filtered_combined_data_Foley_Runs$genes$SYMBOL))))
# of genes theoretically left over at removing of all duplicated
length(filtered_combined_data_Foley_Runs$genes$SYMBOL) -(length(dup_genes_Foley_Runs_NA) + length(which(duplicated(na.omit(filtered_combined_data_Foley_Runs$genes$SYMBOL)))))

#removing all NAs
Combined_data_Foley_Runs_NetAct <-Combined_data_Foley_Runs_NetAct[-dup_genes_Foley_Runs_NA, ]
dim(Combined_data_Foley_Runs_NetAct)

#removing all duplicated symbols
dup_genes_Foley_Runs_SYMBOL <- which(duplicated(Combined_data_Foley_Runs_NetAct$genes$SYMBOL))
length(dup_genes_Foley_Runs_SYMBOL)
Combined_data_Foley_Runs_NetAct <- Combined_data_Foley_Runs_NetAct[-dup_genes_Foley_Runs_SYMBOL, ]
dim(Combined_data_Foley_Runs_NetAct)

####Plotting filtered Data####
lcpm.cutoff_Foley_Runs <- log2(10/M_Foley_Runs + 2/L_Foley_Runs)
lcpm.cutoff_Foley_Runs
nsamples <- ncol(combined_data_Foley_Runs$counts)
col <- brewer.pal(nsamples, "Paired")

if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/FilteredData_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                  ".pdf", sep = ""))}
par(mfrow = c(1,2))
plot(density(lcpm_Foley_Runs[, 1]), col = col[1], lwd=2, ylim=c(0,0.46), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_Foley_Runs, lty=3)
title(main="A. Raw data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_Foley_Runs[, i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

combined_data_Foley_Runs <- filtered_combined_data_Foley_Runs
lcpm_Foley_Runs <- cpm(combined_data_Foley_Runs, log = TRUE)

plot(density(lcpm_Foley_Runs[, 1]), col = col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_Foley_Runs, lty=3)
title(main="B. Filtered data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_Foley_Runs[, i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
if (SavePlots == TRUE) {dev.off() }


####Calculate normalization factors and creates new DGEList object with these new normalization factors inserted. THEN GENERATING BOX PLOT####
colnames(lcpm_Foley_Runs) <- sample_name_Foley_Runs
combined_data_Foley_Runs$samples$norm.factors

if (SavePlots == TRUE) {pdf(paste0(Directory_Plots_RNASeq, "/TMM_Normalized_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), ".pdf"))}


par(mfrow = c(1,2))
boxplot(lcpm_Foley_Runs, las = 2, col = col, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")

combined_data_Foley_Runs = calcNormFactors(combined_data_Foley_Runs, method = "TMM") #calculates normalization factors and creates new DGEList object
#these normalization factors are important because if there is a sample where expression level is different accross the board, that may have been caused
#by an eternal factor that is not of biological interest. TMM helps normalize for that, using the assumption that should have a similar range
#and distribution of expression values.
lcpm_Foley_Runs <- cpm(combined_data_Foley_Runs, log = TRUE) 
colnames(lcpm_Foley_Runs) <- sample_name_Foley_Runs
head(lcpm_Foley_Runs)

combined_data_Foley_Runs$samples$norm.factors
boxplot(lcpm_Foley_Runs, las = 2, col = col, main="")
title(main="B.Normalised with TMM data", ylab="Log-cpm")
if (SavePlots == TRUE) {dev.off() }

####PCA and PCoA Plots####
if (SavePlots == TRUE) {pdf(paste0(Directory_Plots_RNASeq, "/PCA_bygroup_1vs2_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                   ".pdf"))}
color_groups_Foley_Runs <- combined_data_Foley_Runs$samples$group
levels(color_groups_Foley_Runs) <- brewer.pal(nlevels(color_groups_Foley_Runs), "Set1")
color_groups_Foley_Runs <- as.character(color_groups_Foley_Runs)
color_groups_Foley_Runs

combined_data_Foley_Runs$samples$names <- sample_name_Foley_Runs
combined_data_Foley_Runs$samples

par(mfrow = c(1,1))
#plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$names, col = color_groups_Foley_Runs, main = "PCoA")
plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$group, col = color_groups_Foley_Runs, gene.selection="common", 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))
if (SavePlots == TRUE) {dev.off() }

par(mfrow = c(1,2))
plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$names, col = color_groups_Foley_Runs, dim = c(2,3), 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCoA")) 
plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$names, col = color_groups_Foley_Runs, gene.selection="common", dim = c(2,3),
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))

par(mfrow = c(1,1))
plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$names, col = color_groups_Foley_Runs, gene.selection="common", 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))


glMDSPlot(lcpm_Foley_Runs, labels=paste(combined_data_Foley_Runs$samples$names, combined_data_Foley_Runs$samples$group, sep="_"), 
          groups = combined_data_Foley_Runs$samples[,c(2,5)], 
          gene.selection="common", main = paste0(sub("_[^_]+$*" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"), launch=TRUE, path = Directory_Plots_RNASeq,
          html = paste0("InteractivePCA_ScreePlot_1vs2_Foley_Runs", "_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")))

#scree.plot(plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$names, col = color_groups_Foley_Runs, gene.selection="common", main = "PCA"))
#a <- (plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$names, col = color_groups_Foley_Runs, main = "PCoA"))
#screeplot(a, xlim = c(-8,13), ylim=c(0,0.46))
####PCA Plots for other phenotypic data and meta-data####
#color_age_Foley_Runs <- combined_data_Foley_Runs$samples$age
#levels(color_age_Foley_Runs) <- brewer.pal(nlevels(color_age_Foley_Runs), "Set2")
#color_age_Foley_Runs <- as.character(color_age_Foley_Runs)
#color_age_Foley_Runs
#plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$age, col = color_groups_Foley_Runs)
#plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$age, col = color_groups_Foley_Runs, gene.selection="common")

#color_x2dg_Foley_Runs <- combined_data_Foley_Runs$samples$x2dg
#levels(color_x2dg_Foley_Runs) <- brewer.pal(nlevels(color_x2dg_Foley_Runs), "Set3")
#color_x2dg_Foley_Runs <- as.character(color_x2dg_Foley_Runs)
#color_x2dg_Foley_Runs
#plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$x2dg, col = color_groups_Foley_Runs)
#plotMDS(lcpm_Foley_Runs, labels = combined_data_Foley_Runs$samples$x2dg, col = color_groups_Foley_Runs, gene.selection="common")

####Building Design Matrix####

design_matrix_Foley_Runs_no.int <- model.matrix(~0+combined_data_Foley_Runs$samples$group) #no intercept
colnames(design_matrix_Foley_Runs_no.int) = gsub("combined_data_Foley_Runs$samples$group", "", colnames(design_matrix_Foley_Runs_no.int), fixed = TRUE)
design_matrix_Foley_Runs_no.int

design_matrix_Foley_Runs_int <- model.matrix(~combined_data_Foley_Runs$samples$group) #with intercept
colnames(design_matrix_Foley_Runs_int) = gsub("combined_data_Foley_Runs$samples$group", "", colnames(design_matrix_Foley_Runs_int), fixed = TRUE)
design_matrix_Foley_Runs_int


if (TwoGroups == TRUE) {
  contr.matrix_no.int <- makeContrasts(
    paste(levels(group_name_Foley_Runs)[1]," - ", levels(group_name_Foley_Runs)[2], sep = ""),
    levels = colnames(design_matrix_Foley_Runs_no.int))
  colnames(contr.matrix_no.int) <- paste(levels(group_name_Foley_Runs)[1],"vs", levels(group_name_Foley_Runs)[2], sep = "")
  contr.matrix_no.int
} else {
  contrast_vector <- c()
  column_name_vector <- c()
  for (i in 1:(length(levels(group_name_Foley_Runs))-1)) {
    
    for (j in (i+1):length(levels(group_name_Foley_Runs))) {
      
      contrast_vector <- append(contrast_vector, paste(levels(group_name_Foley_Runs)[i]," - ", levels(group_name_Foley_Runs)[j], sep = ""))
      column_name_vector <- append(column_name_vector, paste(levels(group_name_Foley_Runs)[i], "vs", levels(group_name_Foley_Runs)[j], sep = ""))
    }
    
  }
  
  contr.matrix_no.int <- makeContrasts(
    contrasts = contrast_vector,
    levels = colnames(design_matrix_Foley_Runs_no.int)
  )
  colnames(contr.matrix_no.int) <- column_name_vector
  contr.matrix_no.int
}


#contr.matrix_no.int <- makeContrasts(
#AvsB = "A- B",
#AvsC = A - C,
#BvsC = B - C,
#levels = colnames(design_matrix_Foley_Runs_no.int))
#contr.matrix_no.int
#colnames(contr.matrix_no.int) <- c("a", "b", "c")

#contr.matrix_no.int <- makeContrasts(
#contrasts = contrast_vector,
#levels = colnames(design_matrix_Foley_Runs_no.int))
#contr.matrix_no.int

#contr.matrix_int <- makeContrasts( #Intercept column has replcaced "A", not sure how contrast matrix works now
# AvsB = A - B,
# AvsC = A - C,
# BvsC = B - C,
#levels = colnames(design_matrix_Foley_Runs_int))
#contr.matrix_int

####VOOM####
adjP = 0.05
if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Voom_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                  ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,2))
v_Foley_Runs <- voom(combined_data_Foley_Runs, design_matrix_Foley_Runs_no.int, plot=TRUE)
str(v_Foley_Runs)
dim(v_Foley_Runs$genes)

vfit_Foley_Runs <- lmFit(v_Foley_Runs, design_matrix_Foley_Runs_no.int)
summary(decideTests(vfit_Foley_Runs))
str(vfit_Foley_Runs)

vfit_Foley_Runs_contrasts <- contrasts.fit(vfit_Foley_Runs, contrasts=contr.matrix_no.int) 
summary(decideTests(vfit_Foley_Runs_contrasts))

efit_Foley_Runs <- eBayes(vfit_Foley_Runs_contrasts)
summary(decideTests(efit_Foley_Runs))
de_Foley_Runs <- decideTests(efit_Foley_Runs)

plotSA(efit_Foley_Runs)
title(main="voom: Mean-variance trend-Normalized")

if (SavePlots == TRUE) {dev.off()}


####Filter further with logfold change and Venn Diagram####
#Figure out the criteria data has to choose appropriate lfc
lfc = 0.1
tfit_Foley_Runs <- treat(vfit_Foley_Runs_contrasts, lfc=lfc)
dt_Foley_Runs <- decideTests(tfit_Foley_Runs)
summary(dt_Foley_Runs)


topTable_p <-  topTable(efit_Foley_Runs, n=Inf); head(topTable_p)
topTreat_p <- topTreat(efit_Foley_Runs, n=Inf); head(topTreat_p)
topTreat_lfc <- topTreat(tfit_Foley_Runs, n=Inf); head(topTreat_lfc)

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
  de.common_Foley_Runs_2_comparisons <- which(dt_Foley_Runs[,1]!=0 & dt_Foley_Runs[,2]!=0)
  print(length(de.common_Foley_Runs_2_comparisons))
  
  de.common_Foley_Runs_3_comparisons <- which(dt_Foley_Runs[,1]!=0 & dt_Foley_Runs[,2]!=0 & dt_Foley_Runs[,3]!=0)
  print(length(de.common_Foley_Runs_3_comparisons))
  
  if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Venn_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                    ".pdf", sep = ""), width=13, height=10)}
  par(mfrow = c(1,1))
  vennDiagram(dt_Foley_Runs[,1:3], circle.col=c("turquoise", "salmon", "orange"))
  if (SavePlots == TRUE) {dev.off()}
}

if (SaveCSVs == TRUE) {
  write.fit(tfit_Foley_Runs, dt_Foley_Runs, file = paste(Directory_CSVs_DEGenes, "/tfit_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                         "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".txt", sep = ""))
  
  tfit_Foley_Runs.txt = read.table(file = paste(Directory_CSVs_DEGenes, "/tfit_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                                gsub(".", "", as.character(lfc), fixed = TRUE),
                                                ".txt", sep = ""))
  
  write.csv(tfit_Foley_Runs.txt, file = paste(Directory_CSVs_DEGenes, "/tfit_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                              gsub(".", "", as.character(lfc), fixed = TRUE),
                                              ".csv", sep = ""))
}

####Mean-difference plots and arranging DE genes####
if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/Mean-Difference-Plots_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                         gsub(".", "", as.character(lfc), fixed = TRUE),
                                         ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,length(Comparisons_vector)))
for (i in 1:length(Comparisons_vector)) {
  plotMD(efit_Foley_Runs, column = i, status = de_Foley_Runs[,i], xlim = c(-8,13), main = paste0(Comparisons_vector[i], " Mean-Difference Plot"))
  
  glMDPlot(efit_Foley_Runs, coef=i, status=de_Foley_Runs, main=colnames(efit_Foley_Runs)[i], side.main="ENSEMBL", counts=lcpm_Foley_Runs,
           groups=combined_data_Foley_Runs$samples$group,launch=FALSE, 
           path = Directory_Plots_RNASeq, html = paste0("Interactive_MeanDifference_Plot_Foley_Runs", "_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), 
                                                        if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE)))
  
}
if (SavePlots == TRUE){ dev.off()}

#par(mfrow = c(1,length(levels(group_name_Foley_Runs))))
#Ranked_DEGenes_list <- list()
#for (i in 1:(length(levels(group_name_Foley_Runs))-1)) {
#  for (j in (i+1):length(levels(group_name_Foley_Runs))) {
#    print(paste0(levels(group_name_Foley_Runs)[i], "vs", levels(group_name_Foley_Runs)[j]))
#    assign(paste0(levels(group_name_Foley_Runs)[i], "vs", levels(group_name_Foley_Runs)[j]), (topTreat(tfit_Foley_Runs, coef = i, n = Inf)))
#    head(eval(as.name(paste0(levels(group_name_Foley_Runs)[i], "vs", levels(group_name_Foley_Runs)[j]))))
#  }
#  
#} 
####creating list of Toptreat-DEGenes####
# par(mfrow = c(1,length(levels(group_name_Foley_Runs))))
Ranked_DEGenes_list <- list()
Ranked_DEGenes_list_tfit <- list()
k = 0
for (i in 1:(length(levels(group_name_Foley_Runs))-1)) {
  for (j in (i+1):length(levels(group_name_Foley_Runs))) {
    k = k+1 
    Ranked_DEGenes_list[[k]] <- topTable(efit_Foley_Runs, coef = k, n = Inf)
    Ranked_DEGenes_list_tfit[[k]] <- topTreat(tfit_Foley_Runs, coef = k, n = Inf)
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
    write.csv(Ranked_DEGenes_list[[i]], file = paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                      if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
    write.csv(Ranked_DEGenes_list_tfit[[i]], file = paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
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
      
      fwrite(fgsea_results, file = paste0(paste0(Directory_CSVs_DEGenes, "/FGSEA", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), "/FGSEA_Foley_Runs_", 
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
        
        fwrite(fgsea_gskb_results, file = paste0(paste0(Directory_CSVs_DEGenes, "/FGSEA_GSKB", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), "/FGSEA_GSKB_Foley_Runs_", 
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
    
    Foley_Runs_enrichR <- read.csv(paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                          if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
    
    if (EnrichR == TRUE) {
      
      Foley_Runs_enrichR_Down <- Foley_Runs_enrichR[which(Foley_Runs_enrichR$logFC < 0), ]
      Foley_Runs_enrichR_Down_0.05 <- Foley_Runs_enrichR_Down[which(Foley_Runs_enrichR_Down$adj.P.Val < 0.05), ]
      dim(Foley_Runs_enrichR_Down_0.05)
      
      if (length(which(Foley_Runs_enrichR$adj.P.Val < 0.05)) > 800) {
        
        Foley_Runs_enrichR_topDE <- enrichr(as.character(Foley_Runs_enrichR$SYMBOL[1:length(which(Foley_Runs_enrichR$adj.P.Val < 0.05))]), enrichR_dbs)
        
      } else {
        
        Foley_Runs_enrichR_topDE <- enrichr(as.character(Foley_Runs_enrichR$SYMBOL[1:800]), enrichR_dbs)
        
        
      }
      
      if (length(Foley_Runs_enrichR_Down_0.05$adj.P.Val) > 10) {
        Foley_Runs_enriched_Down <- enrichr(as.character(Foley_Runs_enrichR_Down_0.05$SYMBOL), enrichR_dbs)
      }
      
      Foley_Runs_enrichR_Up <- Foley_Runs_enrichR[which(Foley_Runs_enrichR$logFC > 0), ]
      Foley_Runs_enrichR_Up_0.05 <- Foley_Runs_enrichR_Up[which(Foley_Runs_enrichR_Up$adj.P.Val < 0.05), ]
      dim(Foley_Runs_enrichR_Up_0.05)
      
      if (length(Foley_Runs_enrichR_Up_0.05$adj.P.Val) > 10) {
        Foley_Runs_enriched_Up <- enrichr(as.character(Foley_Runs_enrichR_Up_0.05$SYMBOL), enrichR_dbs)
      }
      
      for (k in 1:length(enrichR_dbs)) {
        
        if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_topDE"))) { dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_topDE"), recursive = TRUE) }
        
        write.csv(Foley_Runs_enrichR_topDE[[k]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_topDE"), "/EnrichR_topDE_Foley_Runs_", 
                                                               paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                               if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                               "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        
        if (length(Foley_Runs_enrichR_Down_0.05$adj.P.Val) > 10) {
          if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"))) {dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"), recursive = TRUE) }
          
          write.csv(Foley_Runs_enriched_Down[[j]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"), "/EnrichR_Down_Foley_Runs_", 
                                                                 paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                                 if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                                 "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        }
        
        if (length(Foley_Runs_enrichR_Up_0.05$adj.P.Val) > 10) {
          if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"))) {dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"), recursive = TRUE)}
          
          write.csv(Foley_Runs_enriched_Up[[j]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"), "/EnrichR_Up_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                               if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                               "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        }
      }
      
      
      
      
      
    }
  }
} 




#par(mfrow = c(1,3))
#BvsC_Foley_Runs <- topTreat(tfit_Foley_Runs, coef = 3, n = Inf)
#head(BvsC_Foley_Runs)
#class(BvsC_Foley_Runs)
#dim(BvsC_Foley_Runs)
#BvsC_Foley_Runs_NA_removed <- BvsC_Foley_Runs[!is.na(BvsC_Foley_Runs$SYMBOL),]
#head(BvsC_Foley_Runs_NA_removed$SYMBOL)
#dim(BvsC_Foley_Runs_NA_removed)
#par(mfrow = c(1,1))
#plotMD(tfit_Foley_Runs, column=3, status=dt_Foley_Runs[,3], xlim=c(-8,13), main = "BvsC mean-difference plot")
#glMDPlot(tfit_Foley_Runs, coef=3, status=dt_Foley_Runs, main=colnames(tfit_Foley_Runs)[3], side.main="ENSEMBL", 
#counts=lcpm_Foley_Runs, groups=combined_data_Foley_Runs$samples$group,launch=FALSE)
#main=colnames(tfit_Foley_Runs)[3]

#AvsB_Foley_Runs <- topTreat(tfit_Foley_Runs, coef = 1, n = Inf)
#head(AvsB_Foley_Runs)
#dim(AvsB_Foley_Runs)
#AvsB_Foley_Runs_NA_removed <- AvsB_Foley_Runs[!is.na(AvsB_Foley_Runs$SYMBOL),]
#head(AvsB_Foley_Runs_NA_removed$SYMBOL)
#plotMD(tfit_Foley_Runs, column = 1, status = dt_Foley_Runs[,1], xlim = c(-8,13), main = "AvsB mean-difference plot")

#AvsC_Foley_Runs <- topTreat(tfit_Foley_Runs, coef = 2, n = Inf)
#head(AvsC_Foley_Runs)
#AvsC_Foley_Runs_NA_removed <- AvsC_Foley_Runs[!is.na(AvsC_Foley_Runs$SYMBOL),]
#head(AvsC_Foley_Runs_NA_removed$SYMBOL)
#plotMD(tfit_Foley_Runs, column = 2, status = dt_Foley_Runs[,2], xlim = c(-8,13), main = "AvsC mean-difference plot")

####HEATMAPS####

par(mfrow = c(1,1))
for (i in 1:length(Comparisons_vector)) {
  if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/HeatMap_lcpm_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                           if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, ".pdf", sep = ""), width=16.25, height=12.5)}
  
  j <- which(v_Foley_Runs$genes$ENSEMBL %in% Ranked_DEGenes_list[[i]]$ENSEMBL[1:100])
  mycol <- colorpanel(1000, "blue", "white", "red")
  heatmap.2(lcpm_Foley_Runs[j,], scale = "row", labRow = v_Foley_Runs$genes$SYMBOL[j], 
            labCol = sample_name_Foley_Runs, col = mycol, trace = "none", Colv = FALSE, #Colv reorders the columns of heatmap
            density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'row', distfun = function(x) as.dist(1-cor(t(x), method = "s")),
            main = paste0(Comparisons_vector[i], " Heat Map"))
  if (SavePlots == TRUE) {dev.off()}
}          


#BvsC_topgenes <- BvsC_Foley_Runs_NA_removed$ENSEMBL[1:80]
#i <- which(v_Foley_Runs$genes$ENSEMBL %in% BvsC_topgenes)
#mycol <- colorpanel(1000, "blue", "white", "red")
#heatmap.2(lcpm_Foley_Runs[i,], scale = "row", labRow = v_Foley_Runs$genes$SYMBOL[i], labCol = sample_name_Foley_Runs, col = mycol, trace = "none",
#         density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'both', distfun = function(x) as.dist(1-cor(t(x), method = "s")), main = "BvsC HeatMap")


#}

##Bar plot of DEGenes (P Value)
# experiment_name <- "Foley_Runs"
# bar_plot_names <- gsub(filter_type, "", list.files("/Users/gordid/Desktop/Foley_Runs/Ranked_DEGenes/")) %>% gsub(pattern = paste0("Ranked-DEGenes_", experiment_name), replacement = "") %>%
#   gsub(pattern = "_[^_]+$", replacement = "") %>% gsub(pattern = "^_", replacement = "")
# bar_plot_names
# 
# Ranked_DEGenes_data_Foley_Runs <- list()
# Total_DEGenes_Foley_Runs <- data.frame()
# setwd("/Users/gordid/Desktop/Foley_Runs/Ranked_DEGenes/")
# for (i in 1:length(list.files())) {
#   Ranked_DEGenes_data_Foley_Runs[[i]] <- read.csv(list.files()[i])
#   print(list.files()[i])
#   Total_DEGenes_Foley_Runs[i,1] <- bar_plot_names[i]
#   Total_DEGenes_Foley_Runs[i,2] <- length(which(Ranked_DEGenes_data_Foley_Runs[[i]]$P.Value[which(Ranked_DEGenes_data_Foley_Runs[[i]]$logFC > 0)] < 0.05))
#   Total_DEGenes_Foley_Runs[i,3] <- length(which(Ranked_DEGenes_data_Foley_Runs[[i]]$P.Value[which(Ranked_DEGenes_data_Foley_Runs[[i]]$logFC < 0)] < 0.05))
# }
# names(Ranked_DEGenes_data_Foley_Runs) <- bar_plot_names
# colnames(Total_DEGenes_Foley_Runs) <- c("Comparisons", "DEG_Up_Counts", "DEG_Down_Counts")
# Total_DEGenes_Foley_Runs <- Total_DEGenes_Foley_Runs[order(rowSums(Total_DEGenes_Foley_Runs[,c(2,3)])),]
# Total_DEGenes_Foley_Runs$Comparisons <- factor(as.factor(Total_DEGenes_Foley_Runs$Comparisons), levels = unique(as.factor(Total_DEGenes_Foley_Runs$Comparisons)) )
# Total_DEGenes_Foley_Runs
# 
# Total_DEGenes_Foley_Runs <- melt(Total_DEGenes_Foley_Runs, id.var="Comparisons")
# Total_DEGenes_Foley_Runs
# setwd("/Users/gordid/Desktop/Foley_Runs/")
# 
# if (SavePlots == TRUE) {pdf(file = paste0("/Users/gordid/Desktop/Foley_Runs/All_results", "/Total_DEGenes_P05_perComparison_BarGraph_Foley_Runs",
#                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".pdf"), width=30, height=12.5)}
# 
# ggplot(Total_DEGenes_Foley_Runs, aes(x = Comparisons, y = value, fill = variable)) +
#   geom_bar(stat = "identity") + xlab("Sample Comparisons") + ylab("DEG Counts (P < 0.05)") + coord_flip() + scale_fill_manual(values=c("red", "blue"))
# 
# if (SavePlots == TRUE) {dev.off()}
# 
# ##Bar plot of DEGenes (AdjP Value)
# Ranked_DEGenes_data_Foley_Runs <- list()
# Total_DEGenes_Foley_Runs <- data.frame()
# setwd("/Users/gordid/Desktop/Foley_Runs/Ranked_DEGenes/")
# for (i in 1:length(list.files())) {
#   Ranked_DEGenes_data_Foley_Runs[[i]] <- read.csv(list.files()[i])
#   print(list.files()[i])
#   Total_DEGenes_Foley_Runs[i,1] <- bar_plot_names[i]
#   Total_DEGenes_Foley_Runs[i,2] <- length(which(Ranked_DEGenes_data_Foley_Runs[[i]]$adj.P.Val[which(Ranked_DEGenes_data_Foley_Runs[[i]]$logFC > 0)] < 0.05))
#   Total_DEGenes_Foley_Runs[i,3] <- length(which(Ranked_DEGenes_data_Foley_Runs[[i]]$adj.P.Val[which(Ranked_DEGenes_data_Foley_Runs[[i]]$logFC < 0)] < 0.05))
# }
# names(Ranked_DEGenes_data_Foley_Runs) <- bar_plot_names
# colnames(Total_DEGenes_Foley_Runs) <- c("Comparisons", "DEG_Up_Counts", "DEG_Down_Counts")
# Total_DEGenes_Foley_Runs <- Total_DEGenes_Foley_Runs[order(rowSums(Total_DEGenes_Foley_Runs[,c(2,3)])),]
# Total_DEGenes_Foley_Runs$Comparisons <- factor(as.factor(Total_DEGenes_Foley_Runs$Comparisons), levels = unique(as.factor(Total_DEGenes_Foley_Runs$Comparisons)) )
# Total_DEGenes_Foley_Runs
# 
# Total_DEGenes_Foley_Runs <- melt(Total_DEGenes_Foley_Runs, id.var="Comparisons")
# Total_DEGenes_Foley_Runs
# setwd("/Users/gordid/Desktop/Foley_Runs/")
# 
# if (SavePlots == TRUE) {pdf(file = paste0("/Users/gordid/Desktop/Foley_Runs/All_results", "/Total_DEGenes_AdjP05_perComparison_BarGraph_Foley_Runs",
#                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".pdf"), width=30, height=12.5)}
# 
# ggplot(Total_DEGenes_Foley_Runs, aes(x = Comparisons, y = value, fill = variable)) +
#   geom_bar(stat = "identity") + xlab("Sample Comparisons") + ylab("DEG Counts (AdjP < 0.05)") + coord_flip() + scale_fill_manual(values=c("red", "blue"))
# 
# if (SavePlots == TRUE) {dev.off()}


####NetAct####



setwd("/Users/gordid/Desktop/Foley_Runs")

rownames(Combined_data_Foley_Runs_NetAct$counts) <- Combined_data_Foley_Runs_NetAct$genes$SYMBOL
phenoData_Foley_Runs = new("AnnotatedDataFrame", data = data.frame(celltype_Foley_Runs = group_name_Foley_Runs))
rownames(phenoData_Foley_Runs) = sample_name_Foley_Runs
head(Combined_data_Foley_Runs_NetAct$counts)
phenoData_Foley_Runs
colnames(Combined_data_Foley_Runs_NetAct$counts) <- sample_name_Foley_Runs
head(Combined_data_Foley_Runs_NetAct$counts)

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
for (i in 1:(length(levels(group_name_Foley_Runs))-1)) {
  for (j in (i+1):length(levels(group_name_Foley_Runs))) {
    compare_list <- append(compare_list, paste0(levels(group_name_Foley_Runs)[i], "-", levels(group_name_Foley_Runs)[j]))
    
  }
}
compare_list

if (TwoGroups == TRUE) {
  DErslt = RNAseqDegs_limma(Combined_data_Foley_Runs_NetAct$counts, phenoData_Foley_Runs, compare_list)
  str(DErslt)
  
  e = DErslt$e
  
  neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_Foley_Runs)
  
  save(neweset, DErslt, file = paste(Directory_CSVs_NetAct_RACIPE, "/macro_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                     "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     ".RData", sep = ""))
  
  write.csv(DErslt$table, file = paste(Directory_CSVs_NetAct_RACIPE, "/DEGenes_Foley_Runs_RACIPE_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                       "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                       ".csv", sep = ""))
  
  if (Calcs == TRUE) {
    gsearslt_Foley_Runs <- TF_GSEA(specimen, DErslt, minSize = 5, nperm = 10000, qval = T)
    write.csv(gsearslt_Foley_Runs, file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                ".csv", sep = ""))
  }
  setwd(Directory_CSVs_NetAct_RACIPE)
  gsearslt_Foley_Runs <- read.csv(file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                               "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                               ".csv", sep = ""))
  
  if (NetAct == TRUE) {
    qval <- 0.05
    Tfs_Foley_Runs <- as.character(gsearslt_Foley_Runs$tf[gsearslt_Foley_Runs$qvals < qval])
    tfs <- sort(unique(as.character(c(Tfs_Foley_Runs))))
    acts_mat = TF_Activity(tfs, specimen, neweset, DErslt, useDatabaseSign = FALSE)$all_activities
    tfs
  }
} else {
  DErslt = MultiRNAseqDegs_limma(Combined_data_Foley_Runs_NetAct$counts, phenoData_Foley_Runs, compare_list)
  str(DErslt)
  e = DErslt$Overall$e
  
  neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_Foley_Runs)
  
  save(neweset, DErslt, file = paste(Directory_CSVs_NetAct_RACIPE, "/macro_Foley_Runs", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".RData", sep = ""))
  
  NetAct_gene_list_Foley_Runs <- list()
  gsearslt_list <- list()
  tfs_list <- list()
  TFs_vector <- c()
  for (i in 1:length(compare_list)) {
    NetAct_gene_list_Foley_Runs[[i]] <- DErslt[[i]]$table
    write.csv(NetAct_gene_list_Foley_Runs[[i]], file = paste(Directory_CSVs_NetAct_RACIPE, "/DEGenes_Foley_Runs_RACIPE", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                             filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))
    
    if (Calcs == TRUE) {
      gsearslt_Foley_Runs_multi <- TF_GSEA(specimen, DErslt[[i]], minSize = 5, nperm = 10000, qval = T)
      write.csv(gsearslt_Foley_Runs_multi, file = paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                        filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),".csv", sep = ""))
    }
    setwd(Directory_CSVs_NetAct_RACIPE)
    gsearslt_list[[i]] <- read.csv(paste(Directory_CSVs_NetAct_RACIPE, "/Gsears_Foley_Runs", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                         filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE), ".csv", sep = ""))
    
    qval <- 0.05
    tfs_list[[i]] <- as.character(gsearslt_list[[i]]$tf[gsearslt_list[[i]]$qvals < qval])
    TFs_vector <- append(TFs_vector, tfs_list[[i]])
  }
  names(NetAct_gene_list_Foley_Runs) <- Comparisons_vector
  names(gsearslt_list) <- Comparisons_vector
  names(tfs_list) <- Comparisons_vector
  
  tfs <- sort(unique(as.character(TFs_vector)))
  lengths(tfs_list)
  length(tfs)
  tfs
  acts_mat = TF_Activity(tfs, specimen, neweset, DErslt$Overall)$all_activities
}

####ExpressionHeatmap####
if (RACIPE == TRUE) {
  
  setwd("/Users/gordid/Desktop/Foley_Runs")
  if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_NetAct_RACIPE, "/TFsHeatMap_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), "_qval", 
                                           gsub(".", "", as.character(qval), fixed = TRUE),
                                           ".pdf", sep = ""), width = 15, height = 15)}
  
  Combine_heatmap(acts_mat, neweset) #Heatmap
  
  if (SavePlots == TRUE) {dev.off()}
  
  
  miTh = 0.6
  nbins = 3
  useCor = TRUE
  tf_links <- TF_Filter(acts_mat, specimen, miTh = miTh, nbins = nbins, useCor = useCor, DPI = FALSE, removeSignalling = FALSE)
  dim(tf_links)
  plot_network_v(tf_links)
  
  ####Generating Network####
  if (SavePlots == TRUE) {
    setwd(Directory_Plots_NetAct_RACIPE)
    network_plot_Foley_Runs <- plot_network_v(tf_links)
    htmlwidgets::saveWidget(network_plot_Foley_Runs,file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs_html_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                                 "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                                 "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                                                 if (useCor == TRUE) {"_useCor" }, 
                                                                 ".html", sep = ""))
    
    webshot::webshot(paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs_html_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                           "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                           if (useCor == TRUE) {"_useCor" }, 
                           ".html", sep = ""),
                     file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), 
                                  "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                  "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins), 
                                  if (useCor == TRUE) {"_useCor" }, 
                                  ".png", sep = ""), delay = 10,  vwidth = 1500, vheight = 1500)
    
    
  } else {
    plot_network_v(tf_links)
  }
  
  x <- union(tf_links$Source,tf_links$Target)
  x
  
  
  #if (SavePlots == TRUE) {
  # orca(plot_network_v(tf_links), file = paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Foley_Runs", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
  #    "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
  #                              if (outliers == TRUE) {print("_outliers-removed")}, ".png", sep = "")) 
  #} else {
  #  plot_network_v(tf_links)
  #}
  
  write.csv(tf_links, file = paste(Directory_CSVs_NetAct_RACIPE, "/TFNetwork_InteractionsTable_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                   "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                   if (useCor == TRUE) {"_useCor" }, 
                                   ".csv", sep = ""))
  
  #write.table(tf_links, "net_emt.tpo", quote = FALSE, row.names = FALSE)
  write.table(tf_links, file = paste(Directory_Plots_NetAct_RACIPE, "/net_emt_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), 
                                     "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                     "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                     if (useCor == TRUE) {"_useCor" }, 
                                     ".tpo", sep = ""), quote = FALSE, row.names = FALSE)
  
  ####PERTURBATION####
  if (SavePlots == TRUE) {
    setwd(Directory_Plots_NetAct_RACIPE)
    library(sRACIPE)
    
    racipe <- sRACIPE::simulateGRC(paste(Directory_Plots_NetAct_RACIPE, "/net_emt_Foley_Runs_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                         "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                         "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE), "_nbins", as.character(nbins),
                                         if (useCor == TRUE) {"_useCor" }, 
                                         ".tpo", sep = ""), quote = FALSE, row.names = FALSE, plots = TRUE, plotToFile = TRUE)
    racipePer <- sRACIPE::knockdownAnalysis(racipe,reduceProduction = 10)
    setwd("/Users/gordid/Desktop/Foley_Runs")
  } else {
    print("No Perturbation")
  }
}

#gsearslt_1 = read.csv(paste("Foley_Runs_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = read.csv(paste("Foley_Runs_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = read.csv(paste("Foley_Runs_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#tfs_1 = as.character(gsearslt_1$tf[gsearslt_1$qvals <0.2])
#tfs_2 = as.character(gsearslt_2$tf[gsearslt_2$qvals <0.2])
#tfs_3 = as.character(gsearslt_3$tf[gsearslt_3$qvals <0.2])
#tfs = sort(unique(as.character(c(tfs_1, tfs_2, tfs_3))))
#lengths(list(tfs_1,tfs_2,tfs_3, tfs))
#tfs
# for normal use ####
#e = DErslt$Overall$e
#neweset = ExpressionSet(assayData = as.matrix(e), phenoData = phenoData_Foley_Runs)
#save(neweset, DErslt, file = "macro.RData")

#Lupus_Foley_Runs_BvC_genes <- (DErslt$`B-C`$table)
#head(Lupus_Foley_Runs_BvC_genes)

#write.csv(Lupus_Foley_Runs_BvC_genes, file = paste("Foley_Runs_BvC_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#Lupus_Foley_Runs_AvC_genes <- (DErslt$`A-C`$table)
#write.csv(Lupus_Foley_Runs_BvC_genes, file = paste("Foley_Runs_AvC_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#Lupus_Foley_Runs_AvB_genes <- (DErslt$`A-B`$table)
#write.csv(Lupus_Foley_Runs_BvC_genes, file = paste("Foley_Runs_AvB_genes_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))


#gsearslt_1 = TF_GSEA(mgs, DErslt$`B-C`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_1, file = paste("Foley_Runs_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = TF_GSEA(mgs, DErslt$`A-B`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_2, file = paste("Foley_Runs_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = TF_GSEA(mgs, DErslt$`A-C`, minSize=5, nperm = 10000, qval = T)
#write.csv(gsearslt_3, file = paste("Foley_Runs_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))

#gsearslt_1 = read.csv(paste("Foley_Runs_BvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_2 = read.csv(paste("Foley_Runs_AvB_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#gsearslt_3 = read.csv(paste("Foley_Runs_AvC_gsears_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
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
#write.table(tf_links, paste("Foley_Runs_net_emt_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
#           quote = FALSE, row.names = FALSE)

####PERTURBATION
#library(sRACIPE)
#racipe <- sRACIPE::simulateGRC(paste("Foley_Runs_net_emt_", filter_type, "lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".tpo", sep = ""),
#                              quote = FALSE, row.names = FALSE, plots = TRUE, plotToFile = TRUE)

#racipePer <- sRACIPE::knockdownAnalysis(racipe,reduceProduction = 10)


