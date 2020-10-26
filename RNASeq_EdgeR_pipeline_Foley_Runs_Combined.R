#read in data file on experiment
setwd("/Users/gordid/Desktop/Foley_Runs_Combined/")

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
Rename_file_extensions("/Users/gordid/Desktop/Foley_Runs_Combined/genes_results/")


#reading in design CSV
experiment_Foley_Runs_Combined <- read.csv(paste0("/Users/gordid/Desktop/Foley_Runs_Combined/", "Foley_Runs_Combined_experimental_design", ".csv"))
experiment_Foley_Runs_Combined

#creating vector of data files
files_Foley_Runs_Combined_unsorted <- list.files("/Users/gordid/Desktop/Foley_Runs_Combined/genes_results/")
files_Foley_Runs_Combined_unsorted

#Putting data files character vector into the order of the sample IDS in the design experiment CSV
files_Foley_Runs_Combined <- c()

for (i in 1:length(files_Foley_Runs_Combined_unsorted)) {
  for (j in 1:length(experiment_Foley_Runs_Combined$Sample.Name)) {
    if (grepl(as.character(experiment_Foley_Runs_Combined$Sample.Name)[j], files_Foley_Runs_Combined_unsorted[i]) == TRUE) {
      files_Foley_Runs_Combined[j] <- files_Foley_Runs_Combined_unsorted[i]
    }
  }
}
files_Foley_Runs_Combined
#Files need to be in the same order as the samples in the excel file. Find a way to call by name instead of index.

####Save plots?####
SavePlots = FALSE
SaveCSVs = FALSE
EnrichR = FALSE
FGSEA = FALSE

####What Organism?####
Mouse = TRUE 

####Comparing just 2 groups?####
TwoGroups = TRUE

if (TwoGroups == TRUE) {
  
  Group1 = "CD6weeksSedH"
  Group2 = "WD6weekSedH"
  
  WhichGroups = c(which(experiment_Foley_Runs_Combined$Group.Name == Group1),which(experiment_Foley_Runs_Combined$Group.Name == Group2))
  
  experiment_Foley_Runs_Combined <- experiment_Foley_Runs_Combined[WhichGroups,]
  
  files_Foley_Runs_Combined <- files_Foley_Runs_Combined[WhichGroups]
}


####Are there Outliers?####
outliers = TRUE

if (outliers == TRUE) {
  Outlier1 = "17979H"
  Outlier2 = "17976H"
  
  outliers_index <- c(which(experiment_Foley_Runs_Combined$Sample.Name == Outlier1),which(experiment_Foley_Runs_Combined$Sample.Name == Outlier2))
  
  outliers_name <- as.character(experiment_Foley_Runs_Combined$Sample.Name[outliers_index])
  
  experiment_Foley_Runs_Combined <- experiment_Foley_Runs_Combined[-outliers_index,]
  
  files_Foley_Runs_Combined <- files_Foley_Runs_Combined[-outliers_index]
  
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
setwd("/Users/gordid/Desktop/Foley_Runs_Combined/genes_results/")
file_1_Foley_Runs_Combined <- read.delim(files_Foley_Runs_Combined[1], header=TRUE)
head(file_1_Foley_Runs_Combined)
str(file_1_Foley_Runs_Combined)
dim(file_1_Foley_Runs_Combined)
setwd("/Users/gordid/Desktop/Foley_Runs_Combined")

####Renaming columns and creating groups####

group_name_Foley_Runs_Combined <- as.character(experiment_Foley_Runs_Combined$Group.Name)
group_name_Foley_Runs_Combined <- factor(as.factor(group_name_Foley_Runs_Combined), levels = unique(as.factor(group_name_Foley_Runs_Combined)) ) #preserves order of factors
group_name_Foley_Runs_Combined
levels(group_name_Foley_Runs_Combined)


sample_name_Foley_Runs_Combined <- as.character(experiment_Foley_Runs_Combined$Sample.Name)
sample_name_Foley_Runs_Combined

condition_Foley_Runs_Combined <- factor(as.factor(as.character(experiment_Foley_Runs_Combined$Condition)), levels = unique(as.factor(as.character(experiment_Foley_Runs_Combined$Condition))) )
condition_Foley_Runs_Combined
length(levels(condition_Foley_Runs_Combined))

diet_Foley_Runs_Combined <- factor(as.factor(as.character(experiment_Foley_Runs_Combined$Diet)), levels = unique(as.factor(as.character(experiment_Foley_Runs_Combined$Diet))) )
diet_Foley_Runs_Combined

diet_duration_Foley_Runs_Combined <- factor(as.factor(as.character(experiment_Foley_Runs_Combined$Diet.Duration)), levels = unique(as.factor(as.character(experiment_Foley_Runs_Combined$Diet.Duration))) )
diet_duration_Foley_Runs_Combined

####creating vector of gene names and vector of possible sample comparisons####

Comparisons_vector <- c()
for (i in 1:(length(levels(group_name_Foley_Runs_Combined))-1)) {
  for (j in (i+1):length(levels(group_name_Foley_Runs_Combined))) {
    
    Comparisons_vector <- append(Comparisons_vector, paste0(levels(group_name_Foley_Runs_Combined)[i], "vs", levels(group_name_Foley_Runs_Combined)[j]))
    
  }
}

Comparisons_vector
length(Comparisons_vector)

#Creating new directorys####

#plots/RNASeq folder
Directory_Plots_RNASeq <- paste0("/Users/gordid/Desktop/Foley_Runs_Combined/All_results/Plots/RNASeq/",
                                 (if(TwoGroups == TRUE) { 
                                   paste0(levels(group_name_Foley_Runs_Combined)[1], "vs", levels(group_name_Foley_Runs_Combined)[2],
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

Directory_Plots_NetAct_RACIPE <- paste0("/Users/gordid/Desktop/Foley_Runs_Combined/All_results/Plots/NetAct_RACIPE/",
                                        (if(TwoGroups == TRUE) { 
                                          paste0(levels(group_name_Foley_Runs_Combined)[1], "vs", levels(group_name_Foley_Runs_Combined)[2],
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
Directory_CSVs_DEGenes <- paste0("/Users/gordid/Desktop/Foley_Runs_Combined/All_results/CSVs/DEGenes/", 
                                 (if(TwoGroups == TRUE) {
                                   paste0(levels(group_name_Foley_Runs_Combined)[1], "vs", levels(group_name_Foley_Runs_Combined)[2],
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
Directory_CSVs_NetAct_RACIPE <- paste0("/Users/gordid/Desktop/Foley_Runs_Combined/All_results/CSVs/NetAct_RACIPE/",
                                       (if(TwoGroups == TRUE) { 
                                         paste0(levels(group_name_Foley_Runs_Combined)[1], "vs", levels(group_name_Foley_Runs_Combined)[2],
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
setwd("/Users/gordid/Desktop/Foley_Runs_Combined/genes_results/")
combined_data_Foley_Runs_Combined <- readDGE(files_Foley_Runs_Combined, columns = c(1, 5))
dim(combined_data_Foley_Runs_Combined) #showing dimensions of "counts" matrix
class(combined_data_Foley_Runs_Combined) #showing class of the whole object
summary(combined_data_Foley_Runs_Combined) #showing names, size, class, and mode of each element of the object, not sure why class is "none" for counts
str(combined_data_Foley_Runs_Combined) #showing the structure of the object
head(combined_data_Foley_Runs_Combined$counts)
colnames(combined_data_Foley_Runs_Combined) <- sample_name_Foley_Runs_Combined
setwd("/Users/gordid/Desktop/Foley_Runs_Combined")

####Assigning experimental info to categories in the DGEList obect####
combined_data_Foley_Runs_Combined$samples$group <- group_name_Foley_Runs_Combined
combined_data_Foley_Runs_Combined$samples

#combined_data_Foley_Runs_Combined$samples$remarks <- remarks_Foley_Runs_Combined

####Create new genes matrix by retrieving gene information from the mus library. The gene information of their SYMBOL and TXCHROM#### 
#number of gene IDs
gene_id_Foley_Runs_Combined <- rownames(combined_data_Foley_Runs_Combined)
head(gene_id_Foley_Runs_Combined)
length(gene_id_Foley_Runs_Combined)
#number of duplicated ENSEMBL IDS PRE select
dup_genes_Foley_Runs_Combined_PRESelect <- which(duplicated(gene_id_Foley_Runs_Combined))
length(dup_genes_Foley_Runs_Combined_PRESelect)

columns(Mus.musculus)
if (Mouse == TRUE) {
  genes_Foley_Runs_Combined <- biomaRt::select(Mus.musculus, keys = gene_id_Foley_Runs_Combined, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
} else {
  genes_Foley_Runs_Combined <- biomaRt::select(Homo.sapiens, keys = gene_id_Foley_Runs_Combined, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
  
}
dim(genes_Foley_Runs_Combined)
head(genes_Foley_Runs_Combined)
length(genes_Foley_Runs_Combined$ENSEMBL) - length(gene_id_Foley_Runs_Combined)

####symbols are duplicated too often and so I dont remove genes based on duplicated symbols (cause of NA mostly, will remove NAs later)####
head((which(duplicated(genes_Foley_Runs_Combined$SYMBOL))))
length(which(duplicated(genes_Foley_Runs_Combined$SYMBOL)))
genes_Foley_Runs_Combined[head((which(duplicated(genes_Foley_Runs_Combined$SYMBOL)))),]

#number of duplicated ENSEMBL IDS POST select
dup_genes_Foley_Runs_Combined <- which(duplicated(genes_Foley_Runs_Combined$ENSEMBL))
length(dup_genes_Foley_Runs_Combined)
head(dup_genes_Foley_Runs_Combined)
genes_Foley_Runs_Combined[c(head(which(duplicated(genes_Foley_Runs_Combined$ENSEMBL))), head(which(duplicated(genes_Foley_Runs_Combined$ENSEMBL))) - 1),]

####remove duplicated genes based on duplicated ENSEMBL ID####
genes_Foley_Runs_Combined_final <- genes_Foley_Runs_Combined[-dup_genes_Foley_Runs_Combined, ]
dim(genes_Foley_Runs_Combined_final)
dim(combined_data_Foley_Runs_Combined$counts)

#Number of removed Symbol IDS must be equal to or less than the number of removed ENSEMBL IDs because during the select process, some ENSEMBL ID's may match to two different Symbol IDs.
length(which(duplicated(genes_Foley_Runs_Combined$SYMBOL))) - length(which(duplicated(genes_Foley_Runs_Combined_final$SYMBOL)))
length(dup_genes_Foley_Runs_Combined)

#assign data frame of gene annotations to the combined_data object
combined_data_Foley_Runs_Combined$genes <- genes_Foley_Runs_Combined_final 
summary(combined_data_Foley_Runs_Combined)

#RAW COUNTS DATA
if (SaveCSVs == TRUE) {
  
  write.csv(cbind2(combined_data_Foley_Runs_Combined$counts, combined_data_Foley_Runs_Combined$genes), file = paste0(Directory_CSVs_DEGenes,
                                                                                                   "/Raw_Counts_Data_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), ".csv"))
  
}

####CPM normalization on count data####
cpm_Foley_Runs_Combined <- cpm(combined_data_Foley_Runs_Combined)
lcpm_Foley_Runs_Combined <-cpm(combined_data_Foley_Runs_Combined, log = TRUE)

L_Foley_Runs_Combined <- mean(combined_data_Foley_Runs_Combined$samples$lib.size) * 1e-6 #calculation of L for the prior count.
M_Foley_Runs_Combined <- median(combined_data_Foley_Runs_Combined$samples$lib.size) * 1e-6
c(L_Foley_Runs_Combined, M_Foley_Runs_Combined)
head(lcpm_Foley_Runs_Combined)

table(rowSums(combined_data_Foley_Runs_Combined$counts==0)==length(colnames(lcpm_Foley_Runs_Combined)))

###FILTER THE GENES####
if (CPMFilter == TRUE) {
  kept_expression_Foley_Runs_Combined <- filterByExpr(combined_data_Foley_Runs_Combined, group = combined_data_Foley_Runs_Combined$samples$group) #filterByExpr function
  # keeps genes with 10 read counts or more in a minimum number of samples, where the minimum number of samples is chosen according 
  #to the minimum group sample size. 
}  else {
  #filter by counts, rather than by cpm
  kept_expression_Foley_Runs_Combined <- (rowSums(combined_data_Foley_Runs_Combined$counts > 10) >= 1)
}
length(kept_expression_Foley_Runs_Combined)
head(kept_expression_Foley_Runs_Combined)

#removing all the lowly expressed genes from the combined_data object. When you subset a DGEList and specify keep.lib.sizes=FALSE,
#the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
filtered_combined_data_Foley_Runs_Combined <- combined_data_Foley_Runs_Combined[kept_expression_Foley_Runs_Combined, , keep.lib.sizes = FALSE] 
dim(filtered_combined_data_Foley_Runs_Combined)
dim(combined_data_Foley_Runs_Combined)

####Removing duplicated gene symbols for the NetAct analysis####
Combined_data_Foley_Runs_Combined_NetAct <- filtered_combined_data_Foley_Runs_Combined
rownames(Combined_data_Foley_Runs_Combined_NetAct$counts) <- Combined_data_Foley_Runs_Combined_NetAct$genes$SYMBOL
dim(Combined_data_Foley_Runs_Combined_NetAct)
#Number of Symbols that are NA
dup_genes_Foley_Runs_Combined_NA <- which(is.na(filtered_combined_data_Foley_Runs_Combined$genes$SYMBOL))
length(dup_genes_Foley_Runs_Combined_NA)
#Number of duplicated Symbols that are not NA
length(which(duplicated(na.omit(filtered_combined_data_Foley_Runs_Combined$genes$SYMBOL))))
#Number of duplicated symbols + all NAs
length(dup_genes_Foley_Runs_Combined_NA) + length(which(duplicated(na.omit(filtered_combined_data_Foley_Runs_Combined$genes$SYMBOL))))
# of genes theoretically left over at removing of all duplicated
length(filtered_combined_data_Foley_Runs_Combined$genes$SYMBOL) -(length(dup_genes_Foley_Runs_Combined_NA) + length(which(duplicated(na.omit(filtered_combined_data_Foley_Runs_Combined$genes$SYMBOL)))))

#removing all NAs
Combined_data_Foley_Runs_Combined_NetAct <-Combined_data_Foley_Runs_Combined_NetAct[-dup_genes_Foley_Runs_Combined_NA, ]
dim(Combined_data_Foley_Runs_Combined_NetAct)

#removing all duplicated symbols
dup_genes_Foley_Runs_Combined_SYMBOL <- which(duplicated(Combined_data_Foley_Runs_Combined_NetAct$genes$SYMBOL))
length(dup_genes_Foley_Runs_Combined_SYMBOL)
Combined_data_Foley_Runs_Combined_NetAct <- Combined_data_Foley_Runs_Combined_NetAct[-dup_genes_Foley_Runs_Combined_SYMBOL, ]
dim(Combined_data_Foley_Runs_Combined_NetAct)

####Plotting filtered Data####
lcpm.cutoff_Foley_Runs_Combined <- log2(10/M_Foley_Runs_Combined + 2/L_Foley_Runs_Combined)
lcpm.cutoff_Foley_Runs_Combined
nsamples <- ncol(combined_data_Foley_Runs_Combined$counts)
col <- brewer.pal(nsamples, "Paired")

if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/FilteredData_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                  ".pdf", sep = ""))}
par(mfrow = c(1,2))
plot(density(lcpm_Foley_Runs_Combined[, 1]), col = col[1], lwd=2, ylim=c(0,0.46), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_Foley_Runs_Combined, lty=3)
title(main="A. Raw data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_Foley_Runs_Combined[, i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

combined_data_Foley_Runs_Combined <- filtered_combined_data_Foley_Runs_Combined
lcpm_Foley_Runs_Combined <- cpm(combined_data_Foley_Runs_Combined, log = TRUE)

plot(density(lcpm_Foley_Runs_Combined[, 1]), col = col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_Foley_Runs_Combined, lty=3)
title(main="B. Filtered data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_Foley_Runs_Combined[, i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
if (SavePlots == TRUE) {dev.off() }


####Calculate normalization factors and creates new DGEList object with these new normalization factors inserted. THEN GENERATING BOX PLOT####
colnames(lcpm_Foley_Runs_Combined) <- sample_name_Foley_Runs_Combined
combined_data_Foley_Runs_Combined$samples$norm.factors

if (SavePlots == TRUE) {pdf(paste0(Directory_Plots_RNASeq, "/TMM_Normalized_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), ".pdf"))}


par(mfrow = c(1,2))
boxplot(lcpm_Foley_Runs_Combined, las = 2, col = col, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")

combined_data_Foley_Runs_Combined = calcNormFactors(combined_data_Foley_Runs_Combined, method = "TMM") #calculates normalization factors and creates new DGEList object
#these normalization factors are important because if there is a sample where expression level is different accross the board, that may have been caused
#by an eternal factor that is not of biological interest. TMM helps normalize for that, using the assumption that should have a similar range
#and distribution of expression values.
lcpm_Foley_Runs_Combined <- cpm(combined_data_Foley_Runs_Combined, log = TRUE) 
colnames(lcpm_Foley_Runs_Combined) <- sample_name_Foley_Runs_Combined
head(lcpm_Foley_Runs_Combined)

combined_data_Foley_Runs_Combined$samples$norm.factors
boxplot(lcpm_Foley_Runs_Combined, las = 2, col = col, main="")
title(main="B.Normalised with TMM data", ylab="Log-cpm")
if (SavePlots == TRUE) {dev.off() }

####PCA and PCoA Plots####
if (SavePlots == TRUE) {pdf(paste0(Directory_Plots_RNASeq, "/PCA_bygroup_1vs2_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                   ".pdf"))}
color_groups_Foley_Runs_Combined <- combined_data_Foley_Runs_Combined$samples$group
levels(color_groups_Foley_Runs_Combined) <- brewer.pal(nlevels(color_groups_Foley_Runs_Combined), "Set1")
color_groups_Foley_Runs_Combined <- as.character(color_groups_Foley_Runs_Combined)
color_groups_Foley_Runs_Combined

combined_data_Foley_Runs_Combined$samples$names <- sample_name_Foley_Runs_Combined
combined_data_Foley_Runs_Combined$samples

par(mfrow = c(1,1))
#plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$names, col = color_groups_Foley_Runs_Combined, main = "PCoA")
plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$group, col = color_groups_Foley_Runs_Combined, gene.selection="common", 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))
if (SavePlots == TRUE) {dev.off() }

par(mfrow = c(1,2))
plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$names, col = color_groups_Foley_Runs_Combined, dim = c(2,3), 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCoA")) 
plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$names, col = color_groups_Foley_Runs_Combined, gene.selection="common", dim = c(2,3),
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))

par(mfrow = c(1,1))
plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$names, col = color_groups_Foley_Runs_Combined, gene.selection="common", 
        main = paste0(gsub("_[^_]+$" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"))


glMDSPlot(lcpm_Foley_Runs_Combined, labels=paste(combined_data_Foley_Runs_Combined$samples$names, combined_data_Foley_Runs_Combined$samples$group, sep="_"), 
          groups = combined_data_Foley_Runs_Combined$samples[,c(2,5)], 
          gene.selection="common", main = paste0(sub("_[^_]+$*" ,"", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")), " ", "PCA"), launch=TRUE, path = Directory_Plots_RNASeq,
          html = paste0("InteractivePCA_ScreePlot_1vs2_Foley_Runs_Combined", "_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse="")))

#scree.plot(plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$names, col = color_groups_Foley_Runs_Combined, gene.selection="common", main = "PCA"))
#a <- (plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$names, col = color_groups_Foley_Runs_Combined, main = "PCoA"))
#screeplot(a, xlim = c(-8,13), ylim=c(0,0.46))
####PCA Plots for other phenotypic data and meta-data####
#color_age_Foley_Runs_Combined <- combined_data_Foley_Runs_Combined$samples$age
#levels(color_age_Foley_Runs_Combined) <- brewer.pal(nlevels(color_age_Foley_Runs_Combined), "Set2")
#color_age_Foley_Runs_Combined <- as.character(color_age_Foley_Runs_Combined)
#color_age_Foley_Runs_Combined
#plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$age, col = color_groups_Foley_Runs_Combined)
#plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$age, col = color_groups_Foley_Runs_Combined, gene.selection="common")

#color_x2dg_Foley_Runs_Combined <- combined_data_Foley_Runs_Combined$samples$x2dg
#levels(color_x2dg_Foley_Runs_Combined) <- brewer.pal(nlevels(color_x2dg_Foley_Runs_Combined), "Set3")
#color_x2dg_Foley_Runs_Combined <- as.character(color_x2dg_Foley_Runs_Combined)
#color_x2dg_Foley_Runs_Combined
#plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$x2dg, col = color_groups_Foley_Runs_Combined)
#plotMDS(lcpm_Foley_Runs_Combined, labels = combined_data_Foley_Runs_Combined$samples$x2dg, col = color_groups_Foley_Runs_Combined, gene.selection="common")

####Building Design Matrix####

design_matrix_Foley_Runs_Combined_no.int <- model.matrix(~0+combined_data_Foley_Runs_Combined$samples$group) #no intercept
colnames(design_matrix_Foley_Runs_Combined_no.int) = gsub("combined_data_Foley_Runs_Combined$samples$group", "", colnames(design_matrix_Foley_Runs_Combined_no.int), fixed = TRUE)
design_matrix_Foley_Runs_Combined_no.int

design_matrix_Foley_Runs_Combined_int <- model.matrix(~combined_data_Foley_Runs_Combined$samples$group) #with intercept
colnames(design_matrix_Foley_Runs_Combined_int) = gsub("combined_data_Foley_Runs_Combined$samples$group", "", colnames(design_matrix_Foley_Runs_Combined_int), fixed = TRUE)
design_matrix_Foley_Runs_Combined_int


if (TwoGroups == TRUE) {
  contr.matrix_no.int <- makeContrasts(
    paste(levels(group_name_Foley_Runs_Combined)[1]," - ", levels(group_name_Foley_Runs_Combined)[2], sep = ""),
    levels = colnames(design_matrix_Foley_Runs_Combined_no.int))
  colnames(contr.matrix_no.int) <- paste(levels(group_name_Foley_Runs_Combined)[1],"vs", levels(group_name_Foley_Runs_Combined)[2], sep = "")
  contr.matrix_no.int
} else {
  contrast_vector <- c()
  column_name_vector <- c()
  for (i in 1:(length(levels(group_name_Foley_Runs_Combined))-1)) {
    
    for (j in (i+1):length(levels(group_name_Foley_Runs_Combined))) {
      
      contrast_vector <- append(contrast_vector, paste(levels(group_name_Foley_Runs_Combined)[i]," - ", levels(group_name_Foley_Runs_Combined)[j], sep = ""))
      column_name_vector <- append(column_name_vector, paste(levels(group_name_Foley_Runs_Combined)[i], "vs", levels(group_name_Foley_Runs_Combined)[j], sep = ""))
    }
    
  }
  
  contr.matrix_no.int <- makeContrasts(
    contrasts = contrast_vector,
    levels = colnames(design_matrix_Foley_Runs_Combined_no.int)
  )
  colnames(contr.matrix_no.int) <- column_name_vector
  contr.matrix_no.int
}


#contr.matrix_no.int <- makeContrasts(
#AvsB = "A- B",
#AvsC = A - C,
#BvsC = B - C,
#levels = colnames(design_matrix_Foley_Runs_Combined_no.int))
#contr.matrix_no.int
#colnames(contr.matrix_no.int) <- c("a", "b", "c")

#contr.matrix_no.int <- makeContrasts(
#contrasts = contrast_vector,
#levels = colnames(design_matrix_Foley_Runs_Combined_no.int))
#contr.matrix_no.int

#contr.matrix_int <- makeContrasts( #Intercept column has replcaced "A", not sure how contrast matrix works now
# AvsB = A - B,
# AvsC = A - C,
# BvsC = B - C,
#levels = colnames(design_matrix_Foley_Runs_Combined_int))
#contr.matrix_int

####VOOM####
adjP = 0.05
if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Voom_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                  ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,2))
v_Foley_Runs_Combined <- voom(combined_data_Foley_Runs_Combined, design_matrix_Foley_Runs_Combined_no.int, plot=TRUE)
str(v_Foley_Runs_Combined)
dim(v_Foley_Runs_Combined$genes)

vfit_Foley_Runs_Combined <- lmFit(v_Foley_Runs_Combined, design_matrix_Foley_Runs_Combined_no.int)
summary(decideTests(vfit_Foley_Runs_Combined))
str(vfit_Foley_Runs_Combined)

vfit_Foley_Runs_Combined_contrasts <- contrasts.fit(vfit_Foley_Runs_Combined, contrasts=contr.matrix_no.int) 
summary(decideTests(vfit_Foley_Runs_Combined_contrasts))

efit_Foley_Runs_Combined <- eBayes(vfit_Foley_Runs_Combined_contrasts)
summary(decideTests(efit_Foley_Runs_Combined))
de_Foley_Runs_Combined <- decideTests(efit_Foley_Runs_Combined)

plotSA(efit_Foley_Runs_Combined)
title(main="voom: Mean-variance trend-Normalized")

if (SavePlots == TRUE) {dev.off()}


####Filter further with logfold change and Venn Diagram####
#Figure out the criteria data has to choose appropriate lfc
lfc = 0.1
tfit_Foley_Runs_Combined <- treat(vfit_Foley_Runs_Combined_contrasts, lfc=lfc)
dt_Foley_Runs_Combined <- decideTests(tfit_Foley_Runs_Combined)
summary(dt_Foley_Runs_Combined)


topTable_p <-  topTable(efit_Foley_Runs_Combined, n=Inf); head(topTable_p)
topTreat_p <- topTreat(efit_Foley_Runs_Combined, n=Inf); head(topTreat_p)
topTreat_lfc <- topTreat(tfit_Foley_Runs_Combined, n=Inf); head(topTreat_lfc)

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
  de.common_Foley_Runs_Combined_2_comparisons <- which(dt_Foley_Runs_Combined[,1]!=0 & dt_Foley_Runs_Combined[,2]!=0)
  print(length(de.common_Foley_Runs_Combined_2_comparisons))
  
  de.common_Foley_Runs_Combined_3_comparisons <- which(dt_Foley_Runs_Combined[,1]!=0 & dt_Foley_Runs_Combined[,2]!=0 & dt_Foley_Runs_Combined[,3]!=0)
  print(length(de.common_Foley_Runs_Combined_3_comparisons))
  
  if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Venn_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                    ".pdf", sep = ""), width=13, height=10)}
  par(mfrow = c(1,1))
  vennDiagram(dt_Foley_Runs_Combined[,1:3], circle.col=c("turquoise", "salmon", "orange"))
  if (SavePlots == TRUE) {dev.off()}
}

if (SaveCSVs == TRUE) {
  write.fit(tfit_Foley_Runs_Combined, dt_Foley_Runs_Combined, file = paste(Directory_CSVs_DEGenes, "/tfit_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                         "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".txt", sep = ""))
  
  tfit_Foley_Runs_Combined.txt = read.table(file = paste(Directory_CSVs_DEGenes, "/tfit_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                                gsub(".", "", as.character(lfc), fixed = TRUE),
                                                ".txt", sep = ""))
  
  write.csv(tfit_Foley_Runs_Combined.txt, file = paste(Directory_CSVs_DEGenes, "/tfit_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                              gsub(".", "", as.character(lfc), fixed = TRUE),
                                              ".csv", sep = ""))
}

####Mean-difference plots and arranging DE genes####
if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/Mean-Difference-Plots_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), "_lfc", 
                                         gsub(".", "", as.character(lfc), fixed = TRUE),
                                         ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,length(Comparisons_vector)))
for (i in 1:length(Comparisons_vector)) {
  plotMD(efit_Foley_Runs_Combined, column = i, status = de_Foley_Runs_Combined[,i], xlim = c(-8,13), main = paste0(Comparisons_vector[i], " Mean-Difference Plot"))
  
  glMDPlot(efit_Foley_Runs_Combined, coef=i, status=de_Foley_Runs_Combined, main=colnames(efit_Foley_Runs_Combined)[i], side.main="ENSEMBL", counts=lcpm_Foley_Runs_Combined,
           groups=combined_data_Foley_Runs_Combined$samples$group,launch=FALSE, 
           path = Directory_Plots_RNASeq, html = paste0("Interactive_MeanDifference_Plot_Foley_Runs_Combined", "_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""), 
                                                        if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE)))
  
}
if (SavePlots == TRUE){ dev.off()}

#par(mfrow = c(1,length(levels(group_name_Foley_Runs_Combined))))
#Ranked_DEGenes_list <- list()
#for (i in 1:(length(levels(group_name_Foley_Runs_Combined))-1)) {
#  for (j in (i+1):length(levels(group_name_Foley_Runs_Combined))) {
#    print(paste0(levels(group_name_Foley_Runs_Combined)[i], "vs", levels(group_name_Foley_Runs_Combined)[j]))
#    assign(paste0(levels(group_name_Foley_Runs_Combined)[i], "vs", levels(group_name_Foley_Runs_Combined)[j]), (topTreat(tfit_Foley_Runs_Combined, coef = i, n = Inf)))
#    head(eval(as.name(paste0(levels(group_name_Foley_Runs_Combined)[i], "vs", levels(group_name_Foley_Runs_Combined)[j]))))
#  }
#  
#} 
####creating list of Toptreat-DEGenes####
# par(mfrow = c(1,length(levels(group_name_Foley_Runs_Combined))))
Ranked_DEGenes_list <- list()
Ranked_DEGenes_list_tfit <- list()
k = 0
for (i in 1:(length(levels(group_name_Foley_Runs_Combined))-1)) {
  for (j in (i+1):length(levels(group_name_Foley_Runs_Combined))) {
    k = k+1 
    Ranked_DEGenes_list[[k]] <- topTable(efit_Foley_Runs_Combined, coef = k, n = Inf)
    Ranked_DEGenes_list_tfit[[k]] <- topTreat(tfit_Foley_Runs_Combined, coef = k, n = Inf)
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
    write.csv(Ranked_DEGenes_list[[i]], file = paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                      if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
    write.csv(Ranked_DEGenes_list_tfit[[i]], file = paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
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
      
      fwrite(fgsea_results, file = paste0(paste0(Directory_CSVs_DEGenes, "/FGSEA", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), "/FGSEA_Foley_Runs_Combined_", 
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
        
        fwrite(fgsea_gskb_results, file = paste0(paste0(Directory_CSVs_DEGenes, "/FGSEA_GSKB", if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}), "/FGSEA_GSKB_Foley_Runs_Combined_", 
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

    Foley_Runs_Combined_enrichR <- read.csv(paste0(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                          if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
    
    if (EnrichR == TRUE) {        
      Foley_Runs_Combined_enrichR_top800 <- enrichr(as.character(Foley_Runs_Combined_enrichR$SYMBOL[1:800]), enrichR_dbs)
      head(Foley_Runs_Combined_enrichR_top800[3])
      
      Foley_Runs_Combined_enrichR_Down <- Foley_Runs_Combined_enrichR[which(Foley_Runs_Combined_enrichR$logFC < 0), ]
      Foley_Runs_Combined_enrichR_Down_0.05 <- Foley_Runs_Combined_enrichR_Down[which(Foley_Runs_Combined_enrichR_Down$adj.P.Val < 0.05), ]
      dim(Foley_Runs_Combined_enrichR_Down_0.05)
      
      if (length(Foley_Runs_Combined_enrichR_Down_0.05$adj.P.Val) > 10) {
        Foley_Runs_Combined_enriched_Down <- enrichr(as.character(Foley_Runs_Combined_enrichR_Down_0.05$SYMBOL), enrichR_dbs)
      }
      
      Foley_Runs_Combined_enrichR_Up <- Foley_Runs_Combined_enrichR[which(Foley_Runs_Combined_enrichR$logFC > 0), ]
      Foley_Runs_Combined_enrichR_Up_0.05 <- Foley_Runs_Combined_enrichR_Up[which(Foley_Runs_Combined_enrichR_Up$adj.P.Val < 0.05), ]
      dim(Foley_Runs_Combined_enrichR_Up_0.05)
      
      if (length(Foley_Runs_Combined_enrichR_Up_0.05$adj.P.Val) > 10) {
        Foley_Runs_Combined_enriched_Up <- enrichr(as.character(Foley_Runs_Combined_enrichR_Up_0.05$SYMBOL), enrichR_dbs)
      }
      
      for (k in 1:length(enrichR_dbs)) {
        
        if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_top800"))) { dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_top800"), recursive = TRUE) }
        
        write.csv(Foley_Runs_Combined_enrichR_top800[[k]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_top800"), "/EnrichR_top800_Foley_Runs_Combined_", 
                                                                paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                                if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                                "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        
        if (length(Foley_Runs_Combined_enrichR_Down_0.05$adj.P.Val) > 10) {
          if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"))) {dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"), recursive = TRUE) }
          
          write.csv(Foley_Runs_Combined_enriched_Down[[j]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_Down"), "/EnrichR_Down_Foley_Runs_Combined_", 
                                                                 paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                                 if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                                 "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        }
        
        if (length(Foley_Runs_Combined_enrichR_Up_0.05$adj.P.Val) > 10) {
          if(!dir.exists(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"))) {dir.create(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"), recursive = TRUE)}
          
          write.csv(Foley_Runs_Combined_enriched_Up[[j]], file = paste0(paste0(Directory_CSVs_DEGenes, "/EnrichR_Up"), "/EnrichR_Up_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                                               if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])},"_", enrichR_dbs[k],
                                                               "_adjP", gsub(".", "", as.character(adjP), fixed = TRUE), ".csv"))
        }
      }
      
      
      
      
      
    }
  }
} 




#par(mfrow = c(1,3))
#BvsC_Foley_Runs_Combined <- topTreat(tfit_Foley_Runs_Combined, coef = 3, n = Inf)
#head(BvsC_Foley_Runs_Combined)
#class(BvsC_Foley_Runs_Combined)
#dim(BvsC_Foley_Runs_Combined)
#BvsC_Foley_Runs_Combined_NA_removed <- BvsC_Foley_Runs_Combined[!is.na(BvsC_Foley_Runs_Combined$SYMBOL),]
#head(BvsC_Foley_Runs_Combined_NA_removed$SYMBOL)
#dim(BvsC_Foley_Runs_Combined_NA_removed)
#par(mfrow = c(1,1))
#plotMD(tfit_Foley_Runs_Combined, column=3, status=dt_Foley_Runs_Combined[,3], xlim=c(-8,13), main = "BvsC mean-difference plot")
#glMDPlot(tfit_Foley_Runs_Combined, coef=3, status=dt_Foley_Runs_Combined, main=colnames(tfit_Foley_Runs_Combined)[3], side.main="ENSEMBL", 
#counts=lcpm_Foley_Runs_Combined, groups=combined_data_Foley_Runs_Combined$samples$group,launch=FALSE)
#main=colnames(tfit_Foley_Runs_Combined)[3]

#AvsB_Foley_Runs_Combined <- topTreat(tfit_Foley_Runs_Combined, coef = 1, n = Inf)
#head(AvsB_Foley_Runs_Combined)
#dim(AvsB_Foley_Runs_Combined)
#AvsB_Foley_Runs_Combined_NA_removed <- AvsB_Foley_Runs_Combined[!is.na(AvsB_Foley_Runs_Combined$SYMBOL),]
#head(AvsB_Foley_Runs_Combined_NA_removed$SYMBOL)
#plotMD(tfit_Foley_Runs_Combined, column = 1, status = dt_Foley_Runs_Combined[,1], xlim = c(-8,13), main = "AvsB mean-difference plot")

#AvsC_Foley_Runs_Combined <- topTreat(tfit_Foley_Runs_Combined, coef = 2, n = Inf)
#head(AvsC_Foley_Runs_Combined)
#AvsC_Foley_Runs_Combined_NA_removed <- AvsC_Foley_Runs_Combined[!is.na(AvsC_Foley_Runs_Combined$SYMBOL),]
#head(AvsC_Foley_Runs_Combined_NA_removed$SYMBOL)
#plotMD(tfit_Foley_Runs_Combined, column = 2, status = dt_Foley_Runs_Combined[,2], xlim = c(-8,13), main = "AvsC mean-difference plot")

####HEATMAPS####

par(mfrow = c(1,1))
for (i in 1:length(Comparisons_vector)) {
  if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/HeatMap_lcpm_Foley_Runs_Combined_", paste(trimws(basename(Directory_Plots_RNASeq)), collapse=""),
                                           if (TwoGroups == FALSE) {paste0("_", Comparisons_vector[i])}, ".pdf", sep = ""), width=16.25, height=12.5)}
  
  j <- which(v_Foley_Runs_Combined$genes$ENSEMBL %in% Ranked_DEGenes_list[[i]]$ENSEMBL[1:100])
  mycol <- colorpanel(1000, "blue", "white", "red")
  heatmap.2(lcpm_Foley_Runs_Combined[j,], scale = "row", labRow = v_Foley_Runs_Combined$genes$SYMBOL[j], 
            labCol = sample_name_Foley_Runs_Combined, col = mycol, trace = "none", Colv = FALSE, #Colv reorders the columns of heatmap
            density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'row', distfun = function(x) as.dist(1-cor(t(x), method = "s")),
            main = paste0(Comparisons_vector[i], " Heat Map"))
  if (SavePlots == TRUE) {dev.off()}
}          


#BvsC_topgenes <- BvsC_Foley_Runs_Combined_NA_removed$ENSEMBL[1:80]
#i <- which(v_Foley_Runs_Combined$genes$ENSEMBL %in% BvsC_topgenes)
#mycol <- colorpanel(1000, "blue", "white", "red")
#heatmap.2(lcpm_Foley_Runs_Combined[i,], scale = "row", labRow = v_Foley_Runs_Combined$genes$SYMBOL[i], labCol = sample_name_Foley_Runs_Combined, col = mycol, trace = "none",
#         density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'both', distfun = function(x) as.dist(1-cor(t(x), method = "s")), main = "BvsC HeatMap")


#}

##Bar plot of DEGenes (P Value)
# experiment_name <- "Foley_Runs_Combined"
# bar_plot_names <- gsub(filter_type, "", list.files("/Users/gordid/Desktop/Foley_Runs_Combined/Ranked_DEGenes/")) %>% gsub(pattern = paste0("Ranked-DEGenes_", experiment_name), replacement = "") %>%
#   gsub(pattern = "_[^_]+$", replacement = "") %>% gsub(pattern = "^_", replacement = "")
# bar_plot_names
# 
# Ranked_DEGenes_data_Foley_Runs_Combined <- list()
# Total_DEGenes_Foley_Runs_Combined <- data.frame()
# setwd("/Users/gordid/Desktop/Foley_Runs_Combined/Ranked_DEGenes/")
# for (i in 1:length(list.files())) {
#   Ranked_DEGenes_data_Foley_Runs_Combined[[i]] <- read.csv(list.files()[i])
#   print(list.files()[i])
#   Total_DEGenes_Foley_Runs_Combined[i,1] <- bar_plot_names[i]
#   Total_DEGenes_Foley_Runs_Combined[i,2] <- length(which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$P.Value[which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$logFC > 0)] < 0.05))
#   Total_DEGenes_Foley_Runs_Combined[i,3] <- length(which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$P.Value[which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$logFC < 0)] < 0.05))
# }
# names(Ranked_DEGenes_data_Foley_Runs_Combined) <- bar_plot_names
# colnames(Total_DEGenes_Foley_Runs_Combined) <- c("Comparisons", "DEG_Up_Counts", "DEG_Down_Counts")
# Total_DEGenes_Foley_Runs_Combined <- Total_DEGenes_Foley_Runs_Combined[order(rowSums(Total_DEGenes_Foley_Runs_Combined[,c(2,3)])),]
# Total_DEGenes_Foley_Runs_Combined$Comparisons <- factor(as.factor(Total_DEGenes_Foley_Runs_Combined$Comparisons), levels = unique(as.factor(Total_DEGenes_Foley_Runs_Combined$Comparisons)) )
# Total_DEGenes_Foley_Runs_Combined
# 
# Total_DEGenes_Foley_Runs_Combined <- melt(Total_DEGenes_Foley_Runs_Combined, id.var="Comparisons")
# Total_DEGenes_Foley_Runs_Combined
# setwd("/Users/gordid/Desktop/Foley_Runs_Combined/")
# 
# if (SavePlots == TRUE) {pdf(file = paste0("/Users/gordid/Desktop/Foley_Runs_Combined/All_results", "/Total_DEGenes_P05_perComparison_BarGraph_Foley_Runs_Combined",
#                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".pdf"), width=30, height=12.5)}
# 
# ggplot(Total_DEGenes_Foley_Runs_Combined, aes(x = Comparisons, y = value, fill = variable)) +
#   geom_bar(stat = "identity") + xlab("Sample Comparisons") + ylab("DEG Counts (P < 0.05)") + coord_flip() + scale_fill_manual(values=c("red", "blue"))
# 
# if (SavePlots == TRUE) {dev.off()}
# 
# ##Bar plot of DEGenes (AdjP Value)
# Ranked_DEGenes_data_Foley_Runs_Combined <- list()
# Total_DEGenes_Foley_Runs_Combined <- data.frame()
# setwd("/Users/gordid/Desktop/Foley_Runs_Combined/Ranked_DEGenes/")
# for (i in 1:length(list.files())) {
#   Ranked_DEGenes_data_Foley_Runs_Combined[[i]] <- read.csv(list.files()[i])
#   print(list.files()[i])
#   Total_DEGenes_Foley_Runs_Combined[i,1] <- bar_plot_names[i]
#   Total_DEGenes_Foley_Runs_Combined[i,2] <- length(which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$adj.P.Val[which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$logFC > 0)] < 0.05))
#   Total_DEGenes_Foley_Runs_Combined[i,3] <- length(which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$adj.P.Val[which(Ranked_DEGenes_data_Foley_Runs_Combined[[i]]$logFC < 0)] < 0.05))
# }
# names(Ranked_DEGenes_data_Foley_Runs_Combined) <- bar_plot_names
# colnames(Total_DEGenes_Foley_Runs_Combined) <- c("Comparisons", "DEG_Up_Counts", "DEG_Down_Counts")
# Total_DEGenes_Foley_Runs_Combined <- Total_DEGenes_Foley_Runs_Combined[order(rowSums(Total_DEGenes_Foley_Runs_Combined[,c(2,3)])),]
# Total_DEGenes_Foley_Runs_Combined$Comparisons <- factor(as.factor(Total_DEGenes_Foley_Runs_Combined$Comparisons), levels = unique(as.factor(Total_DEGenes_Foley_Runs_Combined$Comparisons)) )
# Total_DEGenes_Foley_Runs_Combined
# 
# Total_DEGenes_Foley_Runs_Combined <- melt(Total_DEGenes_Foley_Runs_Combined, id.var="Comparisons")
# Total_DEGenes_Foley_Runs_Combined
# setwd("/Users/gordid/Desktop/Foley_Runs_Combined/")
# 
# if (SavePlots == TRUE) {pdf(file = paste0("/Users/gordid/Desktop/Foley_Runs_Combined/All_results", "/Total_DEGenes_AdjP05_perComparison_BarGraph_Foley_Runs_Combined",
#                                           "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), ".pdf"), width=30, height=12.5)}
# 
# ggplot(Total_DEGenes_Foley_Runs_Combined, aes(x = Comparisons, y = value, fill = variable)) +
#   geom_bar(stat = "identity") + xlab("Sample Comparisons") + ylab("DEG Counts (AdjP < 0.05)") + coord_flip() + scale_fill_manual(values=c("red", "blue"))
# 
# if (SavePlots == TRUE) {dev.off()}
