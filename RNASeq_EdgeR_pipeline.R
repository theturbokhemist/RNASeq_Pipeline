#read in data file on experiment
#setwd("/Users/gordid/Desktop/RNA_seq/RNA-seq analysis practice/Lucas_Lupus_Data")
setwd("/Users/gordid/Desktop/Lucas")
experiment_lucas = read.csv("Chang_RNAseq sample submission.csv", header = TRUE)
experiment_lucas

#Bypassing missing data description: group names and number of files#

#creating vector of data files
files_lucas <- c("1_GT18-15666_GGTCACGA-CATAATAC_S12.genes.results.withGeneName.txt", "2_GT18-15667_CTGCTTCC-GATCTATC_S8.genes.results.withGeneName.txt",
                 "3_GT18-15668_TCATCCTT-AGCTCGCT_S9.genes.results.withGeneName.txt", "4_GT18-15669_AGGTTATA-CGGAACTG_S5.genes.results.withGeneName.txt",
                 "5_GT18-15670_GAACCGCG-TAAGGTCA_S3.genes.results.withGeneName.txt", "6_GT18-15671_CTCACCAA-TTGCCTAG_S6.genes.results.withGeneName.txt",
                 "7_GT18-15672_TCTGTTGG-CCATTCGA_S1.genes.results.withGeneName.txt", "8_GT18-15673_ACACTAAG-ATATGGAT_S2.genes.results.withGeneName.txt",
                 "9_GT18-15674_GTGTCGGA-GCGCAAGC_S10.genes.results.withGeneName.txt", "10_GT18-15675_TTCCTGTT-AAGATACT_S4.genes.results.withGeneName.txt",
                 "11_GT18-15676_CCTTCACC-GGAGCGTC_S11.genes.results.withGeneName.txt", "12_GT18-15677_GCCACAGG-ATGGCATG_S7.genes.results.withGeneName.txt")

####What Organism?####
Mouse = TRUE 

####Are there Outliers?####
outliers = FALSE

####Are you using CPM or counts to filter?####
CPMFilter = TRUE
if (CPMFilter == TRUE) {
  filter_type = "_cpm-filter"
} else {
  filter_type = "_counts-filter"
}

####Comparing just 2 groups####

TwoGroups = TRUE
WhichGroups = c(1:4,5:8)
#Better if you could put in comparisons you are interested in

if (TwoGroups == TRUE) {
  files_lucas = files_lucas[WhichGroups]
  files_lucas
}

####Save plots?####
SavePlots = FALSE

#reading in one of the data files
file_1 = read.delim(files_lucas[1], header=TRUE)
head(file_1)
str(file_1)
dim(file_1)

####Renaming columns and creating groups####

sample_name_lucas_simple <- gsub("_GT18.*", "", files_lucas)
sample_name_lucas_simple

group_name_lucas <- as.character(experiment_lucas$Group.name)
group_name_lucas <- as.factor(group_name_lucas)
group_name_lucas

if (TwoGroups == TRUE) {
  group_name_lucas <- as.character(experiment_lucas$Group.name)[WhichGroups]
  group_name_lucas <- as.factor(group_name_lucas)
  group_name_lucas
}
levels(group_name_lucas)

sample_names_lucas_final1 <- paste(sample_name_lucas_simple, group_name_lucas)
sample_names_lucas_final1
sample_names_lucas_final <- gsub(" ", "", sample_names_lucas_final1)
sample_names_lucas_final
class(sample_names_lucas_final)

####creating vector of gene names and vector of possible sample comparisons####
genes_lucas_file <- file_1$GeneName
head(genes_lucas_file)
length(genes_lucas_file)

Comparisons_vector <- c()
for (i in 1:(length(levels(group_name_lucas))-1)) {
  for (j in (i+1):length(levels(group_name_lucas))) {

    Comparisons_vector <- append(Comparisons_vector, paste0(levels(group_name_lucas)[i], "vs", levels(group_name_lucas)[j]))
    
  }
}

Comparisons_vector
length(Comparisons_vector)

#Creating new directorys####

#plots/RNASeq folder
Directory_Plots_RNASeq <- paste("/Users/gordid/Desktop/Lucas/All_results/Plots/RNASeq/Regular", filter_type,
                              (if(TwoGroups == TRUE) { 
                                paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                              } else {
                                "_all-data"
                              }), sep = "")

if (dir.exists(Directory_Plots_RNASeq) == FALSE) {
dir.create(paste("/Users/gordid/Desktop/Lucas/All_results/Plots/RNASeq/Regular", filter_type,
                 (if(TwoGroups == TRUE) { 
                   paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                   } else {
                     "_all-data"
                   }), sep = ""))
} else {
  "Already made this directory!"
}

#plots/NetAct_RACIPE folder

Directory_Plots_NetAct_RACIPE <- paste("/Users/gordid/Desktop/Lucas/All_results/Plots/NetAct_RACIPE/Regular", filter_type,
                                (if(TwoGroups == TRUE) { 
                                  paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                                } else {
                                  "_all-data"
                                }), sep = "")

if (dir.exists(Directory_Plots_NetAct_RACIPE) == FALSE) {
  dir.create(paste("/Users/gordid/Desktop/Lucas/All_results/Plots/NetAct_RACIPE/Regular", filter_type,
                   (if(TwoGroups == TRUE) { 
                     paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                   } else {
                     "_all-data"
                   }), sep = ""))
} else {
  "Already made this directory!"
}


#CSVs/RNASeq folder
Directory_CSVs_DEGenes <- paste("/Users/gordid/Desktop/Lucas/All_results/CSVs/DEGenes/Regular", filter_type,
                                (if(TwoGroups == TRUE) { 
                                  paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                                } else {
                                  "_all-data"
                                }), sep = "")

if (dir.exists(Directory_CSVs_DEGenes) == FALSE) {
  dir.create(paste("/Users/gordid/Desktop/Lucas/All_results/CSVs/DEGenes/Regular", filter_type,
                   (if(TwoGroups == TRUE) { 
                     paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                   } else {
                     "_all-data"
                   }), sep = ""))
} else {
  "Already made this directory!"
}

#CSVs/NetAct_RACIPE folder
Directory_CSVs_NetAct_RACIPE <- paste("/Users/gordid/Desktop/Lucas/All_results/CSVs/NetAct_RACIPE/Regular", filter_type,
                                (if(TwoGroups == TRUE) { 
                                  paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                                } else {
                                  "_all-data"
                                }), sep = "")

if (dir.exists(Directory_CSVs_NetAct_RACIPE) == FALSE) {
  dir.create(paste("/Users/gordid/Desktop/Lucas/All_results/CSVs/NetAct_RACIPE/Regular", filter_type,
                   (if(TwoGroups == TRUE) { 
                     paste("_", levels(group_name_lucas)[1], "vs", levels(group_name_lucas)[2], sep = "")
                   } else {
                     "_all-data"
                   }), sep = ""))
} else {
  "Already made this directory!"
}

####CREATE DGEList OBJECT####

combined_data_lucas <- readDGE(files_lucas, columns = c(1, 5))
dim(combined_data_lucas) #showing dimensions of "counts" matrix
class(combined_data_lucas) #showing class of the whole object
summary(combined_data_lucas) #showing names, size, class, and mode of each element of the object, not sure why class is "none" for counts
str(combined_data_lucas) #showing the structure of the object
head(combined_data_lucas$counts)

####rename the  sample ID's in our "counts" element by subtracting the last 27 characters of the sample names####
samplenames_lucas <- substring(colnames(combined_data_lucas$counts), 1, nchar(colnames(combined_data_lucas$counts))-27)
samplenames_lucas
colnames(combined_data_lucas) <- samplenames_lucas

####Assigning experimental info to categories in the DGEList obect####
if (TwoGroups == TRUE) {
  group_name_lucas <- as.factor(as.character(experiment_lucas$Group.name)[WhichGroups])
  age_lucas <- as.factor(as.character(experiment_lucas$Age)[WhichGroups])
  x2dg_lucas <- as.factor(as.character(experiment_lucas$X2DG)[WhichGroups])
  remarks_lucas <- as.factor(as.character(experiment_lucas$Remarks)[WhichGroups])
  } else {
    
    age_lucas <- as.factor(experiment_lucas$Age)
    x2dg_lucas <- as.factor(experiment_lucas$X2DG)
    remarks_lucas <- as.factor(experiment_lucas$Remarks)
  }

combined_data_lucas$samples$group <- group_name_lucas
combined_data_lucas$samples

combined_data_lucas$samples$age <- age_lucas

combined_data_lucas$samples$x2dg <- x2dg_lucas

combined_data_lucas$samples$remarks <- remarks_lucas

####Create new genes matrix by retrieving gene information from the mus library. The gene information of their SYMBOL and TXCHROM#### 
columns(Mus.musculus)
gene_id_lucas <- rownames(combined_data_lucas)
head(gene_id_lucas)
length(gene_id_lucas)
if (Mouse == TRUE) {
  genes_lucas <- select(Mus.musculus, keys = gene_id_lucas, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
  } else {
  genes_lucas <- select(Homo.sapiens, keys = gene_id_lucas, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
    
  }
dim(genes_lucas)
head(genes_lucas)

####symbols are duplicated too often and so I dont remove genes based on duplicated symbols (cause of NA mostly, will remove NAs later)####
head(which(duplicated(genes_lucas$SYMBOL)))
length(which(duplicated(genes_lucas$SYMBOL)))
genes_lucas[c(243, 246, 546, 741, 1111, 1220, 1364, 1449, 1473, 1885, 2149, 2187, 2311, 2607, 2648, 2807 ), ]

####remove duplicated genes based on duplicated ENSEMBL ID####
dup_genes_lucas <- which(duplicated(genes_lucas$ENSEMBL))
head(dup_genes_lucas)
length(dup_genes_lucas)
head(genes_lucas[dup_genes_lucas, ]) 
genes_lucas <- genes_lucas[-dup_genes_lucas, ]
dim(genes_lucas)
dim(combined_data_lucas$counts)

combined_data_lucas$genes <- genes_lucas #assign data frame of gene annotations to the combined_data object
summary(combined_data_lucas)

####CPM normalization on count data####
cpm_lucas <- cpm(combined_data_lucas)
lcpm_lucas <-cpm(combined_data_lucas, log = TRUE)

L_lucas <- mean(combined_data_lucas$samples$lib.size) * 1e-6 #calculation of L for the prior count.
M_lucas <- median(combined_data_lucas$samples$lib.size) * 1e-6
c(L_lucas, M_lucas)
head(lcpm_lucas)

table(rowSums(combined_data_lucas$counts==0)==12)

###FILTER THE GENES####
if (CPMFilter == TRUE) {
  kept_expression_lucas <- filterByExpr(combined_data_lucas, group = combined_data_lucas$samples$group) #filterByExpr function
  # keeps genes with 10 read counts or more in a minimum number of samples, where the minimum number of samples is chosen according 
  #to the minimum group sample size. 
  }  else {
    #filter by counts, rather than by cpm
    kept_expression_lucas <- (rowSums(combined_data_lucas$counts > 10) >= 1)
}

length(kept_expression_lucas)
head(kept_expression_lucas)

#removing all the lowly expressed genes from the combined_data object. When you subset a DGEList and specify keep.lib.sizes=FALSE,
#the lib.size for each sample will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
filtered_combined_data_lucas <- combined_data_lucas[kept_expression_lucas, , keep.lib.sizes = FALSE] 
#filter_cpm <- dim(filtered_combined_data_lucas$counts)
filter_counts <- dim(filtered_combined_data_lucas$counts)



####Removing duplicated gene symbols for the NetAct analysis####
Combined_data_lucas_NetAct <- filtered_combined_data_lucas
filter_counts
dim(Combined_data_lucas_NetAct$counts)
length(Combined_data_lucas_NetAct$genes$SYMBOL)
rownames(Combined_data_lucas_NetAct$counts) <- Combined_data_lucas_NetAct$genes$SYMBOL

length(rownames(Combined_data_lucas_NetAct$counts))
head(Combined_data_lucas_NetAct$genes$SYMBOL)
head(Combined_data_lucas_NetAct$counts)

length(na.omit(rownames(Combined_data_lucas_NetAct$counts)))

dup_genes_lucas_2 <- which(duplicated(Combined_data_lucas_NetAct$genes$SYMBOL))
length(dup_genes_lucas_2)
Combined_data_lucas_NetAct <- Combined_data_lucas_NetAct[-dup_genes_lucas_2, ]
dim(Combined_data_lucas_NetAct$counts)

which(is.na(rownames(Combined_data_lucas_NetAct$counts)))

Combined_data_lucas_NetAct <- Combined_data_lucas_NetAct[-which(is.na(rownames(Combined_data_lucas_NetAct$counts))), ]
dim(Combined_data_lucas_NetAct$counts)
filter_cpm_2 <- dim(Combined_data_lucas_NetAct$counts)
dim(combined_data_lucas$counts)
#na.omit(rownames(Combined_data_lucas_NetAct$counts))


####Plotting filtered Data####
lcpm.cutoff_lucas <- log2(10/M_lucas + 2/L_lucas)
lcpm.cutoff_lucas
nsamples <- ncol(combined_data_lucas$counts)
col <- brewer.pal(nsamples, "Paired")

if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/FilteredData_Lucas", filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), ".pdf", sep = ""))}
par(mfrow = c(1,2))
plot(density(lcpm_lucas[, 1]), col = col[1], lwd=2, ylim=c(0,0.46), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_lucas, lty=3)
title(main="A. Raw data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_lucas[, i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

combined_data_lucas <- filtered_combined_data_lucas
lcpm_lucas <- cpm(combined_data_lucas, log = TRUE)

plot(density(lcpm_lucas[, 1]), col = col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="" )
abline(v=lcpm.cutoff_lucas, lty=3)
title(main="B. Filtered data", xlab="Log-cpm")
for (i in 2:nsamples) {
  den = density(lcpm_lucas[, i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
if (SavePlots == TRUE) {dev.off() }


####Calculate normalization factors and creates new DGEList object with these new normalization factors inserted. THEN GENERATING BOX PLOT####
colnames(lcpm_lucas) <- sample_names_lucas_final
combined_data_lucas$samples$norm.factors

if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/TMM_Normalized_Lucas", filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), ".pdf", sep = ""))}
par(mfrow = c(1,2))
boxplot(lcpm_lucas, las = 2, col = col, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")

combined_data_lucas = calcNormFactors(combined_data_lucas, method = "TMM") #calculates normalization factors and creates new DGEList object
#these normalization factors are important because if there is a sample where expression level is different accross the board, that may have been caused
#by an eternal factor that is not of biological interest. TMM helps normalize for that, using the assumption that should have a similar range
#and distribution of expression values.
lcpm_lucas <- cpm(combined_data_lucas, log = TRUE) 
colnames(lcpm_lucas) <- sample_names_lucas_final
head(lcpm_lucas)

combined_data_lucas$samples$norm.factors
boxplot(lcpm_lucas, las = 2, col = col, main="")
title(main="B.Normalised with TMM data", ylab="Log-cpm")
if (SavePlots == TRUE) {dev.off() }


####PCA and PCoA Plots####
if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/PCoA-PCA_1vs2_Lucas", filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), ".pdf", sep = ""))}
color_groups_lucas <- combined_data_lucas$samples$group
levels(color_groups_lucas) <- brewer.pal(nlevels(color_groups_lucas), "Set1")
color_groups_lucas <- as.character(color_groups_lucas)
color_groups_lucas

combined_data_lucas$samples$names <- sample_names_lucas_final
combined_data_lucas$samples

par(mfrow = c(1,2))
plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$names, col = color_groups_lucas, main = "PCoA")
plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$names, col = color_groups_lucas, gene.selection="common", main = "PCA")
if (SavePlots == TRUE) {dev.off() }

par(mfrow = c(1,2))
plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$names, col = color_groups_lucas, dim = c(2,3), main = "PCoA")
plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$names, col = color_groups_lucas, gene.selection="common", dim = c(2,3), main = "PCA")

par(mfrow = c(1,1))
plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$group, col = color_groups_lucas, gene.selection="common", main = "PCA")

glMDSPlot(lcpm_lucas, labels=paste(combined_data_lucas$samples$names, combined_data_lucas$samples$age, sep="_"), groups = combined_data_lucas$samples[,c(2,5)], main = "PCoA", launch=FALSE)
glMDSPlot(lcpm_lucas, labels=paste(combined_data_lucas$samples$names, combined_data_lucas$samples$age, sep="_"), groups = combined_data_lucas$samples[,c(2,5)],
          gene.selection="common", main = "PCA", launch=TRUE)

#scree.plot(plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$names, col = color_groups_lucas, gene.selection="common", main = "PCA"))
#a <- (plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$names, col = color_groups_lucas, main = "PCoA"))
#screeplot(a, xlim = c(-8,13), ylim=c(0,0.46))
####PCA Plots for other phenotypic data and meta-data####
#color_age_lucas <- combined_data_lucas$samples$age
#levels(color_age_lucas) <- brewer.pal(nlevels(color_age_lucas), "Set2")
#color_age_lucas <- as.character(color_age_lucas)
#color_age_lucas
#plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$age, col = color_groups_lucas)
#plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$age, col = color_groups_lucas, gene.selection="common")

#color_x2dg_lucas <- combined_data_lucas$samples$x2dg
#levels(color_x2dg_lucas) <- brewer.pal(nlevels(color_x2dg_lucas), "Set3")
#color_x2dg_lucas <- as.character(color_x2dg_lucas)
#color_x2dg_lucas
#plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$x2dg, col = color_groups_lucas)
#plotMDS(lcpm_lucas, labels = combined_data_lucas$samples$x2dg, col = color_groups_lucas, gene.selection="common")

####Building Design Matrix####

design_matrix_lucas_no.int <- model.matrix(~0+combined_data_lucas$samples$group) #no intercept
colnames(design_matrix_lucas_no.int) = gsub("combined_data_lucas$samples$group", "", colnames(design_matrix_lucas_no.int), fixed = TRUE)
design_matrix_lucas_no.int

design_matrix_lucas_int <- model.matrix(~combined_data_lucas$samples$group) #with intercept
colnames(design_matrix_lucas_int) = gsub("combined_data_lucas$samples$group", "", colnames(design_matrix_lucas_int), fixed = TRUE)
design_matrix_lucas_int


if (TwoGroups == TRUE) {
 contr.matrix_no.int <- makeContrasts(
 paste(levels(group_name_lucas)[1]," - ", levels(group_name_lucas)[2], sep = ""),
 levels = colnames(design_matrix_lucas_no.int))
 colnames(contr.matrix_no.int) <- paste(levels(group_name_lucas)[1],"vs", levels(group_name_lucas)[2], sep = "")
} else {
  contrast_vector <- c()
  column_name_vector <- c()
  for (i in 1:(length(levels(group_name_lucas))-1)) {
    
    for (j in (i+1):length(levels(group_name_lucas))) {
      
      contrast_vector <- append(contrast_vector, paste(levels(group_name_lucas)[i]," - ", levels(group_name_lucas)[j], sep = ""))
      column_name_vector <- append(column_name_vector, paste(levels(group_name_lucas)[i], "vs", levels(group_name_lucas)[j], sep = ""))
    }
    
  }
  
  contr.matrix_no.int <- makeContrasts(
    contrasts = contrast_vector,
    levels = colnames(design_matrix_lucas_no.int)
  )
colnames(contr.matrix_no.int) <- column_name_vector
contr.matrix_no.int
}
  

#contr.matrix_no.int <- makeContrasts(
  #AvsB = "A- B",
  #AvsC = A - C,
  #BvsC = B - C,
  #levels = colnames(design_matrix_lucas_no.int))
#contr.matrix_no.int
#colnames(contr.matrix_no.int) <- c("a", "b", "c")

#contr.matrix_no.int <- makeContrasts(
 #contrasts = contrast_vector,
  #levels = colnames(design_matrix_lucas_no.int))
#contr.matrix_no.int

#contr.matrix_int <- makeContrasts( #Intercept column has replcaced "A", not sure how contrast matrix works now
# AvsB = A - B,
# AvsC = A - C,
# BvsC = B - C,
#levels = colnames(design_matrix_lucas_int))
#contr.matrix_int

####VOOM####
if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Voom_Lucas", filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,2))
v_lucas <- voom(combined_data_lucas, design_matrix_lucas_no.int, plot=TRUE)
str(v_lucas)
dim(v_lucas$genes)

vfit_lucas <- lmFit(v_lucas, design_matrix_lucas_no.int)
summary(decideTests(vfit_lucas))
str(vfit_lucas)

vfit_lucas_contrasts <- contrasts.fit(vfit_lucas, contrasts=contr.matrix_no.int) 
summary(decideTests(vfit_lucas_contrasts))

efit_lucas <- eBayes(vfit_lucas_contrasts)
summary(decideTests(efit_lucas))

plotSA(efit_lucas)
title(main="voom: Mean-variance trend-Normalized")
if (SavePlots == TRUE) {dev.off()}

summary(decideTests(efit_lucas))

####Filter further with logfold change and Venn Diagram####
#Figure out the criteria data has to choose appropriate lfc
tfit_lucas <- treat(vfit_lucas_contrasts, lfc=0.1)
dt_lucas <- decideTests(tfit_lucas)
summary(dt_lucas)
lfc = 0.1

#which genes are DE in multiple comparisons    
#Figure out how to add this step to the pipeline
if (TwoGroups == FALSE) {
de.common_lucas_2_comparisons <- which(dt_lucas[,1]!=0 & dt_lucas[,2]!=0)
length(de.common_lucas_2_comparisons)

de.common_lucas_3_comparisons <- which(dt_lucas[,1]!=0 & dt_lucas[,2]!=0 & dt_lucas[,3]!=0)
length(de.common_lucas_3_comparisons)

if (SavePlots == TRUE) {pdf(paste(Directory_Plots_RNASeq, "/Venn_Lucas_lfc",gsub(".", "", as.character(lfc), fixed = TRUE), filter_type, gsub(".*filter", "", Directory_Plots_RNASeq),
                                  ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,1))
vennDiagram(dt_lucas[,1:3], circle.col=c("turquoise", "salmon", "orange"))
if (SavePlots == TRUE) {dev.off()}
}

write.fit(tfit_lucas, dt_lucas, file = paste(Directory_CSVs_DEGenes, "/tfit_Lucas_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), filter_type, gsub(".*filter", "", Directory_Plots_RNASeq),
          if (outliers == TRUE) {print("_outliers-removed")}, ".txt", sep = ""))

tfit_lucas.txt = read.table(file = paste(Directory_CSVs_DEGenes, "/tfit_Lucas_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), filter_type, gsub(".*filter", "", Directory_Plots_RNASeq),
                                                if (outliers == TRUE) {print("_outliers-removed")}, ".txt", sep = ""))

write.csv(tfit_lucas.txt, file = paste(Directory_CSVs_DEGenes, "/tfit_Lucas_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), filter_type, gsub(".*filter", "", Directory_Plots_RNASeq),
                                       if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))


####Mean-difference plots and arranging DE genes####
if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/Mean-Difference-Plots_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                         filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), if (outliers == TRUE) {print("_outliers-removed")}, ".pdf", sep = ""), width=13, height=10)}
par(mfrow = c(1,length(Comparisons_vector)))
for (i in 1:length(Comparisons_vector)) {
  plotMD(tfit_lucas, column = i, status = dt_lucas[,i], xlim = c(-8,13), main = paste0(Comparisons_vector[i], " Mean-Difference Plot"))
  glMDPlot(tfit_lucas, coef=i, status=dt_lucas, main=colnames(tfit_lucas)[i], side.main="ENSEMBL", counts=lcpm_lucas, groups=combined_data_lucas$samples$group,launch=FALSE)
  
}
if (SavePlots == TRUE){ dev.off()}

#par(mfrow = c(1,length(levels(group_name_lucas))))
#Ranked_DEGenes_list <- list()
#for (i in 1:(length(levels(group_name_lucas))-1)) {
#  for (j in (i+1):length(levels(group_name_lucas))) {
#    print(paste0(levels(group_name_lucas)[i], "vs", levels(group_name_lucas)[j]))
#    assign(paste0(levels(group_name_lucas)[i], "vs", levels(group_name_lucas)[j]), (topTreat(tfit_lucas, coef = i, n = Inf)))
#    head(eval(as.name(paste0(levels(group_name_lucas)[i], "vs", levels(group_name_lucas)[j]))))
#  }
#  
#} 

par(mfrow = c(1,length(levels(group_name_lucas))))
Ranked_DEGenes_list <- list()
k = 0
for (i in 1:(length(levels(group_name_lucas))-1)) {
  for (j in (i+1):length(levels(group_name_lucas))) {
    k = k+1 
    Ranked_DEGenes_list[[k]] <- topTreat(tfit_lucas, coef = k, n = Inf)
  }
} 
names(Ranked_DEGenes_list) <- Comparisons_vector
summary(Ranked_DEGenes_list)
dim(Ranked_DEGenes_list[[1]])

for (i in 1:length(Comparisons_vector)) {
  Ranked_DEGenes_list[[i]] <- Ranked_DEGenes_list[[i]][!is.na(Ranked_DEGenes_list[[i]]$SYMBOL),]
}
str(Ranked_DEGenes_list)
dim(Ranked_DEGenes_list[[1]])

#removing duplicated Symbols#
for (i in 1:length(Comparisons_vector)) {
  Ranked_DEGenes_list[[i]] <- Ranked_DEGenes_list[[i]][-which(duplicated(Ranked_DEGenes_list[[i]]$SYMBOL)),]
}
str(Ranked_DEGenes_list)
dim(Ranked_DEGenes_list[[1]])

for (i in 1:length(Comparisons_vector)) {
  if (TwoGroups == FALSE) {
  write.csv(Ranked_DEGenes_list[[i]], file = paste(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Lucas", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                    filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
} else {
  write.csv(Ranked_DEGenes_list[[i]], file = paste(Directory_CSVs_DEGenes, "/Ranked-DEGenes_Lucas", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                                   filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
}
}
  
#par(mfrow = c(1,3))
#BvsC_lucas <- topTreat(tfit_lucas, coef = 3, n = Inf)
#head(BvsC_lucas)
#class(BvsC_lucas)
#dim(BvsC_lucas)
#BvsC_lucas_NA_removed <- BvsC_lucas[!is.na(BvsC_lucas$SYMBOL),]
#head(BvsC_lucas_NA_removed$SYMBOL)
#dim(BvsC_lucas_NA_removed)
#par(mfrow = c(1,1))
#plotMD(tfit_lucas, column=3, status=dt_lucas[,3], xlim=c(-8,13), main = "BvsC mean-difference plot")
#glMDPlot(tfit_lucas, coef=3, status=dt_lucas, main=colnames(tfit_lucas)[3], side.main="ENSEMBL", counts=lcpm_lucas, groups=combined_data_lucas$samples$group,launch=FALSE)
#main=colnames(tfit_lucas)[3]

#AvsB_lucas <- topTreat(tfit_lucas, coef = 1, n = Inf)
#head(AvsB_lucas)
#dim(AvsB_lucas)
#AvsB_lucas_NA_removed <- AvsB_lucas[!is.na(AvsB_lucas$SYMBOL),]
#head(AvsB_lucas_NA_removed$SYMBOL)
#plotMD(tfit_lucas, column = 1, status = dt_lucas[,1], xlim = c(-8,13), main = "AvsB mean-difference plot")

#AvsC_lucas <- topTreat(tfit_lucas, coef = 2, n = Inf)
#head(AvsC_lucas)
#AvsC_lucas_NA_removed <- AvsC_lucas[!is.na(AvsC_lucas$SYMBOL),]
#head(AvsC_lucas_NA_removed$SYMBOL)
#plotMD(tfit_lucas, column = 2, status = dt_lucas[,2], xlim = c(-8,13), main = "AvsC mean-difference plot")
 
####HEATMAPS####

par(mfrow = c(1,1))
for (i in 1:length(Comparisons_vector)) {
  if (SavePlots == TRUE) {pdf(file = paste(Directory_Plots_RNASeq, "/HeatMap_Lucas", paste0("_", Comparisons_vector[i]), "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                                           filter_type, gsub(".*filter", "", Directory_Plots_RNASeq), if (outliers == TRUE) {print("_outliers-removed")}, ".pdf", sep = ""), width=16.25, height=12.5)}
  
  j <- which(v_lucas$genes$ENSEMBL %in% Ranked_DEGenes_list[[i]]$ENSEMBL[1:100])
             mycol <- colorpanel(1000, "blue", "white", "red")
             heatmap.2(lcpm_lucas[j,], scale = "row", labRow = v_lucas$genes$SYMBOL[j], labCol = sample_names_lucas_final, col = mycol, trace = "none",
                       density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'both', distfun = function(x) as.dist(1-cor(t(x), method = "s")),
                       main = paste0(Comparisons_vector[i], " Heat Map"))
             if (SavePlots == TRUE) {dev.off()}
}          
  
#BvsC_topgenes <- BvsC_lucas_NA_removed$ENSEMBL[1:80]
#i <- which(v_lucas$genes$ENSEMBL %in% BvsC_topgenes)
#mycol <- colorpanel(1000, "blue", "white", "red")
#heatmap.2(lcpm_lucas[i,], scale = "row", labRow = v_lucas$genes$SYMBOL[i], labCol = sample_names_lucas_final, col = mycol, trace = "none",
#         density.info = "none", margin = c(8,6), lhei = c(2,10), dendrogram = 'both', distfun = function(x) as.dist(1-cor(t(x), method = "s")), main = "BvsC HeatMap")

#write.csv(BvsC_lucas_NA_removed, paste("BvsC_lucas_NA_removed", filter_type, "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#write.csv(AvsB_lucas_NA_removed, paste("AvsB_lucas_NA_removed", filter_type, "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
#write.csv(AvsC_lucas_NA_removed, paste("AvsC_lucas_NA_removed", filter_type, "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE), if (outliers == TRUE) {print("_outliers-removed")}, ".csv", sep = ""))
