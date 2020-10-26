

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")


#installing some packages, not sure why all this extra code
#if (!requireNamespace("BiocManager"))
  #install.packages("BiocManager")
#BiocManager::install()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mus.musculus", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.8")

install.packages("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Glimma", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db", version = "3.8")
library("org.Mm.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Homo.sapiens", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")



install.packages("Hmisc")

install.packages("pheatmap")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea", version = "3.8")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")
library("ComplexHeatmap")

install.packages("infotheo")
library("infotheo")

install.packages("Rtsne")
library("Rtsne")

library("edgeR")
library("Glimma")
library("Mus.musculus")
library("RColorBrewer")
library("org.Mm.eg.db")
library("Homo.sapiens")
library(gplots)
library("biomaRt")
install.packages("qpcR")
library("qpcR")

install.packages("circlize")

install.packages("devtools")
library(devtools)
install_github("TheJacksonLaboratory/sRACIPE_dev")
