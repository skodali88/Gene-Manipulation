
#script to manipulate gene expression data
#setwd("C:/Users/kodal/Desktop/Bioinformagician")

#download packages
install.packages("dplyr")
install.packages("tidyverse")
if(!requireNamespace("BiocManager",quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("GEOquery") 
#load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

#read in the data
dat <- read.csv(file= "../Bioinformagician/GSE183947_fpkm.csv")
#check dimensions of data frame
dim(dat)

#get metadata
gse <- getGEO(GEO ='GSE183947',GSEMatrix = TRUE)

#got an error: Error in system(cmd, intern = intern, wait = wait | intern, show.output.on.console = wait,  : 
  '/c' not found (made sure latest version of GEOquery is being used)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")

#re-run metadata
gse <- getGEO(GEO ='GSE183947',GSEMatrix = TRUE)

#check the object installed
gse

#get phenodata for first object of first element in that list
metadata <- pData(phenoData(gse[[1]]))
#to look at how the metadata looks like
head(metadata)

#select needed columns
metadata.subset <- select(metadata, c(1,10,11,17))
View(metadata.subset)

metadata %>%
select(1,10,11,17) %>%
head()

#rename the columns using pipe operator function(%>%)
metadata %>%
select(1,10,11,17) %>%
rename(tissue = characteristics_ch1) %>%
rename(metastasis = characteristics_ch1.1) %>%
head()

#mutate the columns using pipe operator(%>%)
metadata %>%
select(1,10,11,17) %>%
rename(tissue = characteristics_ch1) %>%
rename(metastasis = characteristics_ch1.1) %>%
mutate(tissue = gsub("tissue: ","",tissue)) %>%
mutate(metastasis = gsub("metastasis: ","",metastasis)) %>%
head() 

#save to new variable called metadata.modified
metadata.modified <- metadata %>%
select(1,10,11,17) %>%
rename(tissue = characteristics_ch1) %>%
rename(metastasis = characteristics_ch1.1) %>%
mutate(tissue = gsub("tissue: ","",tissue)) %>%
mutate(metastasis = gsub("metastasis: ","",metastasis))
View(metadata.modified)

#see how our gene expression data looks like
head(dat)

#reshaping data our gene expression data which is in a wide format to long format
dat %>%
rename(gene = X) %>%
#-gene indicates i dont want to touch the gene column 
gather(key = 'samples',value = 'FPKM',-gene) %>%
head()

#save to new variable called dat.long
dat.long <- dat %>%
rename(gene = X) %>%
gather(key = 'samples',value = 'FPKM',-gene)


#join dataframes= dat.long + metadata.modified
dat.long %>%
left_join(., metadata.modified, by = c("samples" = "description")) %>%
head()

#saving it to data.long (now we have data of both FPKM and metadata)
dat.long <- dat.long %>%
left_join(., metadata.modified, by = c("samples" = "description"))

#explore data ( I want to extract genes for BRCA1 and BRCA2 and compare expression in tumor vs normal samples)
dat.long %>%
filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
#group by our gene to calculate mean and median exp of each gene in each type of tissue and cal mean and median
group_by(gene, tissue) %>%
summarize(mean_FPKM = mean(FPKM), 
          median_FPKM = median(FPKM)) %>%
#to arrange in ascending order
dat.long %>%
filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
#group by our gene to calculate mean and median exp of each gene in each type of tissue and cal mean and median
group_by(gene, tissue) %>%
summarize(mean_FPKM = mean(FPKM), 
          median_FPKM = median(FPKM)) %>%
#to arrange in ascending order
arrange(mean_FPKM)
#to arrange in descending order
arrange(-mean_FPKM)
View(dat.long)







