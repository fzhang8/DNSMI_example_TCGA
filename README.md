Required R packages: TCGAbiolinks, SummarizedExperiment, GenomicRanges, stringr, parallel


To make data set RDT1 from extracting data from TCGA:

1, create a folder called Folder2 and paths "Folder2/examples/TCGA/" and "Folder2/results/examples/Rdata/UCEC/".

2, Download all the files to Folder2.

3, In Rstudio, set the working directory to Folder2, and run the following commands:


   source("./tcga_ucec prepare clinical.R")
    
   source("./tcga_ucec prepare dna methylation.R")
    
   source("./tcga_ucec prepare transcripts.R")
    
   Note that it can take quite a long time from hours to days to download the above database depending on your connection speed.

   The raw database is stored under the path Folder2/examples/TCGA, and it is **64.4 GB**
   
   The "GDCdownload" command is just needed once. Once the database is downloaded, that command can be commented out.
   
   The usable .Rdata file is stored under the path Folder2/results/examples/Rdata/UCEC/
   
   clinical_UCEC.Rdata is 135 KB
   
   dnamethy_UCEC.Rdata is **4.29 GB**
   
   transcripts_UCEC.Rdata is 392 MB
   
   So, arrange sufficient storage before downloading.
  
 4, After the downloading is completed, run the command: source("./make GXY V2.R")
 
 The result is stored under the path: Folder2/results/examples/Rdata/UCEC/obs RDT0.5.Rdata, in that .Rdata file "obs05" is the RDT1 dataset and it is a list with 3 elements: y, a vector of 269, transcripts, 269 x 2107, genes, 269 x 10135. "w05" is the NSM matrix, 10135 x 2107, generated from "obs05".
 
 5, ***Notation clarification***
 
 Throughout the code, the letter "w" refers to the meaning of NSM matrix in the manuscript.
