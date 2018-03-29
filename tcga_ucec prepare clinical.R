library(TCGAbiolinks)
library(SummarizedExperiment)



#### clincal
query.ucec.clinical <- GDCquery(project = "TCGA-UCEC",data.category = "Clinical",legacy = FALSE)


GDCdownload(query = query.ucec.clinical,directory = "./examples/TCGA")

ucec.clinical <- GDCprepare_clinic(query = query.ucec.clinical,
																	clinical.info = "patient",
																	directory = "./examples/TCGA")


#### legacy clinical
query.ucec.clinical.legcy <- GDCquery(project = "TCGA-UCEC",data.category = "Clinical",
																			legacy = TRUE, data.type = "Clinical data" ,file.type = "txt")
																			
																			
GDCdownload(query = query.ucec.clinical.legcy, directory = "./examples/TCGA")


ucec.clinical.legacy <- GDCprepare(query = query.ucec.clinical.legcy,directory = "./examples/TCGA")

ucec.clinical.legacy <- ucec.clinical.legacy$clinical_patient_ucec


save.image(file = "./results/examples/Rdata/UCEC/clinical_UCEC.Rdata")

