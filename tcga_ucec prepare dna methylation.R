library(TCGAbiolinks)
library(SummarizedExperiment)

### DNA methylation
query.ucec.dnamethy <- GDCquery(project = "TCGA-UCEC",
													data.category = "DNA Methylation",
													legacy = FALSE,
													platform = "Illumina Human Methylation 450")

GDCdownload(query = query.ucec.dnamethy, directory = "./examples/TCGA")

methy.UCEC <- GDCprepare(query = query.ucec.dnamethy,
												directory = "./examples/TCGA",
												summarizedExperiment = FALSE)
												
methy.UCEC.sumy <- GDCprepare(query = query.ucec.dnamethy,
												directory = "./examples/TCGA",
												summarizedExperiment = TRUE)

methy.measure <- assay(methy.UCEC.sumy)

methy.compound.info <- rowRanges(methy.UCEC.sumy)

methy.patient.info <- colData(methy.UCEC.sumy)


save.image(file = "./results/examples/Rdata/UCEC/dnamethy_UCEC.Rdata")


