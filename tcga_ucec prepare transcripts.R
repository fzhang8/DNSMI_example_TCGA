library(TCGAbiolinks)
library(SummarizedExperiment)


### transcripts
query.ucec.transcripts <- GDCquery(project = "TCGA-UCEC",
															data.category = "Transcriptome Profiling",
															data.type = "Gene Expression Quantification",
															legacy = FALSE,
															workflow.type = "HTSeq - FPKM-UQ")
															
GDCdownload(query = query.ucec.transcripts, directory = "./examples/TCGA")


trans.UCEC <- GDCprepare(query = query.ucec.transcripts,
												directory = "./examples/TCGA",
												summarizedExperiment = FALSE)


trans.UCEC.sumy <- GDCprepare(query = query.ucec.transcripts,
												directory = "./examples/TCGA",
												summarizedExperiment = TRUE)

trans.measure <- assay(trans.UCEC.sumy)

trans.compound.info <- rowRanges(trans.UCEC.sumy)

trans.patient.info <- colData(trans.UCEC.sumy)


save.image(file = "./results/examples/Rdata/UCEC/transcripts_UCEC.Rdata")

