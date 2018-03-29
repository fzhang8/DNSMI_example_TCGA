

cat("loading clinical_UCEC.Rdata","\n")
load("./results/examples/Rdata/UCEC/clinical_UCEC.Rdata")

rm(list = c("query.ucec.clinical","query.ucec.clinical.legcy","ucec.clinical.legacy"))


cat("loading dnamethy_UCEC.Rdata","\n")
load("./results/examples/Rdata/UCEC/dnamethy_UCEC.Rdata")


rm(list = c("methy.UCEC","methy.UCEC.sumy","query.ucec.dnamethy"))


cat("loading transcripts_UCEC.Rdata","\n")
load("./results/examples/Rdata/UCEC/transcripts_UCEC.Rdata")


rm(list = c("query.ucec.transcripts","trans.UCEC","trans.UCEC.sumy"))

library(GenomicRanges)
library(stringr)
library(parallel)

source("./function_genW_v5_parallel.R")


#############################################



methy.patient.7bcode <- rownames(methy.patient.info)[methy.patient.info$definition == "Primary solid Tumor"]
names(methy.patient.7bcode) <- methy.patient.info$patient[methy.patient.info$definition == "Primary solid Tumor"]



trans.patient.7bcode <- rownames(trans.patient.info)[trans.patient.info$definition == "Primary solid Tumor"]
names(trans.patient.7bcode) <- trans.patient.info$patient[trans.patient.info$definition == "Primary solid Tumor"]


clinical.patients <- ucec.clinical$bcr_patient_barcode[ucec.clinical$histological_type == "Endometrioid endometrial adenocarcinoma"]

patients.Y.unmissing <- ucec.clinical$bcr_patient_barcode[!is.na(ucec.clinical$pct_tumor_invasion)]


selected.patients <- intersect(intersect(names(methy.patient.7bcode),
																					names(trans.patient.7bcode)),
																intersect(clinical.patients, patients.Y.unmissing))


selected.methy.pt.7code <- methy.patient.7bcode[selected.patients]


selected.trans.pt.7code <- trans.patient.7bcode[selected.patients]



cat("generating y","\n")
y <- ucec.clinical$pct_tumor_invasion[ucec.clinical$bcr_patient_barcode %in% selected.patients]
names(y) <- ucec.clinical$bcr_patient_barcode[ucec.clinical$bcr_patient_barcode %in% selected.patients]
y <- y[selected.patients]



methy.compound.info <- as.data.frame(methy.compound.info)

selected.methy.compounds <- methy.compound.info$Composite.Element.REF[
																		methy.compound.info$Gene_Symbol != "." & 
																		methy.compound.info$Gene_Type != "." & 
																		methy.compound.info$Transcript_ID != "." & 
																		methy.compound.info$CGI_Coordinate != "." &
																		methy.compound.info$Feature_Type != "." & 
																		methy.compound.info$Position_to_TSS != "." &
																		methy.compound.info$seqnames == "chr10"]


cat("generating G","\n")
G <- t(methy.measure[selected.methy.compounds, selected.methy.pt.7code])

G <- log2(G/(1 - G)) # transfer to M value from beta value


G <- G[,sapply(X = seq(1,ncol(G)),FUN = function(col){all(!is.na(G[,col]))})]

G <- G[,sapply(X = seq(1,ncol(G)),FUN = function(col){!all(G[,col] == 0)})]








########################################

trans.compound.info <- as.data.frame(trans.compound.info)

selected.trans.ensembles <- trans.compound.info$ensembl_gene_id[trans.compound.info$seqnames == "chr10"]



cat("generating X","\n")
X <- t(trans.measure[selected.trans.ensembles, selected.trans.pt.7code])



X <- X[,sapply(X = seq(1,ncol(X)),FUN = function(col){all(!is.na(X[,col]))})]

X <- X[,sapply(X = seq(1,ncol(X)),FUN = function(col){!all(X[,col] == 0)})]







##################################################


subjectsbad <- which(y > 100 | y < 0)

y <- y[-subjectsbad]

X <- X[-subjectsbad,]

G <- G[-subjectsbad,]



obs05 <- list(y = y - mean(y),
						transcripts = scale(X, scale=FALSE), # only centered
						genes = scale(G, scale=FALSE))			# only centered
						

####################################################


cat("generating w...","\n")
cat("Setting parallel cores...")
		setDefaultCluster(NULL)
		setDefaultCluster(cl <- makePSOCKcluster(detectCores(logical = TRUE) - 1)) 
		cat(length(cl))
		clusterEvalQ(expr = {
  		source("./function_OLS.R")
		})
cat(" Done\n")


w05 <- W_svd_v5_parallel(df = obs05,SVDdecomp = FALSE)$w





save(list = c("obs05","w05"),file = "./results/examples/Rdata/UCEC/obs RDT0.5.Rdata")						

save.image(file = "./results/examples/Rdata/UCEC/GXY prepare RDT0.5.Rdata")

cat("Done","\n")




