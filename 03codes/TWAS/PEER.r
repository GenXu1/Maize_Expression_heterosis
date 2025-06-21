#module load anaconda
#conda activate r-peer
library(data.table)
library(peer)
library(impute)
run_peer=function(type,K,rand.seed)
  {
  file.input=paste(type,".txt", sep="")
file.out.pdf=paste(type,"_plotModel.pdf", sep="")
file.out.factors=paste(type,"_factors.txt", sep="")
file.out.weights=paste(type,"_weights.txt", sep="")
file.out.precision=paste(type,"_precision.txt", sep="")
file.out.residuals=paste(type,"_residuals.txt", sep="")
file.out.seed=paste(type,"_seed.txt", sep="")
ex_clean.raw <- fread(file = file.input, head = TRUE, data.table = F,na.str="NA")
ex_clean <- ex_clean.raw[, -1]
ex_clean=impute.knn(as.matrix(ex_clean)) ##impute missing expression data
cat("Finish imputing expression data\n")
ex_clean=ex_clean$data
### create the model object
set.seed(rand.seed)
model = PEER()

### Set 25 factors to model
PEER_setNk(model, K)

### Set expression data
PEER_setPhenoMean(model, as.matrix(ex_clean))
dim(PEER_getPhenoMean(model))

### Train the model, observing convergence
PEER_update(model)
### Plot the posterior variance of the factor weights and convergence diagnostics
pdf(file.out.pdf, width = 8, height = 8)
PEER_plotModel(model)
dev.off()

# name of genes & accessions
genes <- colnames(ex_clean)
samples <- as.character(ex_clean.raw$Taxa)

### write PEER factors
fac <- PEER_getX(model)
fac <- as.data.frame(fac)
row.names(fac) <- samples
colnames(fac) <- paste("PEERfac_", 1:K, sep="")
fwrite(fac, file = file.out.factors, col.names = T, row.names = T, sep = "\t", quote = F)

### write PEER factors weights
wg <- PEER_getW(model)
wg <- as.data.frame(wg)
row.names(wg) <- genes
colnames(wg) <- paste("PEERfac_", 1:K, sep="")
fwrite(wg, file = file.out.weights, col.names = T, row.names = T, sep = "\t", quote = F)

### write precision (inverse variance) of the weights
pre <- PEER_getAlpha(model)
pre <- as.data.frame(pre)
pre <- cbind(paste("PEER", 1:K, sep = ""), pre)
colnames(pre) <- c("factor", "precision")
write.table(pre, file = file.out.precision, col.names = T, row.names = F, sep = "\t", quote = F)

### write residual dataset
res <- PEER_getResiduals(model)
res <- as.data.frame(res)
row.names(res) <- samples
colnames(res) <- genes
fwrite(res, file = file.out.residuals, col.names = T, row.names = T, sep = "\t", quote = F)

# save random seed
write(rand.seed, file = file.out.seed)
}
