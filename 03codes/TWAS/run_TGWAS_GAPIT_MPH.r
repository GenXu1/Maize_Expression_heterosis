##module load anaconda
#conda activate my-R
library(Ropt)
library(data.table)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
#Phenotypic Data
myY <- read.table("MPH_599_hybrid.seeding.phenotype.txt", head = TRUE)
#Numerical genotype format
#--------------------A pair of Genotypic Data and map files-------------------------------
gene=read.table("maize-genome-V4.gene.txt", head = TRUE)
myGD <- fread("Heterosis.599Hybrids.log_transform_imputed_02.txt",data.table=F,sep="\t")
myGM <- gene[,1:3]
ov=intersect(myGM[,1],colnames(myGD)[-1])
myGD=cbind(myGD[,1],myGD[,colnames(myGD)%in%ov])
myGM <-myGM[myGM[,1]%in%colnames(myGD)[-1],]
colnames(myGM)=c("Name","Chromosome","Position")
myGM <-myGM[order(myGM[,2],myGM[,3]),]
colnames(myGD)[1]="Taxa"
d1=data.frame(Name=colnames(myGD)[-1],id=2:ncol(myGD))
d2=merge(myGM,d1,by="Name")
d2=d2[order(d2[,2],d2[,3]),]
myGD=myGD[,c(1,d2[,4])]
#Kinship matrix
myKI <- read.table("599hybrids.kinship2.txt", head = FALSE)
#covaraite variables (such as population structure represented by Q matrix or PC)
myCV <- read.table("599hybrids_PCA_PEER_miss01.txt", head = TRUE)
all(myGD[,1]==myY[,1])
all(myGD[,1]==myKI[,1])
all(myGM[,1]==colnames(myGD)[-1])

#Run TWAS with CMLM
myGAPIT <- GAPIT(Y=myY,
                 GD=myGD, GM=myGM,
                 KI=myKI,CV=myCV,Multiple_analysis=TRUE,
                 kinship.cluster=("average"),
                 kinship.group=("Mean"),
                 group.from=30,
                 group.to=100000,
                 group.by=10,
                 SNP.MAF=0,
                 model=c("MLM","CMLM","FarmCPU","Blink")
)

colnames(myGM)=c("SNP","Chr","Pos")
fwrite(myGM,file="map.txt",row.names=F,col.names=T,quote=F,sep="\t")
fwrite(t(myGD[,-1]),file="Numeric.txt",row.names=F,col.names=F,quote=F,sep="\t")
