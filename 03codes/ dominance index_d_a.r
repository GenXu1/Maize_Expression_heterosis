library(data.table)
snp1=fread("hyb_RNA_hyb_GWAS.SIG.SNPs.recode.vcf",head=T,skip=25,data.table = F)
exp=fread("fpkm.802samp.gr0.consistentHMPline.txt",head=T,data.table = F)
colnames(exp)=gsub("&","/",colnames(exp))
hyb=fread("RNA_Hybrid_GWAS.SIG.Loci.txt",head=T,data.table = F)
res=NULL
for(i in 1:nrow(snp1))
{
  d1=snp1[i,]
 chr=d1[1,1];pos=d1[1,2]
 AA=colnames(d1)[grep("0/0",d1[1,])]
 aa=colnames(d1)[grep("1/1",d1[1,])]
 Aa=colnames(d1)[grep("0/1|1/0",d1[1,])]
 g=hyb[hyb[,2]==chr & hyb[,6]==pos,]

if(length(AA)>0 & length(Aa)>0 & length(aa)>0 & nrow(g)>0)
{
  for(j in 1:nrow(g))
  {
    g1=unlist(strsplit(g[j,1],"_"))[2]
  exp1=exp[which(exp[,1]==g1),]
  if(nrow(exp1)==0){break}
  g_AA=as.numeric(exp1[1,colnames(exp1)%in%AA])
  g_aa=as.numeric(exp1[1,colnames(exp1)%in%aa])
  g_Aa=as.numeric(exp1[1,colnames(exp1)%in%Aa])
  g_AA1=median(g_AA,na.rm=T)
  g_aa1=median(g_aa,na.rm=T)
  g_Aa1=median(g_Aa,na.rm=T)
  m1=0.5*(g_AA1+g_aa1)
  A=abs(0.5*(g_AA1-g_aa1))
  D=g_Aa1-m1
  da=D/A
  
  
  re=c(chr,pos,g1,da,A,D)
  res=rbind(res,re)
  }
  
}
 cat(i,"\n")
}
colnames(res)=c("Chr","Pos","target_Gene","D/A_median","A_median","D_median")
write.table(res,file="RNA_hyb_GWAS.SIG.SNPs.DA_Index_D_A.txt",col.names = T,row.names = F,quote=F,sep="\t")

