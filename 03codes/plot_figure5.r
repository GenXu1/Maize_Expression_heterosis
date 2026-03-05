library(Ropt)
library(data.table)
##trait per se
sig=fread("inputdata/fig5/TWAS_sig_V4.txt", header=T,data.table=F)
sig1=sig[sig[,12]!=0,]
m=sig[sig[,1]=="Hybrid",]
m=m[m[,2]%in%sig1[,2],]
ch=fread("inputdata/fig5/chrom.txt", header=T,data.table=F)
tr=fread("inputdata/fig5/traits.txt", header=F,data.table=F)
col=colours()[c(258,53,96,32,145,28,43)]
col=adjustcolor(col,alpha.f = 0.6)
res=NULL
sig=NULL
for(i in 1:nrow(tr))
{
  tra=tr[i,1]
  inf=qq("inputdata/fig5/GAPIT.Association.GWAS_Results.CMLM.{tra}.csv")
  d=fread(inf, header=T,data.table=F)[,c(2:4,1)]
  ##############
d1=NULL
  for (k in 1:10)
  {
    sub=subset(d,d[,1]==k)
    sub[,2]=sub[,2]+ch[k,3]
    d1=rbind(d1,sub)
  }
  d1$Pos=d1$Pos/1e6
  col1=ifelse(d1[,1]%%2==0,colours()[234],colours()[215])
  col1=adjustcolor(col1,alpha.f = .6)
  d1=data.frame(Trait=rep(tra,nrow(d1)),d1,colour=col1)
  thr=0.05/nrow(d)
  d2=d1[d1[,4]>thr,]
  d3=d1[d1[,4]<=thr,]

  res=rbind(res,d2)
  if(nrow(d3)>0)
  {    d3[,6]=col[i]
  sig=rbind(sig,d3)
  }
}
d4=rbind(res,sig)
thr=-log10(0.05/nrow(d))
d4$cex=0.3
d4[d4[,4]<=(0.05/nrow(d)),7]=0.5
d5=d4[d4[,5]%in%m[,2],]
d5=d5[-log10(d5[,4])>=thr,]
tiff("TWAS_res_exp_trait_per_se_raw_add_gene.tiff",res=600,units = "mm",height = 80,width = 160)
par(mfrow=c(1,1),mar=c(2,4,2,2))
plot(d4[,3],-log10(d4[,4]),col=d4[,6],pch=16,cex=d4[,7],bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="",ylab="",main="",ylim=c(0,12))
abline(h=thr,lty=2)
axis(2,las=2,tck=-.02,cex.axis=1.2,at=c(0,2,4,6,8,10))
axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
#box() 
#legend("top",tr[,1],col=col,pch=16,ncol=7,bty="n",cex=0.5)
points(d5[,3],-log10(d5[,4]),col="red",cex=0.6)
dev.off()
d6=d5[d5[,5]=="Zm00001d026147",]
write.table(d5,"TWAS_sig_exp_trait_per_se_add_gene.txt",quote=F,sep="\t",row.names=F,col.names=T)

###trait heterosis
sig=fread("inputdata/fig5/TWAS_sig_V4.txt", header=T,data.table=F)
sig1=sig[sig[,12]!=0,]
m=sig[sig[,1]=="MPH",]
m=m[m[,2]%in%sig1[,2],]
ch=fread("inputdata/fig5/chrom.txt", header=T,data.table=F)
tr=fread("inputdata/fig5/traits.txt", header=F,data.table=F)
col=colours()[c(258,53,96,32,145,28,43)]
col=adjustcolor(col,alpha.f = 0.6)
res=NULL
sig=NULL
for(i in 1:nrow(tr))
{
  tra=tr[i,1]
  inf=qq("inputdata/fig5/GAPIT.Association.GWAS_Results.CMLM.{tra}.csv")
  d=fread(inf, header=T,data.table=F)[,c(2:4,1)]
  ##############
d1=NULL
  for (k in 1:10)
  {
    sub=subset(d,d[,1]==k)
    sub[,2]=sub[,2]+ch[k,3]
    d1=rbind(d1,sub)
  }
  d1$Pos=d1$Pos/1e6
  col1=ifelse(d1[,1]%%2==0,colours()[234],colours()[215])
  col1=adjustcolor(col1,alpha.f = .6)
  d1=data.frame(Trait=rep(tra,nrow(d1)),d1,colour=col1)
  thr=0.05/nrow(d)
  d2=d1[d1[,4]>thr,]
  d3=d1[d1[,4]<=thr,]

  res=rbind(res,d2)
  if(nrow(d3)>0)
  {    d3[,6]=col[i]
  sig=rbind(sig,d3)
  }
}
d4=rbind(res,sig)
thr=-log10(0.05/nrow(d))
d4$cex=0.3
d4[d4[,4]<=(0.05/nrow(d)),7]=0.5
d5=d4[d4[,5]%in%m[,2],]
d5=d5[-log10(d5[,4])>=thr,]
tiff("TWAS_res_exp_trait_MPH_add_gene.tiff",res=600,units = "mm",height = 40,width = 160)
par(mfrow=c(1,1),mar=c(2,4,2,2))
plot(d4[,3],-log10(d4[,4]),col=d4[,6],pch=16,cex=d4[,7],bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="",ylab="",main="",ylim=c(0,10))
abline(h=thr,lty=2)
axis(2,las=2,tck=-.02,cex.axis=1.2,at=c(0,2,4,6,8,10))
axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
#box() 
points(d5[,3],-log10(d5[,4]),col="red",cex=1)
legend("top",tr[,1],col=col,pch=16,ncol=7,bty="n",cex=0.5)
dev.off()
write.table(d5,"TWAS_sig_exp_trait_MPH_add_gene.txt",quote=F,sep="\t",row.names=F,col.names=T)

##plot LUC signals
g=fread("inputdata/fig5/candidate_genes.txt", header=T,data.table=F)[,1]
d0=fread("inputdata/fig5/LUC_res.txt", header=T,data.table=F)
d=d0[grep("/",d0[,2]),]
re=NULL
for(i in g)
{
  d1=d[grep(i,d[,1]),]
  re=rbind(re,d1)
}
d2=t(re[,-c(1:2)])
colnames(d2)=re[,1]
col=c("gray","goldenrod1")
col=adjustcolor(col,alpha.f = 0.6)
pdf("LUC.pdf",width=8.6,height=2)
par(mfrow=c(1,5),mar=c(4,4,2,2))
re=NULL
for(i in seq(1,9,by=2))
{
  vioplot(d2[,c(i,i+1)], col=col, border =col,ylab="LUC/REN",cex.axis = 1,cex.lab = 1.5,main="",las=1,axes=F,horizontal =F,xlab="",cex=0.6,wex=0.7)
  points(rep(jitter(rep(1,nrow(d2)),factor =2)),d2[,i],pch=16,cex=1,col=adjustcolor("darkred", alpha.f = 0.8))
  points(rep(jitter(rep(2,nrow(d2)),factor =1)),d2[,i+1],pch=16,cex=1,col=adjustcolor("darkred", alpha.f = 0.8))
  p1=wilcox.test(d2[,i],d2[,i+1])$p.value
  r=c(colnames(d2)[i],p1)
  re=rbind(re,r)
}
dev.off()
