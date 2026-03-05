library(data.table)
library(Ropt)
ch=fread("inputdata/fig4/chr_length_v4.txt", header=T,data.table=F)
g=fread("inputdata/fig4/maize-genome-V4.gene.txt", header=T,data.table=F)
g1=NULL
for (k in 1:10)
{
  sub=subset(g,g[,2]==k)
  sub[,3]=sub[,3]+ch[k,3]
  sub[,4]=sub[,4]+ch[k,3]
  g1=rbind(g1,sub)
}
g1[,3]=g1[,3]/1e6;g1[,4]=g1[,4]/1e6
g2=g1[g1[,1]%in%c("Zm00001d025281","Zm00001d039267","Zm00001d048623"),]
N=c(174,521,96)##number of target gene
g3=data.frame(g2,N)

ch=fread("inputdata/fig4/chr_length_v4.txt", header=T,data.table=F)
pos=fread("inputdata/fig4/Pos_MPH_gene_enrichment_in_eQTL_hotspot.txt", header=T,data.table=F)
pos1=pos[pos[,11]<=0.01,]
mk0=qq("{pos1[,1]}-{pos1[,2]}")

neg=fread("inputdata/fig4/Neg_MPH_gene_enrichment_in_eQTL_hotspot.txt", header=T,data.table=F)
neg1=neg[neg[,11]<=0.01,]
mk1=qq("{neg1[,1]}-{neg1[,2]}")

d=fread("inputdata/fig4/RNA_MPH_GWAS.SIG.Loci.eQTL_number.3sig.hotspot.eQTL.grp100_genes.merged.txt", header=T,data.table=F)
mk=qq("{d[,1]}-{d[,2]}")
d=data.frame(d,mk)
pos=(d[,3]+d[,2])/2
d=cbind(d[,1],pos,d[,7],d[,9])
d=as.data.frame(d)
colnames(d)=c("Chr","Pos","Gene_N","Loci")
col=colours()[228]
col=adjustcolor(col,alpha.f = 0.6)
d=data.frame(d,col)
d[,2]=as.numeric(d[,2])
d[,3]=as.numeric(d[,3])
dd=NULL
for (k in 1:10)
{
  sub=subset(d,d[,1]==k)
  sub[,2]=sub[,2]+ch[k,3]
  dd=rbind(dd,sub)
}
dd$Pos=dd$Pos/1e6
pos2=dd[dd[,4]%in%mk0,]
neg2=dd[dd[,4]%in%mk1,]
dd=dd[dd[,2]<=2100,]
pdf("merged_eQTL_hotspot_exp_MPH_highlight_enrich_pos_neg_gene_V2.pdf",height = 2.5,width =8.6)
  par(mar=c(4,4,1,2),mfrow=c(1,1))
 plot(dd[,2],dd[,3],col=col,pch=16,bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="",ylab="#Gene",ylim=c(0,2000),cex=.6)
 axis(2,las=2,tck=-.03,cex.axis=1.2)
 axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
abline(v=c(ch[,3]/1e6,2100),col="gray")
 box()  
 #abline(h=90,col="red",lty=2)
 points(pos2[,2],pos2[,3],col=adjustcolor("goldenrod1",alpha.f = 0.8),pch=16,cex=.8)
 points(neg2[,2],neg2[,3],col=adjustcolor("darkgreen",alpha.f = 0.8),pch=16,cex=.8)
 points(g3[,3],g3[,6],col=colours()[32],pch=25,bg=colours()[32],cex=.8)
 dev.off()

###plot case
d0=fread("inputdata/fig4/MPH_gene_summary_Gen_06_19_2025.txt", header=T,data.table=F)[,c(1,2,5,7)]
g23=fread("inputdata/fig4/eQTL_hotspot_Zm00001d048623_target_genes.txt", header=F,data.table=F)
g67=fread("inputdata/fig4/eQTL_hotspot_Zm00001d039267_target_genes.txt", header=F,data.table=F)
g81=fread("inputdata/fig4/eQTL_hotspot_Zm00001d025281_target_genes.txt", header=F,data.table=F)
colnames(g23)=colnames(g67)=colnames(g81)="Gene"
g23=merge(g23,d0,by="Gene")
g67=merge(g67,d0,by="Gene")
g81=merge(g81,d0,by="Gene")

ch=fread("inputdata/fig4/chr_length_v4.txt", header=T,data.table=F)
g=fread("inputdata/fig4/maize-genome-V4.gene.txt", header=T,data.table=F)
g1=NULL
for (k in 1:10)
{
  sub=subset(g,g[,2]==k)
  sub[,3]=sub[,3]+ch[k,3]
  sub[,4]=sub[,4]+ch[k,3]
  g1=rbind(g1,sub)
}
g1[,3]=g1[,3]/1e6;g1[,4]=g1[,4]/1e6
g23=merge(g23,g1,by="Gene")
g67=merge(g67,g1,by="Gene")
g81=merge(g81,g1,by="Gene")

ch=fread("inputdata/fig4/chr_length_v4.txt", header=T,data.table=F)
gp81=fread("inputdata/fig4/Zm00001d025281_peaks.narrowPeak.genes.txt", header=F,data.table=F)
gp67=fread("inputdata/fig4/Zm00001d039267_peaks.narrowPeak.genes.txt", header=F,data.table=F)
gp23=fread("inputdata/fig4/Zm00001d048623_peaks.narrowPeak.genes.txt", header=F,data.table=F)
gp81=gp81[gp81[,11]>0,]
gp67=gp67[gp67[,11]>0,]
gp23=gp23[gp23[,11]>0,]
ov81=intersect(g81[,1],gp81[,15])
ov67=intersect(g67[,1],gp67[,15])
ov23=intersect(g23[,1],gp23[,15])

g67_p=g67[g67[,4]=="positive",]
g67_n=g67[g67[,4]=="negative",]
g67a=g67[g67[,1]%in%ov67,]

g23_p=g23[g23[,4]=="positive",]
g23_n=g23[g23[,4]=="negative",]
g23a=g23[g23[,1]%in%ov23,]

g81_p=g81[g81[,4]=="positive",]
g81_n=g81[g81[,4]=="negative",]
g81a=g81[g81[,1]%in%ov81,]
pdf("three_case_gene_targets_V3.pdf",height = 2.5,width =10)
  par(mar=c(3,2,1,1),mfrow=c(1,3))
 plot(g67[,6],g67[,2],col=adjustcolor("gray",alpha.f = .6),pch=16,cex=1,bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="",ylab="MPH(%)",ylim=c(-60,140))
 axis(2,las=2,tck=-.03,cex.axis=1.2)
 axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
abline(v=c(ch[,3]/1e6,2100),col="gray")
 box()  
 #abline(h=90,col="red",lty=2)
 points(g67_p[,6],g67_p[,2],col=adjustcolor("goldenrod1",alpha.f = .6),pch=16,cex=1.2)
 points(g67_n[,6],g67_n[,2],col=adjustcolor("darkgreen",alpha.f = .6),pch=16,cex=1.2)
 
 points(g67a[,6],g67a[,2],col=adjustcolor("red",alpha.f = .4),pch=25,bg=adjustcolor("red",alpha.f = .4),cex=0.6)
 abline(h=0,col="black",lty=2)
 
 plot(g23[,6],g23[,2],col=adjustcolor("gray",alpha.f = .6),pch=16,cex=1,bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="",ylab="",ylim=c(-60,140))
 axis(2,las=2,tck=-.03,cex.axis=1.2)
 axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
 abline(v=c(ch[,3]/1e6,2100),col="gray")
 box()  
 #abline(h=90,col="red",lty=2)
 points(g23_p[,6],g23_p[,2],col=adjustcolor("goldenrod1",alpha.f = .6),pch=16,cex=1.2)
 points(g23_n[,6],g23_n[,2],col=adjustcolor("darkgreen",alpha.f = .6),pch=16,cex=1.2)
 points(g23a[,6],g23a[,2],col=adjustcolor("red",alpha.f = .4),pch=25,bg=adjustcolor("red",alpha.f = .4),cex=0.6)
 abline(h=0,col="black",lty=2)
 
 plot(g81[,6],g81[,2],col=adjustcolor("gray",alpha.f = .6),pch=16,cex=1,bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="",ylab="",ylim=c(-60,140))
 axis(2,las=2,tck=-.03,cex.axis=1.2)
 axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
 abline(v=c(ch[,3]/1e6,2100),col="gray")
 box()  
 #abline(h=90,col="red",lty=2)
 points(g81_p[,6],g81_p[,2],col=adjustcolor("goldenrod1",alpha.f = .6),pch=16,cex=1.2)
 points(g81_n[,6],g81_n[,2],col=adjustcolor("darkgreen",alpha.f = .6),pch=16,cex=1.2)
 
 points(g81a[,6],g81a[,2],col=adjustcolor("red",alpha.f = .4),pch=25,bg=adjustcolor("red",alpha.f = .4),cex=0.6)
 abline(h=0,col="black",lty=2)
  dev.off()
