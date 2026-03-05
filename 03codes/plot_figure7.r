library(Ropt)
library(data.table)
exp=fread("inputdata/fig7/fpkm.allsamples.add2nd.txt", header=T,data.table=F)
exp=exp[apply(exp[,-1],1,mean)>=1,]
g=fread("inputdata/fig7/anthocyanin_pathway.txt", header=T,data.table=F)[,1:2]
lig=fread("inputdata/fig7/lignin_syn_reduced.txt", header=T,data.table=F)
x=fread("inputdata/fig7/MPH_gene_summary_Gen_06_19_2025.txt", header=T,data.table=F)
pos=x[x[,7]=="Positive",]
neg=x[x[,7]=="Negative",]
d=fread("inputdata/fig7/R1_MPH_target_Gene.txt", header=T,data.table=F)
d1=fread("inputdata/fig7/R1_interacted_gene_inb.txt", header=T,data.table=F)
d2=fread("inputdata/fig7/R1_interacted_gene_hyb.txt", header=T,data.table=F)
d1=unique(c(d1[,1],d1[,2]))
d1=d1[d1!="Zm00001d026147"]
d2=unique(c(d2[,1],d2[,2]))
d2=d2[d2!="Zm00001d026147"]
re=fread("inputdata/fig7/ZmR1.gene.FPKM_MPH_V2.txt", header=T,data.table=F)
re=re[re[,1]%in%exp[,1],]
re1=re
re1[,-1]=apply(re1[,-1],2,as.numeric)
re1[,-1]=apply(re1[,-1],2,function(x) {x*100})

JJ=apply(re1[,c(4,8)],1,mean)
JJO1=apply(re1[,c(3,6)],1,mean)
JJO2=apply(re1[,c(2,9)],1,mean)
JJO3=apply(re1[,c(5,7)],1,mean)
re2=data.frame(re1[,1],JJ,JJO1,JJO2,JJO3)
colnames(re2)[1]="Gene"

re3=re2[re2[,1]%in%d[,1],]
re4=re2[re2[,1]%in%g[,2],]
re5=re2[re2[,1]%in%lig[,1],]
ov=intersect(g[,2],lig[,1])
re6=re2[re2[,1]%in%ov,]

pdf("OE_hyb_MPH_V6.pdf",width=8.6,height=2)
par(mfrow=c(1,3),mar=c(2,2,2,2))

####R1 targets OE1
plot(re2[,2],re2[,3],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE1",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)
points(re3[,2],re3[,3],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
cor.test(re2[,2],re2[,3])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(re4[,2],re4[,3],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,2],re5[,3],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,2],re6[,3],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))
g2=re2[re2[,1]=="Zm00001d026147",]

x1=re5[,3]-re5[,2]
n=length(x1[x1>0])/length(x1)
x2=re5[,4]-re5[,2]
n2=length(x2[x2>0])/length(x2)
x3=re5[,5]-re5[,2]
n3=length(x3[x3>0])/length(x3)
points(g2[,2],g2[,3],pch=16,cex=1,col="red")
text(g2[,2],g2[,3],"ZmR1",cex=0.8,pos=4,font=3)

####R1 targets OE2
plot(re2[,2],re2[,4],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE2",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)

cor.test(re2[,2],re2[,4])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(re3[,2],re3[,4],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
points(re4[,2],re4[,4],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,2],re5[,4],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,2],re6[,4],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))

points(g2[,2],g2[,4],pch=16,cex=1,col="red")
text(g2[,2],g2[,4],"ZmR1",cex=0.8,pos=4,font=3)

####R1 targets OE3
plot(re2[,2],re2[,5],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE3",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)

cor.test(re2[,2],re2[,5])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))
points(re3[,2],re3[,5],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
points(re4[,2],re4[,5],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,2],re5[,5],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,2],re6[,5],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))

points(g2[,2],g2[,5],pch=16,cex=1,col="red")
text(g2[,2],g2[,5],"ZmR1",cex=0.8,pos=4,font=3)

dev.off()

