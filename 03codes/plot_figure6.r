library(Ropt)
library(data.table)
library(vioplot)
#d1=fread("ZmR1_Seedling_trait.txt", header=T,data.table=F)
d1=fread("Seedling_R1_heterosis.txt", header=T,data.table=F)
d1[,3:6]=apply(d1[,3:6],2,function(x){x*100})
d2=d1[d1[,1]==1 | d1[,1]==2,]
unique(d2[,2])

#######OE
d2=d1[d1[,1]==5 | d1[,1]==6,]
d3=d2[d2[,2]=="J92/JING724",]
d4=d2[d2[,2]=="J92/OE1",]
d5=d2[d2[,2]=="J92/OE2",]
d6=d2[d2[,2]=="J92/OE3",]

####EMS
d7=d1[d1[,1]==1 | d1[,1]==2,]
d8=d7[d7[,2]=="MO17/R1-EMS1-6(AA)",]
d9=d7[d7[,2]=="MO17/R1-EMS1-11(aa)",]

col=adjustcolor(c("gray","goldenrod1","gray","goldenrod1","goldenrod1","goldenrod1"),alpha.f = 0.8)
pdf("Jing92_Jing724_OE_EMS_seedling_MPH2.pdf",width=8.6,height=4)
par(mfrow=c(2,1),mar=c(2,4,2,2))
vioplot(d8[,3],d9[,3],d3[,3],d4[,3],d5[,3],d6[,3],
        d8[,4],d9[,4],d3[,4],d4[,4],d5[,4],d6[,4],
        d8[,5],d9[,5],d3[,5],d4[,5],d5[,5],d6[,5],
        d8[,6],d9[,6],d3[,6],d4[,6],d5[,6],d6[,6],
        col=col,border = col,las=1,tck=-.03,cex=0.7,xlim=c(1.3,23.7))
abline(v=c(2.5,6.5,8.5,12.5,14.5,18.5,20.5),col="gray")
points(rep(jitter(rep(1,nrow(d8)),factor =3)),d8[,3],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(2,nrow(d9)),factor =2)),d9[,3],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(3,nrow(d3)),factor =1)),d3[,3],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(4,nrow(d4)),factor =1)),d4[,3],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(5,nrow(d5)),factor =1)),d5[,3],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(6,nrow(d6)),factor =1)),d6[,3],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(7,nrow(d8)),factor =1)),d8[,4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(8,nrow(d9)),factor =1)),d9[,4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(9,nrow(d3)),factor =0.6)),d3[,4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(10,nrow(d4)),factor =0.6)),d4[,4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(11,nrow(d5)),factor =0.6)),d5[,4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(12,nrow(d6)),factor =0.6)),d6[,4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(13,nrow(d8)),factor =0.6)),d8[,5],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(14,nrow(d9)),factor =0.6)),d9[,5],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(15,nrow(d3)),factor =0.6)),d3[,5],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(16,nrow(d4)),factor =0.5)),d4[,5],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(17,nrow(d5)),factor =0.5)),d5[,5],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(18,nrow(d6)),factor =0.5)),d6[,5],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(19,nrow(d8)),factor =0.5)),d8[,6],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(20,nrow(d9)),factor =0.4)),d9[,6],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(21,nrow(d3)),factor =0.4)),d3[,6],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(22,nrow(d4)),factor =0.3)),d4[,6],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(23,nrow(d5)),factor =0.3)),d5[,6],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(24,nrow(d6)),factor =0.3)),d6[,6],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
wilcox.test(d3[,3],d4[,3])
wilcox.test(d3[,3],d5[,3])
wilcox.test(d3[,3],d6[,3])
wilcox.test(d8[,3],d9[,3])

wilcox.test(d3[,4],d4[,4])
wilcox.test(d3[,4],d5[,4])
wilcox.test(d3[,4],d6[,4])
wilcox.test(d8[,4],d9[,4])

wilcox.test(d3[,5],d4[,5])
wilcox.test(d3[,5],d5[,5])
wilcox.test(d3[,5],d6[,5])
wilcox.test(d8[,5],d9[,5])

wilcox.test(d3[,6],d4[,6])
wilcox.test(d3[,6],d5[,6])
wilcox.test(d3[,6],d6[,6])
wilcox.test(d8[,6],d9[,6])
###plot adult traits
d=fread("adult_traits/Adult_R1_heterosis.txt", header=T,data.table=F)

tra=c("PH","TSS","MSS","BSS")
d=d[d[,1]%in%tra,]
d=d[d[,4]>(-20),]
d1=d[d[,3]=="MO17/R1-EMS1AA" ,]
d2=d[d[,3]=="MO17/R1-EMS1aa",]
d3=d[d[,3]=="J92/JING724" ,]
d4=d[d[,3]=="J92/OE1",]
  d5=d[d[,3]=="J92/OE2",]  
  d6=d[d[,3]=="J92/OE3",]
  
  vioplot(d1[d1[,1]=="PH",4],d2[d2[,1]=="PH",4],d3[d3[,1]=="PH",4],d4[d4[,1]=="PH",4],d5[d5[,1]=="PH",4],d6[d6[,1]=="PH",4],
          d1[d1[,1]=="TSS",4],d2[d2[,1]=="TSS",4],d3[d3[,1]=="TSS",4],d4[d4[,1]=="TSS",4],d5[d5[,1]=="TSS",4],d6[d6[,1]=="TSS",4],
          d1[d1[,1]=="MSS",4],d2[d2[,1]=="MSS",4],d3[d3[,1]=="MSS",4],d4[d4[,1]=="MSS",4],d5[d5[,1]=="MSS",4],d6[d6[,1]=="MSS",4],
          d1[d1[,1]=="BSS",4],d2[d2[,1]=="BSS",4],d3[d3[,1]=="BSS",4],d4[d4[,1]=="BSS",4],d5[d5[,1]=="BSS",4],d6[d6[,1]=="BSS",4],
         col=col,border = col,las=1,tck=-.03,cex=0.7,ylim=c(-20,100),xlim=c(1.3,23.7))
  abline(v=c(2.5,6.5,8.5,12.5,14.5,18.5,20.5),col="gray")
  
  points(rep(jitter(rep(1,nrow(d1[d1[,1]=="PH",])),factor =3)),d1[d1[,1]=="PH",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(2,nrow(d2[d2[,1]=="PH",])),factor =2)),d2[d2[,1]=="PH",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(3,nrow(d3[d3[,1]=="PH",])),factor =1)),d3[d3[,1]=="PH",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(4,nrow(d4[d4[,1]=="PH",])),factor =1)),d4[d4[,1]=="PH",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
points(rep(jitter(rep(5,nrow(d5[d5[,1]=="PH",])),factor =1)),d5[d5[,1]=="PH",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(6,nrow(d6[d6[,1]=="PH",])),factor =1)),d6[d6[,1]=="PH",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  
  points(rep(jitter(rep(7,nrow(d1[d1[,1]=="TSS",])),factor =1)),d1[d1[,1]=="TSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(8,nrow(d2[d2[,1]=="TSS",])),factor =.6)),d2[d2[,1]=="TSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(9,nrow(d3[d3[,1]=="TSS",])),factor =0.6)),d3[d3[,1]=="TSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(10,nrow(d4[d4[,1]=="TSS",])),factor =0.6)),d4[d4[,1]=="TSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(11,nrow(d5[d5[,1]=="TSS",])),factor =0.6)),d5[d5[,1]=="TSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(12,nrow(d6[d6[,1]=="TSS",])),factor =0.5)),d6[d6[,1]=="TSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  
  
  points(rep(jitter(rep(13,nrow(d1[d1[,1]=="MSS",])),factor =0.5)),d1[d1[,1]=="MSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(14,nrow(d2[d2[,1]=="MSS",])),factor =0.5)),d2[d2[,1]=="MSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(15,nrow(d3[d3[,1]=="MSS",])),factor =0.5)),d3[d3[,1]=="MSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(16,nrow(d4[d4[,1]=="MSS",])),factor =0.4)),d4[d4[,1]=="MSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(17,nrow(d5[d5[,1]=="MSS",])),factor =0.4)),d5[d5[,1]=="MSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(18,nrow(d6[d6[,1]=="MSS",])),factor =0.4)),d6[d6[,1]=="MSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  
  points(rep(jitter(rep(19,nrow(d1[d1[,1]=="BSS",])),factor =0.4)),d1[d1[,1]=="BSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(20,nrow(d2[d2[,1]=="BSS",])),factor =0.4)),d2[d2[,1]=="BSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(21,nrow(d3[d3[,1]=="BSS",])),factor =0.3)),d3[d3[,1]=="BSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(22,nrow(d4[d4[,1]=="BSS",])),factor =0.3)),d4[d4[,1]=="BSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(23,nrow(d5[d5[,1]=="BSS",])),factor =0.3)),d5[d5[,1]=="BSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
  points(rep(jitter(rep(24,nrow(d6[d6[,1]=="BSS",])),factor =0.3)),d6[d6[,1]=="BSS",4],pch=16,cex=0.5,col=adjustcolor("darkred", alpha.f = 0.5))
dev.off()


