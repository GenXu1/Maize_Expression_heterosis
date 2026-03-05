library(data.table)
library(vioplot)
##aldult traits
d1=fread("inputdata/fig1/MPH_597_harvest.phenotype.txt",header = T,data.table = F)
d1=d1[,c(1,grep("Overall",colnames(d1)))]
colnames(d1)=gsub("_Overall","",colnames(d1))
x=d1[,-1]
id=order(apply(x, 2, median,na.rm=T))+1
d1=d1[,c(1,id)]

d2=fread("inputdata/fig1/MPH_718_hybrid.seeding.phenotype.txt",header = T,data.table = F)
x=d2[,-1]
id=order(apply(x, 2, median,na.rm=T))+1
d2=d2[,c(1,id)]

d3=fread("inputdata/fig1/MPH_ratio_FPKM_F1_P_grp1_miss03_b.txt",header = T,data.table = F)##this is the MPH ratio of gene expression
h=fread("inputdata/fig1/House_Keeping_gene.txt", header=T,data.table=F)
h1=h[h[,3]=="Lin F et al., 2014",]
d4=d3[,c(1,match(h1[,2],colnames(d3)))]
x=d4[,-1]
id=order(apply(x, 2, median,na.rm=T))+1
d5=d4[,c(1,id)]
d6=merge(d1,d2,by=1)
d6=merge(d6,d5,by=1)
cols=c("darkgreen","goldenrod1")
cols=adjustcolor(cols,alpha.f = 0.6)
col2=c(rep(cols[1],5),rep(cols[2],7))
pdf("MPH_distribution_adult_seedling_traits_V3.pdf",height = 2.6,width = 3.7)
par(mfrow=c(1,1),mar=c(2,3,2,2))
vioplot(d6[,7:13]*100,col=col2[6:12],border =col2[6:12],ylab="MPH (100%)",names=colnames(d6)[7:13],cex.axis = 1,cex.lab = 1,main="",las=1,ylim=c(-50,300),axes=F,horizontal =F,xlab="",cex=0.4,xlim=c(0.7,7.2))
abline(h=0,col="black",lty=2)
#grid()
#abline(v=5.5)
#box()
#legend("topright",c("Adult traits","Seedling traits","Housekeeping genes"),col=cols,pch=15,bty="n",cex=1.2)
dev.off()

####plot trait density
inf=fread("inputdata/fig1/lines_info.txt", header=F,data.table=F)
p1=unique(inf[,2]);p2=unique(inf[,3]);f1=inf[,1];p=unique(p1,p2)
hyb1=fread("inputdata/fig1/06.blup.hybrid.pheno.harvest.final.txt", header=T,data.table=F)[,c(1,4,6,3,5,2)]
inb1=fread("inputdata/fig1/06.blup.inbred.pheno.harvest.final.txt", header=T,data.table=F)[,c(1,4,6,3,5,2)]

hyb2=fread("inputdata/fig1/blup.hybrid.pheno.seedling2.txt", header=T,data.table=F)[,c(1,4,7,5,3,8,6,2)]
inb2=fread("inputdata/fig1/blup.inbred.pheno.seedling2.txt", header=T,data.table=F)[,c(1,4,7,5,3,8,6,2)]
colnames(hyb1)[1]=colnames(hyb2)[1]
colnames(inb1)[1]=colnames(inb2)[1]
hyb=merge(hyb1,hyb2,by=1)
inb=merge(inb1,inb2,by=1)
cols=c("darkgreen","goldenrod1")
cols=adjustcolor(cols,alpha.f = 0.8)
pdf("Adult_sleedling_density_case_v2.pdf",width=8,height=3)
par(mfrow=c(1,1),mar=c(3,1,2,1))
i=13
  p1a=inb[inb[,1]%in%p1,i]
  p2a=inb[inb[,1]%in%p2,i]
  f1a=hyb[,i]
  p1a=p1a[!is.na(p1a)]
  p2a=p2a[!is.na(p2a)]
  f1a=f1a[!is.na(f1a)]
  plot(0,0,xlim=range(c(p1a,f1a)),ylim=range(c(density(p1a)$y,density(f1a)$y)),xlab="",ylab="",axes=F,pch=16,col="white",main=colnames(hyb)[i],cex.lab=1)
  axis(1,cex.axis=1.1)
  axis(2,cex.axis=1.1)
  # axis(2,cex.axis=1.1,las=1)
  lines(density(c(p1a,p2a)),lty=1,col=cols[1],lwd=2.5)
  # lines(density(p2a),lty=1,col=cols[2],lwd=2.5)
  lines(density(f1a),lty=1,col=cols[2],lwd=2.5) 
  if(i==13)
  {
    legend("topright",c("Parents","F1"),col=cols[c(1,2)],lty=1,lwd=2,bty="n",cex=1.2)
  }

dev.off()


###plot molecular traits
d=fread("inputdata/fig1/MPH_gene_summary_Gen_06_19_2025.txt", header=T,data.table=F)
sd=fread("inputdata/fig1/Gene_expression_MPH_sd.txt", header=T,data.table=F)
d$col="gray"
pos=d1[d1[,2]>0,]
neg=d1[d1[,2]<0,]
d[d[,7]=="Positive",8]="goldenrod1"
d[d[,7]=="Negative",8]="darkgreen"
d=d[order(d[,2]),]
d$ID=1:nrow(d)
h=fread("inputdata/fig1/House_Keeping_gene.txt", header=T,data.table=F)
h1=h[h[,3]=="Lin F et al., 2014",]
d2=d[d[,1]%in%h1[,2],]
sd1=sd[sd[,1]%in%h1[,2],]
sd2=merge(d,sd1,by=1)
sd2=sd2[order(sd2[,2]),c(1,2,9,10)]
pdf("MPH_distribution_mocular_traits.pdf",height = 2.8,width = 5.1)
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(d[,2],pch=16,col=adjustcolor(d[,8],0.8),cex=0.5,ylab="MPH (100%)",main="",ylim=c(-100,150),axes=F,xlab="Gene order")
abline(h=0,lty=2,col="red")
axis(2,las=1)
axis(1)
segments(x0=sd2[,3],y0=sd2[,2]+sd2[,4], x1 = sd2[,3], y1 = sd2[,2]-sd2[,4],col=adjustcolor("black",0.6))
segments(x0=sd2[,3]-60,y0=sd2[,2]+sd2[,4], x1 = sd2[,3]+60, y1 = sd2[,2]+sd2[,4],col=adjustcolor("black",0.6))
segments(x0=sd2[,3]-60,y0=sd2[,2]-sd2[,4], x1 = sd2[,3]+60, y1 = sd2[,2]-sd2[,4],col=adjustcolor("black",0.6))
points(d2$ID,d2[,2],pch=16,col="blue",cex=0.6)
#grid()
#legend("topleft",c("Positive","Negative","n.s.","Housekeeping"),pch=16,cex=1.0,col=c("goldenrod1","darkgreen","gray","blue"),ncol = 1,bty="n")
dev.off()

######################

inf=fread("inputdata/fig1/lines_info.txt", header=F,data.table=F)
p1=unique(inf[,2]);p2=unique(inf[,3]);f1=inf[,1];p=unique(p1,p2)
g=fread("inputdata/fig1/bx.txt", header=T,data.table=F)[c(9:10,11,14,1:8,13),]
d=fread("inputdata/fig1/fpkm.802samp.gr0.consistentHMPline_v2.txt", header=T,data.table=F)
re=data.frame(ID=1:599)
pa1=NULL
for(i in g[,2])
{
  d1=d[d[,1]==i,]
  hy=d1[,grep("/",colnames(d1))]
  pa=d1[,colnames(d1)%in%p]
  hy=data.frame(ID=1:nrow(t(hy)),t(hy))
  pa=data.frame(ID=1:nrow(t(pa)),t(pa))
  colnames(hy)[2]=paste("hyb",i,sep = "_")
  colnames(pa)[2]=paste("par",i,sep = "_")
  re=merge(re,pa,by=1,all.x=T)
  re=merge(re,hy,by=1,all.x=T)
pav=wilcox.test(hy[,2],pa[,2])$p.value
r=c(i,pav)
pa1=rbind(pa1,r)
}

re[,-1]=apply(re[,-1],2,function(x) {x=log2(x+1)})
cols=c("#3ECDF9","mediumorchid1")
cols=adjustcolor(cols,alpha.f = 0.6)
pdf("bx_FPKM2.pdf",height = 2.5,width = 10)
par(mfrow=c(1,1),mar=c(2,3,2,2))
vioplot(re[,-1],col=cols,border =cols,ylab="log2(FPKM)",names=colnames(re)[-1],cex.axis = 1,cex.lab = 1,main="",las=1,ylim=c(0,10),axes=F,horizontal =F,xlab="",cex=0.4)
#grid()
#abline(v=5.5)
#box()
legend("topright",c("Parents","F1"),col=cols,pch=15,bty="n",cex=1,ncol=2)
dev.off()
