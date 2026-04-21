library(data.table)
library(vioplot)
######plot manhatan
ch=fread("inputdata/fig3/chrom.txt", header=T,data.table=F)
ge=fread("inputdata/fig3/maize-genome-V4.gene.txt",header = T,sep="\t",data.table=F)
gen=NULL
for (k in 1:10)
{
  sub=subset(ge,ge[,2]==k)
  sub[,3]=sub[,3]+ch[k,3]
  sub[,4]=sub[,4]+ch[k,3]
  gen=rbind(gen,sub)
}
gen[,3]=gen[,3]/1e6;gen[,4]=gen[,4]/1e6;

i="Zm00001d003817"
  g1=gen[gen[,1]==i,]
  inbf=paste("inputdata/fig3/Inb_",i,".001sig.txt",sep = "")
  hybf=paste("inputdata/fig3/hyb_",i,".001sig.txt",sep = "")
  MPHf=paste("inputdata/fig3/MPH_",i,".001sig.txt",sep = "")
  inb=fread(inbf,header = T,sep="\t",data.table=F)
  hyb=fread(hybf,header = T,sep="\t",data.table=F)
  MPH=fread(MPHf,header = T,sep="\t",data.table=F)
  inb=na.omit(inb);hyb=na.omit(hyb);MPH=na.omit(MPH);
  inb1=inb[,c(1,2,3,9)]
  va1=max(-log10(inb1[,4]))+1
  hyb1=hyb[,c(1,2,3,9)]
  va2=max(-log10(hyb1[,4]))+1
  MPH1=MPH[,c(1,2,3,9)]
  va3=max(-log10(MPH1[,4]))+1
  
  inb2=NULL
  for (k in 1:10)
  {
    sub=subset(inb1,inb1[,1]==k)
    sub[,3]=sub[,3]+ch[k,3]
    inb2=rbind(inb2,sub)
  }
  
  hyb2=NULL
  for (k in 1:10)
  {
    sub=subset(hyb1,hyb1[,1]==k)
    sub[,3]=sub[,3]+ch[k,3]
    hyb2=rbind(hyb2,sub)
  }
  
  MPH2=NULL
  for (k in 1:10)
  {
    sub=subset(MPH1,MPH1[,1]==k)
    sub[,3]=sub[,3]+ch[k,3]
    MPH2=rbind(MPH2,sub)
  }
  
  inb2[,3]=inb2[,3]/1e6
  hyb2[,3]=hyb2[,3]/1e6
  MPH2[,3]=MPH2[,3]/1e6
  col1=ifelse(inb2[,1]%%2==0,"gray81", "gray62")
  col2=ifelse(hyb2[,1]%%2==0,"gray81", "gray62")
  col3=ifelse(MPH2[,1]%%2==0,"gray81", "gray62")
inb3=fread("inputdata/fig3/Zm00001d003817.inb_sig_SNP.txt",header = F,sep="\t",data.table=F)
hyb3=fread("inputdata/fig3/Zm00001d003817.hyb_sig_SNP.txt",header = F,sep="\t",data.table=F)
mph3=fread("inputdata/fig3/Zm00001d003817.MPH_sig_SNP.txt",header = F,sep="\t",data.table=F)
inb3=inb2[inb2[,2]%in%inb3[,2],]
hyb3=hyb2[hyb2[,2]%in%hyb3[,2],]
mph3=MPH2[MPH2[,2]%in%mph3[,2],]

  tiff(paste(i,".man_V2.tiff",sep=""),width = 210, height =70,res=600,units="mm" )
  par(mfrow=c(4,1),mar=c(0,4,0,2))
  plot(inb2[,3],-log10(inb2[,4]),col=col1,pch=16,cex=0.8,bty="l",xlim=c(0,2100),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)",main="")
  abline(h=-log10(1/339202),lty=2)
  axis(2,las=2,tck=-.03,cex.axis=0.8,at=c(3,5,7))
  points(inb3[,3],-log10(inb3[,4]),col=adjustcolor("orange",0.4),pch=16,cex=0.8)
    abline(v=g1[,3],lwd=2,col=adjustcolor("red",alpha.f = 0.4)) 
  box() 
  
  plot(hyb2[,3],-log10(hyb2[,4]),col=col2,pch=16,cex=0.8,bty="l",xlim=c(0,2100),axes=F,cex.lab=0.6,xlab="",ylab="-log10(P)",main="")
  points(hyb3[,3],-log10(hyb3[,4]),col=adjustcolor("orange",0.4),pch=16,cex=0.8)
  abline(h=-log10(1/183864),lty=2)
 axis(2,las=2,tck=-.03,cex.axis=0.8,at=c(3,5,7))

   abline(v=g1[,3],lwd=2,col=adjustcolor("red",alpha.f = 0.4)) 
  box() 
  
  plot(MPH2[,3],-log10(MPH2[,4]),col=col3,pch=16,cex=0.8,bty="l",xlim=c(0,2100),axes=F,cex.lab=1.2,xlab="Chromosome",ylab="-log10(P)",main="")
  points(mph3[,3],-log10(mph3[,4]),col=adjustcolor("orange",0.4),pch=16,cex=0.8)
  abline(h=-log10(1/183864),lty=2)
  axis(1,at=c(153.5208585,429.262855,669.31791,910.6491295,1146.097552,1345.065257,1523.272613,1705.024703,1875.470912,2030.84696),labels=1:10,tck=-0.03,cex.axis=1.2)
  axis(2,las=2,tck=-.03,cex.axis=0.8,at=c(3,5,7))
  abline(v=g1[,3],lwd=2,col=adjustcolor("red",alpha.f = 0.4)) 
  box()
  dev.off()

##plot pie 
inb=fread("RNA_Inbred_GWAS.SIG.Loci.leading_SNP_feature_5k.txt", header=T,data.table=F)
hyb=fread("RNA_hyb_GWAS.SIG.Loci.leading_SNP_feature_5k.txt", header=T,data.table=F)
MPH=fread("RNA_MPH_GWAS.SIG.Loci.leading_SNP_feature_5k.txt", header=T,data.table=F)
inb$type2=inb$Type
hyb$type2=hyb$Type
MPH$type2=MPH$Type

inb[inb[,9]=="Trans" & inb[,2]==inb[,4],9]="Trans1"
inb[inb[,9]=="Trans" & inb[,2]!=inb[,4],9]="Trans2"
hyb[hyb[,9]=="Trans" & hyb[,2]==hyb[,4],9]="Trans1"
hyb[hyb[,9]=="Trans" & hyb[,2]!=hyb[,4],9]="Trans2"
MPH[MPH[,9]=="Trans" & MPH[,2]==MPH[,4],9]="Trans1"
MPH[MPH[,9]=="Trans" & MPH[,2]!=MPH[,4],9]="Trans2"


inb1=as.data.frame(table(inb[,9]))[c(6,1,3,7,2,4,5),]
hyb1=as.data.frame(table(hyb[,9]))[c(6,1,3,7,2,4,5),]
MPH1=as.data.frame(table(MPH[,9]))[c(6,1,3,7,2,4,5),]
all(inb1[,1]==hyb1[,1]) # check if the order of categories is the same
all(inb1[,1]==MPH1[,1]) # check if the order of categories is the same
pc1=round(inb1[,2]*100/sum(inb1[,2]),2)
lb1=paste(pc1,"%",sep="");names(inb1[,2])=lb1
pc2=round(hyb1[,2]*100/sum(hyb1[,2]),2)
lb2=paste(pc2,"%",sep="");names(hyb1[,2])=lb2
pc3=round(MPH1[,2]*100/sum(MPH1[,2]),2)
lb3=paste(pc3,"%",sep="");names(MPH1[,2])=lb3
color=c("#1b9e77", "#d95f02", "#7570b3",
        "#e7298a", "#66a61e", "#E6AB02", "#F2C50C")

pdf("Leading_SNP_feature_1M.pdf",height = 2.5,width = 8.6)
par(mfrow=c(1,3),mar=c(2,3,2,2))
pie(inb1[,2],labels = lb1,main="Inbred",radius=1,clockwise =F,col=color,border ="white")
legend("topleft",c("Up5k","Down5K","Gene Body","Up5k_1.5M","Down5k_1.5M","Trans1","Trans2"),pch=15,
       col=color,ncol =6,bty = "n")
pie(hyb1[,2],labels = lb2,main="Hybrid",radius=0.7,clockwise =F,col=color,border ="white")
pie(MPH1[,2],labels = lb3,main="MPH",radius=0.4,clockwise =F,col=color,border ="white")
dev.off()

####plot permutation results
######plot cis
ob=fread("inputdata/fig3/res_permutation/observed_eQTL_ovlap_with_feature.txt", header=T,data.table=F)
inb_cis_k4=fread("inputdata/fig3/res_permutation/Inbred_cis-H3K4me3_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
inb_cis_k27=fread("inputdata/fig3/res_permutation/Inbred_cis-H3K27_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
hyb_cis_k4=fread("inputdata/fig3/res_permutation/Hybrid_cis-H3K4me3_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
hyb_cis_k27=fread("inputdata/fig3/res_permutation/Hybrid_cis-H3K27_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
mph_cis_k4=fread("inputdata/fig3/res_permutation/MPH_cis-H3K4me3_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
mph_cis_k27=fread("inputdata/fig3/res_permutation/MPH_cis-H3K27_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)

inb_trans_k4=fread("inputdata/fig3/res_permutation/Inbred_trans-H3K4me3_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
inb_trans_k27=fread("inputdata/fig3/res_permutation/Inbred_trans-H3K27_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
hyb_trans_k4=fread("inputdata/fig3/res_permutation/Hybrid_trans-H3K4me3_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
hyb_trans_k27=fread("inputdata/fig3/res_permutation/Hybrid_trans-H3K27_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
mph_trans_k4=fread("inputdata/fig3/res_permutation/MPH_trans-H3K4me3_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)
mph_trans_k27=fread("inputdata/fig3/res_permutation/MPH_trans-H3K27_peaks_Shoot-eQTL_overlap_sampled_feature.txt", header=T,data.table=F)

col=c("goldenrod1","darkgreen","darkorchid1")
col=adjustcolor(col,alpha.f = 0.6)
ob=ob[grep("H3K27_peaks|H3K4me3_peaks",ob[,2]),]
 ob$f1="K27"
 ob[grep("H3K4me3",ob[,2]),5]="K4"
 ob$type="hybrid"
 ob[grep("Inbred",ob[,1]),6]="inbred"
 ob[grep("MPH",ob[,1]),6]="MPH"
 ob$type2="cis"
 ob[grep("trans",ob[,1]),7]="trans"
 ob[,4]=log10(ob[,4])

df <- data.frame(
     value = c(log10(inb_cis_k4[,3]), log10(inb_cis_k27[,3]),log10(hyb_cis_k4[,3]), log10(hyb_cis_k27[,3]),
                             log10(mph_cis_k4[,3]), log10(mph_cis_k27[,3]), 
                             log10(inb_trans_k4[,3]), log10(inb_trans_k27[,3]),log10(hyb_trans_k4[,3]), log10(hyb_trans_k27[,3]),
                            log10(mph_trans_k4[,3]), log10(mph_trans_k27[,3])
                         ),
    
   f1 = rep(c("K4", "K27","K4", "K27","K4", "K27","K4", "K27","K4", "K27","K4", "K27"), times = rep(1000,12)),
   type = rep(c("inbred","inbred","hybrid","hybrid","MPH","MPH","inbred","inbred","hybrid","hybrid","MPH","MPH"), times = rep(1000,12)),
  type2= rep(c(rep("cis",6),rep("trans",6)), times = rep(1000,12))
   )

library(ggplot2)
library(ggbreak)
library(patchwork)

# 确保 df 和 ob 的 factor 顺序
df$type <- factor(df$type, levels = c("inbred", "hybrid", "MPH"))
df$f1 <- factor(df$f1, levels = c("K4", "K27"))
df$type2 <- factor(df$type2, levels = c("cis", "trans"))

ob$type <- factor(ob$type, levels = c("inbred", "hybrid", "MPH"))
ob$f1 <- factor(ob$f1, levels = c("K4", "K27"))
ob$type2 <- factor(ob$type2, levels = c("cis", "trans"))

# 分别提取 cis 和 trans 数据
cis_df <- df[df$type2 == "cis", ]
cis_ob <- ob[ob$type2 == "cis", ]
trans_df <- df[df$type2 == "trans", ]
trans_ob <- ob[ob$type2 == "trans", ]

p_cis <- ggplot(cis_df, aes(x = type, y = value, fill = f1)) +
  geom_violin(trim = FALSE, scale = "width",
              position = position_dodge(0.8),
              alpha = 0.7, width = 0.6,size = 0.2,color = NA) +
  # geom_boxplot(
  #   width = 0.15,
  #   position = position_dodge(0.8),
  #   outlier.size = 0.2,  # 离群点大小
  #   alpha = 0.8, # 透明度
  #   size = 0.2,
  # #  color = "white" # 箱线颜色
  #  )+
  geom_segment(
    data = cis_ob,
    aes(x = as.numeric(type) - 0.2 + as.numeric(f1 == "K27") * 0.4,
        xend = as.numeric(type) - 0.2 + as.numeric(f1 == "K27") * 0.4+0.1 ,
        y = Overlapped_eQTL, yend = Overlapped_eQTL),  # log10(ob[,4]) 是 V4
    color = "red", linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(1.5,2.5),  # 在x=1.5和x=2.5处画线
    color = "grey70",           # 浅灰色
    linewidth = 0.5,            # 中等粗细
    linetype = "solid"          # 实线
  )+
  #facet_wrap(~ f1) +
  scale_fill_manual(values = c("K4" = "#FFC12599", "K27" = "#00640099")) +
#  theme_minimal() +
  scale_y_break(c(2.85, 4.38), space = 0.3, scales = 5.5) +  # 自定义断层位置
  labs(title = "Cis", y = "log10(value)", x = NULL)+
  theme_bw()+
  theme(
    plot.title = element_text(color="black", size=14, face="bold", hjust=0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text = element_text(face = "bold", color = "black", size = 10),
    legend.position = "none",
    # 以下是新增的去除网格线的设置
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)  # 保留外边框
  )


ggsave("permutation_test_overlap_number_cis_eQTL_V4.pdf", p_cis,width= 6 , height= 3.5 , units="in")

############trans
p_trans <- ggplot(trans_df, aes(x = type, y = value, fill = f1)) +
  geom_violin(trim = FALSE, scale = "width",
              position = position_dodge(0.8),
              alpha = 0.7, width = 0.6,size = 0.2,color = NA) +
  # geom_boxplot(
  #   width = 0.15,
  #   position = position_dodge(0.8),
  #   outlier.size = 0.2,  # 离群点大小
  #   alpha = 0.8, # 透明度
  #   size = 0.2,
  #   #  color = "white" # 箱线颜色
  # )+
  geom_segment(
    data = trans_ob,
    aes(x = as.numeric(type) - 0.2 + as.numeric(f1 == "K27") * 0.4,
        xend = as.numeric(type) - 0.2 + as.numeric(f1 == "K27") * 0.4+0.1 ,
        y = Overlapped_eQTL, yend = Overlapped_eQTL),  # log10(ob[,4]) 是 V4
    color = "red", linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(1.5,2.5),  # 在x=1.5和x=2.5处画线
    color = "grey70",           # 浅灰色
    linewidth = 0.5,            # 中等粗细
    linetype = "solid"          # 实线
  )+
  #facet_wrap(~ f1) +
  scale_fill_manual(values = c("K4" = "#FFC12599", "K27" = "#00640099")) +
  #  theme_minimal() +
 # scale_y_break(c(2.85, 4.3), space = 0.3, scales = 5.5) +  # 自定义断层位置
  labs(title = "Trans", y = "log10(value)", x = NULL)+
  theme_bw()+
  theme(
    plot.title = element_text(color="black", size=14, face="bold", hjust=0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text = element_text(face = "bold", color = "black", size = 10),
    legend.position = "none",
    # 以下是新增的去除网格线的设置
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)  # 保留外边框
  )


ggsave("permutation_test_overlap_number_trans_eQTL_V4.pdf", p_trans,width= 6 , height= 3.5 , units="in")

####plot D/A
a0=fread("inputdata/fig3/MPH_gene_summary_Gen_06_19_2025.txt",header = T,sep="\t",data.table=F)[,c(1,2,7)]
a=fread("inputdata/fig3/RNA_hyb_GWAS.SIG.SNPs.cis_trans_DA_Index.txt",header = T,sep="\t",data.table=F)[,-c(7:8)]
a=merge(a,a0,by="Gene")
a=na.omit(a)
a1=a[a[,4]>-2 & a[,4]<2,]
col=c("goldenrod1","darkgreen","gray")
col=adjustcolor(col,alpha.f = 0.6)
cis1=a1[a1[,6]=="cis",]

####################trans
trans1=a1[a1[,6]=="trans",]

pdf("hyb_MPH_cis_trans_eQTL_DA_V2.pdf",width=6,height=5)
par(mfrow=c(2,1),mar=c(2,4,1,1))
####################cis and trans
plot(0,0,ylim=c(0,1.5),xlab="D/A",ylab="",axes=F,pch=16,col="white",main="",cex.lab=1,xlim=c(-2.5,2.5))
axis(1,cex.axis=1.1)
axis(2,cex.axis=1.1,las=1)
lines(density(cis1[,4]),lty=1,col=col[1],lwd=2)
lines(density(trans1[,4]),lty=1,col=col[2],lwd=2)
abline(v=c(-1.25,-0.25,0.25,1.25),lty=2,lwd=1,col=adjustcolor("gray",alpha.f = 0.8))
legend("topleft",c("cis","trans"),col=col[1:2],lty=1,lwd=2,bty = "n")



d0=fread("inputdata/fig3/MPH_gene_summary_Gen_06_19_2025.txt",header = T,sep="\t",data.table=F)[,c(1,2,7)]
d=fread("inputdata/fig3/RNA_MPH_GWAS.SIG.SNPs.cis_trans_DA_Index.txt",header = T,sep="\t",data.table=F)[,-c(7:8)]
d=merge(d,d0,by="Gene")
d=na.omit(d)
d1=d[d[,4]>-2 & d[,4]<2,]
col=c("goldenrod1","darkgreen","gray")
col=adjustcolor(col,alpha.f = 0.6)
cis=d1[d1[,6]=="cis",]

####################trans
trans=d1[d1[,6]=="trans",]

####################cis and trans
plot(0,0,xlim=range(-2.5,2.5),ylim=c(0,0.5),xlab="D/A",ylab="",axes=F,pch=16,col="white",main="",cex.lab=1)
axis(1,cex.axis=1.1)
axis(2,cex.axis=1.1,las=1)
lines(density(cis[,4]),lty=1,col=col[1],lwd=2)
lines(density(trans[,4]),lty=1,col=col[2],lwd=2)
abline(v=c(-1.25,-0.25,0.25,1.25),lty=2,lwd=1,col=adjustcolor("gray",alpha.f = 0.8))
dev.off()