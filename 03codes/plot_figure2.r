library(Ropt)
library(data.table)
library(vioplot)
library(basicTrendline)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(grid)
#######plot genetic load vs gene expression
d0=fread("inputdata/fig2/211inbreds.gene.up5k.load.txt", header=T,data.table=F)
d=fread("inputdata/fig2/all_Gene_inb_expression_med_mean.txt", header=T,data.table=F)
d=na.omit(d)
d1=d[d[,2]>=1,]
d2=merge(d0,d1,by="Gene")
# 计算 5% 间隔的分位数作为分组断点（共21个）
breaks <- quantile(d2[, 3], probs = seq(0, 1, by = 0.005), na.rm = TRUE)

# 创建标签 t1 到 t20
labels <- paste0("t", 1:200)

# 对第3列进行分组，并添加到新列 d2$group 中
d2$group <- cut(d2[, 3], breaks = breaks, include.lowest = TRUE, labels = labels)

median_by_group <- aggregate(d2[, 2], by = list(Group = d2$group), FUN = median, na.rm = TRUE)
colnames(median_by_group)[2] <- "Median"
median_by_group$id=gsub("t", "", median_by_group$Group)
median_by_group$id=as.numeric(median_by_group$id)

p1 <- ggplot(median_by_group, aes(x = id, y = Median)) +
  geom_point(color = "darkgreen", alpha = 0.6, size = 1) +  # 使用透明蓝色点
  stat_smooth(method = "lm",
              formula = y ~ x + I(x^2),
              se = T,
              color = "red",
              size = 1,
              linetype = "dashed") +
  labs(x = "Expression rank", y = "Genetic load") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, face = "bold"))    # 坐标轴标题字体大小

pdf("Genetic_load_up5k_expression_rank.pdf",width=4.5,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,2))
p1
dev.off()

###plot expression CV and expression |DMP|
d1=fread("inputdata/fig2/Hyb_gene_expression_CV.txt", header=T,data.table=F)
d2=fread("inputdata/fig2/MPH_gene_summary_Gen_06_19_2025.txt", header=T,data.table=F)[,c(1,2)]
colnames(d2)=c("Gene","MPH")
d3=merge(d1,d2,by="Gene")
d3[,3]=abs(d3[,3])
p1 <- ggplot(d3, aes(x =CV,  y = MPH)) +
  geom_pointdensity(method = "neighbors",size=0.7) +
  scale_color_viridis_c() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              se = FALSE, 
              color = "red", 
              size = 1, 
              linetype = "dashed")  +
  labs(x = "Coefficient of variation", y = "|DMP|") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.text = element_text(size = 14, color = "black"),       # 坐标轴刻度字体大小
    axis.title = element_text(size = 16, face = "bold")) 
p1

d1=fread("inputdata/fig2/Hyb_gene_expression_CV.txt", header=T,data.table=F)
d2=fread("inputdata/fig2/RNA_hyb_SNP_based_heritability.txt", header=T,data.table=F)
d2=d2[,c(1,4)]
colnames(d2)=c("Gene","heritability")
d3=fread("inputdata/fig2/MPH_gene_summary_Gen_06_19_2025.txt", header=T,data.table=F)[,c(1,2)]
colnames(d3)=c("Gene","MPH")
d4=merge(d1,d2,by="Gene")
d5=merge(d4,d3,by="Gene")
d5[, 4]=abs(d5[, 4])

d6=fread("inputdata/fig2/all_Gene_hyb_expression_med_mean.txt", header=T,data.table=F)[,c(1,2)]
colnames(d6)=c("Gene","Expression")
d6[,2]=log2(d6[,2]+1) ##log2 transform the expression value
d5=merge(d5,d6,by="Gene")
# 计算 5% 间隔的分位数作为分组断点（共21个）
breaks <- quantile(d5[, 4], probs = seq(0, 1, by = 0.005), na.rm = TRUE)

# 创建标签 t1 到 t20
labels <- paste0("t", 1:200)

# 对第3列进行分组，并添加到新列 d2$group 中
d5$group <- cut(d5[, 4], breaks = breaks, include.lowest = TRUE, labels = labels)
d5=d5[order(d5$group),]
j=1
re=NULL
for(i in unique(d5$group)) 
{
  d6=d5[d5$group==i,]
  r=cor.test(d6$heritability,d6$CV)$estimate
  r=round(r,2)
  
  r1=cor.test(d6$heritability,d6$Expression)$estimate
  r1=round(r1,2)
  
  a=c(i,r,r1)
  re=rbind(re,a)
}
re=as.data.frame(re)
colnames(re)=c("Group","cor_CV_H2","log2EXP_H2")
re$MPH_rank=as.numeric(gsub("t","",re$Group))
re[,3]=as.numeric(re[,3])
re[,2]=as.numeric(re[,2])
re[,2]=re[,2]^2
re[,3]=re[,3]^2
##plot MPH rank and CV h2 cor
df=re

colnames(df)[2]=c("CV_H2")
scale_factor <- max(df$CV_H2, na.rm = TRUE) / max(df$log2EXP_H2, na.rm = TRUE)

p1=ggplot(df, aes(x = MPH_rank)) +
  # 左边 y 轴：CV_H2
  geom_point(aes(y = CV_H2, color = "CV_H2")) +
  stat_smooth(aes(y = CV_H2, color = "CV_H2"),
              method = "lm", formula = y ~ x, se = TRUE,
              size = 1, linetype = "dashed") +
  
  # 右边 y 轴：log2EXP_H2（需乘以缩放因子）
  geom_point(aes(y = log2EXP_H2 * scale_factor, color = "log2EXP_H2")) +
  stat_smooth(aes(y = log2EXP_H2 * scale_factor, color = "log2EXP_H2"),
              method = "lm", formula = y ~ x, se = TRUE,
              size = 1, linetype = "dashed") +
  
  # 设置主副 y 轴
  scale_y_continuous(
    name = "CV_H2",
    sec.axis = sec_axis(~ . / scale_factor, name = "log2EXP_H2")
  ) +
  
  # 设置颜色
  scale_color_manual(
    name = "Trait",
    values = c("CV_H2" = adjustcolor("darkgreen",alpha.f = 0.6), "log2EXP_H2" = adjustcolor("blue",alpha.f = 0.6))
  ) +
  
  # 设置主题
  theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkgreen"),
    axis.title.y.right = element_text(color = "blue"),
    legend.position = "top"
  ) +
  # 主题调整
  xlab("MPH_rank")+
theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.text = element_text(size = 12, color = "black"),       # 坐标轴刻度字体大小
    axis.title = element_text(size = 12, face = "bold")) 
pdf("Gene_exp_absMPH_rank_H2_CV_Exp_cor_200rank_V2.pdf",height = 4,width = 5)
par(mfrow=c(1,1),mar=c(2,4,1,1))
p1
dev.off()


###plot hub gene and module genes MPH, modules were identified from hybrid expression data
d1=fread("inputdata/fig2/Gene_in_module.txt",header = T,data.table = F)
d1=na.omit(d1)
d1[,3]=d1[,3]*100

d2=fread("inputdata/fig2/hub_gene_all_modules.txt",header = T,data.table = F)
d2=na.omit(d2)
d2[,4]=d2[,4]*100

d1=d1[!d1[,2]%in%d2[,2],]

res=NULL
for(i in unique(d1[,1]))
{
  d3=d1[d1[,1]==i,]
  d4=d2[d2[,1]==i,]
  r=c(i,median(d3[,3]),nrow(d3),median(d4[,4]),nrow(d4))
  res=rbind(res,r)
}
res=res[order(as.numeric(res[,2])),]

ID=1:max(as.numeric(res[,3]))
ID=as.data.frame(ID)
for(i in 1:nrow(res))
{
  m=res[i,1]
  d1a=d1[d1[,1]==m,]
  d2a=d2[d2[,1]==m,]
  d1b=data.frame(ID=1:nrow(d1a),d1a[,3])
  d2b=data.frame(ID=1:nrow(d2a),d2a[,4])
  ID=merge(ID,d1b,by=1,all.x = T)
  colnames(ID)[ncol(ID)]=m
  ID=merge(ID,d2b,by=1,all.x = T)
  colnames(ID)[ncol(ID)]=qq("{m}_hub")
}
col1=adjustcolor(c("gray","goldenrod1"),alpha.f = 0.6)
pdf("MPH_of_Hub_genes2.pdf",height = 3,width = 8.6)
par(mfrow=c(1,1),mar=c(2,4,1,1))
vioplot(ID[,-1],col=col1,ylab="MPH (%)",las=1,cex=0.5,border =col1)
abline(h=0,lty=2)
legend("topleft",c("Nonhub","Hub"),cex=1.0,col=col1,ncol = 1,bty="n",pch=15)

dev.off()

p=NULL
j=1
for(i in seq(2,38,by=2))
{
  d5=ID[,c(1,i,i+1)]
  pa=wilcox.test(d5[,2],d5[,3])$p.value
  r=c(j,colnames(ID)[i],pa)
  p=rbind(p,r)
  j=j+1
}
colnames(p)=c("ID","Module","P-value")
p=as.data.frame(p)
write.table(p,"Nonhub_hub_wilcox_p.txt",col.names = T,row.names = F,quote=F,sep="\t")
p[,3]=as.numeric(p[,3])
p1=p[p[,3]<=0.05,]
x=match(p1[,2],colnames(ID))
x=c(x,x+1)
x=sort(x)
pdf("MPH_of_Hub_genes_V2.pdf",height = 3,width = 5)
par(mfrow=c(1,1),mar=c(2,4,1,1))
vioplot(ID[,x],col=col1,ylab="MPH (%)",las=1,cex=0.5,border =col1,axes=F,xlim=c(1.2,20.1))
abline(h=0,lty=2)
legend("topleft",c("Nonhub","Hub"),cex=1.0,col=col1,ncol = 1,bty="n",pch=15)

dev.off()
x=ID[,x]
colnames(x)[seq(1,20,by=2)]
###plot gene number in each module 
d1=fread("inputdata/fig2/Gene_in_module.txt",header = T,data.table = F)[,-4]
d1=na.omit(d1)
d1[,3]=d1[,3]*100
d2=fread("inputdata/fig2/MPH_gene_summary_Gen_06_19_2025.txt",header = T,data.table = F)[,c(1,7)]
d3=merge(d1,d2,by="Gene")

re=NULL
for(i in unique(d3[,2]))
{
  d4=d3[d3[,2]==i,]
  p=d4[d4[,4]=="positive",]
  n=d4[d4[,4]=="negative",]
  f=d4[d4[,4]=="FALSE",]
  r=c(i,median(d4[,3],na.rm=T),nrow(d4),nrow(p),nrow(n),nrow(f))
  re=rbind(re,r)
}
re=as.data.frame(re)
colnames(re)=c("Module","MPH_med","Gene_N","Pos_N","Neg_N","Non_N")
re1=re[order(re[,2]),]
write.table(re1,file="Module_gene_number_from_exp_hyb.txt",col.names = T,row.names = F,
            sep="\t",quote=F)

library(ggplot2)
library(PieGlyph)
d1=fread("inputdata/fig2/Module_gene_number_from_exp_hyb.txt",header = T,data.table = F)
d1[,3]=log10(d1[,3])
d1=d1[order(d1[,2]),]
d1$id=1:nrow(d1)
col=adjustcolor(c("goldenrod1","darkgreen","gray"),alpha.f = 0.6)
p <- ggplot(data = d1, aes(x = id, y = MPH_med))+
  geom_pie_glyph(aes(radius = Gene_N), 
                 slices = c('Pos_N', 'Neg_N', 'Non_N'), 
                 colour = 'black')+
  scale_fill_manual(values = c("Pos_N" = col[1], "Neg_N" = col[2], "Non_N" = col[3])) +
  labs(x = "Module", y = "MPH (%)") +
  theme(
    axis.text = element_text(size = 14, color = "black"),       # 坐标轴刻度字体大小
    axis.title = element_text(size = 16, face = "bold"),       # 坐标轴标题字体大小
    panel.grid.major = element_line(color = "gray", size = 0.25),  # 主要网格线
    panel.grid.minor = element_line(color = "lightgray", size = 0.1), # 次要网格线
    panel.background = element_blank()  # 移除背景
  )

 # theme_classic()
pdf("Module_gene_number2.pdf",height = 5,width = 8.6)
par(mfrow=c(1,1),mar=c(2,4,1,1))
p
dev.off()

