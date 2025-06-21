library(data.table)
d=fread("MPH_ratio_FPKM_hybrid_grp_1_miss03_imputed.txt",head=T,data.table=F,sep="\t")
Quantile<- apply(d[,-1],2,
                 function(A){
                   A1=sort(A)
                   si=ceiling(length(A1)*0.05);
                   bi=ceiling(length(A1)*0.95);
                   min_x=A1[si];
                   max_x=A1[bi];
                 out<-(2*(A-min_x)/(max_x-min_x));
                 out[out>2]<-2;out[out< 0]<- 0;return(out)})
d1=cbind(d[,1],Quantile)
d1=as.data.frame(d1)
colnames(d1)[1]="Taxa"

fwrite(d1,file="MPH_ratio_FPKM_hybrid_grp_1_miss03_imputed_02.txt",row.names=F,col.names=T,quote=F,sep="\t")

