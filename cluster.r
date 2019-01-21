#!/usr/bin/env Rscript
rm(list = ls())
args<-commandArgs(TRUE)
input_name<-paste(as.character(args[3]),"_gene-characteristic.tmp",sep="")
gene_chara<-read.table(input_name,sep = "\t",header = FALSE)
distance<-dist(gene_chara[,5],method="euclidean")
mode<-hclust(distance,method="ward.D")
f<-function(x){
result<-cutree(mode,k=x)
gene_chara_temp<-cbind(gene_chara,result)
esti_matrix<-matrix(0,nrow=x,ncol=9)
for (i in 1:nrow(gene_chara_temp)){
  if (!is.na(gene_chara_temp[i,5])){
    esti_matrix[gene_chara_temp[i,17],1]<-esti_matrix[gene_chara_temp[i,17],1]+gene_chara_temp[i,5]
    esti_matrix[gene_chara_temp[i,17],6]<-esti_matrix[gene_chara_temp[i,17],6]+1
  }
  if (!is.na(gene_chara_temp[i,6])){
    esti_matrix[gene_chara_temp[i,17],2]<-esti_matrix[gene_chara_temp[i,17],2]+gene_chara_temp[i,6]
    esti_matrix[gene_chara_temp[i,17],7]<-esti_matrix[gene_chara_temp[i,17],7]+1
  }
  if (!is.na(gene_chara_temp[i,8])){
    esti_matrix[gene_chara_temp[i,17],4]<-esti_matrix[gene_chara_temp[i,17],4]+gene_chara_temp[i,8]
    esti_matrix[gene_chara_temp[i,17],8]<-esti_matrix[gene_chara_temp[i,17],8]+1
  }
  if(!is.na(gene_chara_temp[i,9])){
    esti_matrix[gene_chara_temp[i,17],5]<-esti_matrix[gene_chara_temp[i,17],5]+gene_chara_temp[i,9]
    esti_matrix[gene_chara_temp[i,17],9]<-esti_matrix[gene_chara_temp[i,17],9]+1
  }
}
for (i in 1:nrow(esti_matrix)){
  esti_matrix[i,1]<-esti_matrix[i,1]/esti_matrix[i,6]
  esti_matrix[i,2]<-esti_matrix[i,2]/esti_matrix[i,7]
  esti_matrix[i,4]<-esti_matrix[i,4]/esti_matrix[i,8]
  esti_matrix[i,5]<-esti_matrix[i,5]/esti_matrix[i,9]
}
na_pos<-which(is.na(gene_chara),arr.ind=TRUE)
for (i in 1:nrow(na_pos)){
  gene_chara[na_pos[i,1],na_pos[i,2]]<-esti_matrix[gene_chara_temp[na_pos[i,1],17],na_pos[i,2]-4]
}
return(sum(is.na(gene_chara)))
}

left<-1
right<-as.numeric(args[1])
if(left==right){
  cluster_number<-left
  middle<-0
}
while(right-left>1){
middle<-round((left+right)/2)
  if(f(middle)==0){
    left<-middle
  }
  else{
    right<-middle
  }
}
if(left==middle && f(right)==0){
  cluster_number<-right
} else if(left==middle && f(right)!=0){
  cluster_number<-left
} else {
  cluster_number<-left
}
cat("interpolation number is",cluster_number,"\n")
#---------------------------------------------------------------------------------------#
result<-cutree(mode,k=cluster_number)
gene_chara_temp<-cbind(gene_chara,result)
esti_matrix<-matrix(0,nrow=cluster_number,ncol=9)
for (i in 1:nrow(gene_chara_temp)){
  if (!is.na(gene_chara_temp[i,5])){
    esti_matrix[gene_chara_temp[i,17],1]<-esti_matrix[gene_chara_temp[i,17],1]+gene_chara_temp[i,5]
    esti_matrix[gene_chara_temp[i,17],6]<-esti_matrix[gene_chara_temp[i,17],6]+1
  }
  if (!is.na(gene_chara_temp[i,6])){
    esti_matrix[gene_chara_temp[i,17],2]<-esti_matrix[gene_chara_temp[i,17],2]+gene_chara_temp[i,6]
    esti_matrix[gene_chara_temp[i,17],7]<-esti_matrix[gene_chara_temp[i,17],7]+1
  }
  if (!is.na(gene_chara_temp[i,8])){
    esti_matrix[gene_chara_temp[i,17],4]<-esti_matrix[gene_chara_temp[i,17],4]+gene_chara_temp[i,8]
    esti_matrix[gene_chara_temp[i,17],8]<-esti_matrix[gene_chara_temp[i,17],8]+1
  }
  if(!is.na(gene_chara_temp[i,9])){
    esti_matrix[gene_chara_temp[i,17],5]<-esti_matrix[gene_chara_temp[i,17],5]+gene_chara_temp[i,9]
    esti_matrix[gene_chara_temp[i,17],9]<-esti_matrix[gene_chara_temp[i,17],9]+1
  }
}
for (i in 1:nrow(esti_matrix)){
  esti_matrix[i,1]<-esti_matrix[i,1]/esti_matrix[i,6]
  esti_matrix[i,2]<-esti_matrix[i,2]/esti_matrix[i,7]
  esti_matrix[i,4]<-esti_matrix[i,4]/esti_matrix[i,8]
  esti_matrix[i,5]<-esti_matrix[i,5]/esti_matrix[i,9]
}
na_pos<-which(is.na(gene_chara),arr.ind=TRUE)
for (i in 1:nrow(na_pos)){
  gene_chara[na_pos[i,1],na_pos[i,2]]<-esti_matrix[gene_chara_temp[na_pos[i,1],17],na_pos[i,2]-4]
}
#----------------------------------------------------------------------------------------#
#-------------------------cluster number used to estimate BMR----------------------------#
#----------------------------------------------------------------------------------------#
args2<-as.numeric(args[2])
cluster_number4<-args2
distance4<-dist(gene_chara[,c(5,6,8,9)],method="euclidean")
mode4<-hclust(distance4,method="ward.D")
result4<-cutree(mode4,k=cluster_number4)
gene_class_tmp<-cbind(gene_chara,result4)
gene_class<-gene_class_tmp[,c(1,17)]

out_name<-paste(as.character(args[3]),"_gene-class.tmp",sep="")


sink(out_name)
gene_class
