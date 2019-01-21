#!/usr/bin/env Rscript
rm(list = ls())
args<-commandArgs(TRUE)
number<-as.numeric(args[1])
geneRatio<-as.numeric(args[2])
name<-as.character(args[3])
date<-as.character(args[4])

out_file_1<-paste(date,'_out_file_1.tmp',sep='')
outSummary_all<-read.table(out_file_1,header = TRUE,sep = '\t')

cat("read 1")
p_file_name<-paste(date,'_1_p.tmp',sep='')

pvalues<-read.table(p_file_name,header=F,sep='\t',col.names=c('p_no_N_2s','p_no_N_2s_adj'))

cat("read 2")
out_with_p<-data.frame(outSummary_all,pvalues)
out_with_p<-out_with_p[1,]
if(number==1){
for(f in 1:number){
  out_fi <- paste("_out_file_",f,sep = "")
  out_file <- paste(out_fi,'.tmp',sep="")
  out_file_date<-paste(date,out_file,sep='')
  file<-read.table(out_file_date,header = TRUE,sep = '\t')

  cat(f)
  p_file<-paste(date,f,sep='_')
  p_file_name<-paste(p_file,'_p.tmp',sep='')

  pfile<-read.table(p_file_name,header=F,sep='\t',col.names=c('p_no_N_2s','p_no_N_2s_adj'))


  out_with_p_loop<-data.frame(file,pfile)
  out_with_p<-rbind(out_with_p,out_with_p_loop)
}

out_with_p<-out_with_p[-1,]

write.table(out_with_p[order(out_with_p$p_no_N_2s),],name,row.names=F,quote=F,sep='\t')

}
if(number>1){

for(f in 1:number){
  out_fi <- paste("_out_file_",f,sep = "")
  out_file <- paste(out_fi,'.tmp',sep="")
  out_file_date<-paste(date,out_file,sep='')
  file<-read.table(out_file_date,header = TRUE,sep = '\t')

  cat(f)
  p_file<-paste(date,f,sep='_')
  p_file_name<-paste(p_file,'_p.tmp',sep='')

  pfile<-read.table(p_file_name,header=F,sep='\t',col.names=c('p_no_N_2s','p_no_N_2s_adj'))

  out_with_p_loop<-data.frame(file,pfile)
  out_with_p_every<-rbind(out_with_p,out_with_p_loop)
  out_with_p_every<-out_with_p_every[-1,]
  name<-as.character(args[3])
  name<-paste(f,name,sep="_")
  write.table(out_with_p_every[order(out_with_p_every$p_no_N_2s),],name,row.names=F,quote=F,sep='\t')

}

}
