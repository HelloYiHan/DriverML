#!/usr/bin/env Rscript

rm(list = ls())
gc()
args<-commandArgs(TRUE)
#***************************parameter setting*************************#
sample_num<-as.numeric(args[1])
core_number<-as.numeric(args[2])
monte_t<-as.numeric(args[3])
geneRatio<-as.numeric(args[4])
indelRatio<-as.numeric(args[5])
prior<-as.character(args[6])
date<-as.character(args[7])
eps_test<-as.numeric(args[8])
if(is.na(core_number))core_number<-10
if(is.na(monte_t))monte_t<-10000
if(is.na(geneRatio))geneRatio<-1
if(is.na(indelRatio))indelRatio<-0.05
#*********************************************************************#
fun_LRTnew_num <- function(n,N,eta)
{
  logL <- sum(-N + (n / eta))
  return(logL)
}
fun_LRTnew_den <- function(eta)
{
  logL <- sum(1 / eta)
  return(logL)
}
#######################################################################
#pre date_pre=pre_${date}
date<-paste("pre",date,sep="_")
sample_num<-1
date_pre <- paste(date,"_",sep = "")
ex_file <- paste("subclass-",sample_num,sep = "")
exp_file <- paste(ex_file,'.tmp_exp.tmp',sep = "")
exp_file_date<-paste(date_pre,exp_file,sep="")
ob_file <- paste("subclass-",sample_num,sep = "")
obs_file <- paste(ob_file,'.tmp_obs.tmp',sep = "")
obs_file_date<-paste(date_pre,obs_file,sep="")
geneTableExp <- read.table(exp_file_date, header = FALSE, row.names = 1)
geneTableObs <- read.table(obs_file_date, header = FALSE, row.names = 1)
if(nrow(geneTableExp) != nrow(geneTableObs) | any(rownames(geneTableExp) != rownames(geneTableObs))){
  stop('Invalid row number or row names of the 2 input files!')
}
nGene <- nrow(geneTableExp)
nType <- ncol(geneTableExp)
if(nType == 36){
  ifIndel <- F
}else if(nType == 37){
  nType <- 38
  ifIndel <- T
}else stop('Invalid column number of .exp file!')
if(ncol(geneTableObs) %% nType == 0){
  nPeople <- ncol(geneTableObs) / nType
}else stop('Invalid column number of .obs file!')

if(ifIndel){
  arrGeneTableExp <- cbind(as.matrix(geneTableExp), geneTableExp[, ncol(geneTableExp)])
  dimnames(arrGeneTableExp) <- list(geneName = rownames(geneTableExp), type = c(colnames(geneTableExp)[-ncol(geneTableExp)], 'Fs_indel', 'nFs_ind
el'))
}else{
  arrGeneTableExp <- as.matrix(geneTableExp)
}
arrGeneTableObs <- array(as.matrix(geneTableObs), dim = c(nGene, nType, nPeople))
rm(geneTableObs)
gc()
dimnames(arrGeneTableObs) <- list(geneName = rownames(geneTableExp), type = colnames(arrGeneTableExp), peopleIdx = paste('p', 1:nPeople, sep = ''
))
for(i in 1:nPeople){
  if(any(temp <- arrGeneTableObs[, , i] > arrGeneTableExp))arrGeneTableObs[, , i][temp] <- arrGeneTableExp[temp]
}
M <- arrGeneTableExp[, 1:9]
N <- arrGeneTableExp[, -(1:9)]
obsData <- list(m = arrGeneTableObs[, 1:9, ], n = arrGeneTableObs[, -(1:9), ])
rm(arrGeneTableObs)
gc()
etaEstimate <- apply(obsData$m, 2:3, sum) / apply(M, 2, sum)
meanEtaPrior <- apply(etaEstimate, 1, mean)
varEtaPrior <- apply(etaEstimate, 1, var)
alphaPrior <- (meanEtaPrior * (1 - meanEtaPrior) / varEtaPrior - 1) * meanEtaPrior
betaPrior <- (meanEtaPrior * (1 - meanEtaPrior) / varEtaPrior - 1) * (1 - meanEtaPrior)
etaEstimate <- (apply(obsData$m, 2:3, sum) + alphaPrior) / (apply(M, 2, sum) + alphaPrior + betaPrior) * geneRatio
rm(meanEtaPrior,varEtaPrior,alphaPrior,betaPrior)
gc()
if(sum(is.na(etaEstimate))>0){
etaEstimate[is.na(etaEstimate)] <- min(etaEstimate[!is.na(etaEstimate)])
}
if(ifIndel){
  etaEstimateIndel <- apply(obsData$m, 3, sum) / sum(M)
  meanEtaPriorIndel <- mean(etaEstimateIndel)
  varEtaPriorIndel <- var(etaEstimateIndel)
  alphaPriorIndel <- (meanEtaPriorIndel * (1 - meanEtaPriorIndel ) / varEtaPriorIndel - 1) * meanEtaPriorIndel
  betaPriorIndel <- (meanEtaPriorIndel * (1 - meanEtaPriorIndel) / varEtaPriorIndel - 1) * (1 - meanEtaPriorIndel)
  etaEstimateIndel <- (apply(obsData$m, 3, sum) + alphaPriorIndel) / (sum(M) + alphaPriorIndel + betaPriorIndel)
  etaEstimateIndel <- etaEstimateIndel * indelRatio * geneRatio
  etaEstimate <- rbind(etaEstimate, etaEstimate, etaEstimate, etaEstimateIndel, etaEstimateIndel)
  rownames(etaEstimate) <- colnames(N)
}else{
  etaEstimate <- rbind(etaEstimate, etaEstimate, etaEstimate)
  rownames(etaEstimate) <- colnames(N)
}
rm(etaEstimateIndel,meanEtaPriorIndel,varEtaPriorIndel,alphaPriorIndel,betaPriorIndel)
gc()
nTypeTest <- nType - 9
gene_name <- row.names(geneTableExp)
gene_name <- as.character(gene_name)
#***************************************#
census <- read.table(prior,header=F,sep='\t')
census <- census[,1]
#***************************************#
census <- as.character(census)
label<-rep(0,time=length(gene_name))
for (i in 1:length(gene_name)) {
  for (j in 1:length(census)) {
    if (gene_name[i]==census[j]) {
      label[i]<--1
    }
  }
}
pre_gene_name<-gene_name
rm(census,gene_name)
gc()
#******************************************************************#
a<-rep(1,time=nTypeTest)
para_tmp<-array(NA,c(nGene,nTypeTest))
para_tmp2<-array(NA,c(1,nTypeTest))
para_all<-array(NA,c(nGene,nTypeTest))
for (k in which(label!=0)){
  f<-function(w){
    s<-c(0,0,0)
    for (j in 1:nTypeTest){
      s[1]<-s[1]+w[j]*fun_LRTnew_num(obsData$n[k, j, ],N[k,j],etaEstimate[j, ])
      s[2]<-s[2]+w[j]^2*fun_LRTnew_den(etaEstimate[j, ])
    }
    s[3]<-s[3]+s[1]/sqrt(s[2])
    return(-s[3])
  }
  para<-optim(a,f)
  para_tmp[k,]<-para$par
  para_tmp2<-rbind(para_tmp2,para$par)
}
para_tmp2<-para_tmp2[-1,]
para_mean<-colMeans(para_tmp2)
for(k in 1:nGene){
  if(k %in% which(label!=0)){
    para_all[k,]<-para_tmp[k,]
  }
  else{
    para_all[k,]<-para_mean
  }
}
pre_para_all<-para_all
row.names(pre_para_all)<-pre_gene_name
pre_para_mean<-para_mean
print(pre_para_mean)
rm(para,para_tmp,para_tmp2,para_mean,label,para_all)
gc()

######################################################################
date<-as.character(args[7])
sample_num<-as.numeric(args[1])
date_pre <- paste(date,"_",sep = "")
ex_file <- paste("subclass-",sample_num,sep = "")
exp_file <- paste(ex_file,'.tmp_exp.tmp',sep = "")
exp_file_date<-paste(date_pre,exp_file,sep="")
ob_file <- paste("subclass-",sample_num,sep = "")
obs_file <- paste(ob_file,'.tmp_obs.tmp',sep = "")
obs_file_date<-paste(date_pre,obs_file,sep="")
geneTableExp <- read.table(exp_file_date, header = FALSE, row.names = 1)
geneTableObs <- read.table(obs_file_date, header = FALSE, row.names = 1)
if(nrow(geneTableExp) != nrow(geneTableObs) | any(rownames(geneTableExp) != rownames(geneTableObs))){
  stop('Invalid row number or row names of the 2 input files!')
}
nGene <- nrow(geneTableExp)
nType <- ncol(geneTableExp)
if(nType == 36){
  ifIndel <- F
}else if(nType == 37){ 
  nType <- 38
  ifIndel <- T
}else stop('Invalid column number of .exp file!')
if(ncol(geneTableObs) %% nType == 0){
  nPeople <- ncol(geneTableObs) / nType
}else stop('Invalid column number of .obs file!')

if(ifIndel){
  arrGeneTableExp <- cbind(as.matrix(geneTableExp), geneTableExp[, ncol(geneTableExp)])
  dimnames(arrGeneTableExp) <- list(geneName = rownames(geneTableExp), type = c(colnames(geneTableExp)[-ncol(geneTableExp)], 'Fs_indel', 'nFs_ind
el'))
}else{
  arrGeneTableExp <- as.matrix(geneTableExp)
}
arrGeneTableObs <- array(as.matrix(geneTableObs), dim = c(nGene, nType, nPeople))
rm(geneTableObs)
gc()
dimnames(arrGeneTableObs) <- list(geneName = rownames(geneTableExp), type = colnames(arrGeneTableExp), peopleIdx = paste('p', 1:nPeople, sep = ''
))
for(i in 1:nPeople){
  if(any(temp <- arrGeneTableObs[, , i] > arrGeneTableExp))arrGeneTableObs[, , i][temp] <- arrGeneTableExp[temp]
}
M <- arrGeneTableExp[, 1:9]
N <- arrGeneTableExp[, -(1:9)]
obsData <- list(m = arrGeneTableObs[, 1:9, ], n = arrGeneTableObs[, -(1:9), ])
rm(arrGeneTableObs)
gc()
etaEstimate <- apply(obsData$m, 2:3, sum) / apply(M, 2, sum)
meanEtaPrior <- apply(etaEstimate, 1, mean)
varEtaPrior <- apply(etaEstimate, 1, var)
alphaPrior <- (meanEtaPrior * (1 - meanEtaPrior) / varEtaPrior - 1) * meanEtaPrior
betaPrior <- (meanEtaPrior * (1 - meanEtaPrior) / varEtaPrior - 1) * (1 - meanEtaPrior)
etaEstimate <- (apply(obsData$m, 2:3, sum) + alphaPrior) / (apply(M, 2, sum) + alphaPrior + betaPrior) * geneRatio
rm(meanEtaPrior,varEtaPrior,alphaPrior,betaPrior)
gc()
if(sum(is.na(etaEstimate))>0){
etaEstimate[is.na(etaEstimate)] <- min(etaEstimate[!is.na(etaEstimate)])
}
if(ifIndel){
  etaEstimateIndel <- apply(obsData$m, 3, sum) / sum(M)
  meanEtaPriorIndel <- mean(etaEstimateIndel)
  varEtaPriorIndel <- var(etaEstimateIndel)
  alphaPriorIndel <- (meanEtaPriorIndel * (1 - meanEtaPriorIndel ) / varEtaPriorIndel - 1) * meanEtaPriorIndel
  betaPriorIndel <- (meanEtaPriorIndel * (1 - meanEtaPriorIndel) / varEtaPriorIndel - 1) * (1 - meanEtaPriorIndel)
  etaEstimateIndel <- (apply(obsData$m, 3, sum) + alphaPriorIndel) / (sum(M) + alphaPriorIndel + betaPriorIndel)
  etaEstimateIndel <- etaEstimateIndel * indelRatio * geneRatio
  etaEstimate <- rbind(etaEstimate, etaEstimate, etaEstimate, etaEstimateIndel, etaEstimateIndel)
  rownames(etaEstimate) <- colnames(N)
}else{
  etaEstimate <- rbind(etaEstimate, etaEstimate, etaEstimate)
  rownames(etaEstimate) <- colnames(N)
}
rm(etaEstimateIndel,meanEtaPriorIndel,varEtaPriorIndel,alphaPriorIndel,betaPriorIndel)
gc()
nTypeTest <- nType - 9
gene_name <- row.names(geneTableExp)
gene_name <- as.character(gene_name)

para_all<-array(NA,c(nGene,nTypeTest))

#pre
for (k in 1:nGene){
	if(gene_name[k]%in%pre_gene_name){
		para_all[k,]<-pre_para_all[gene_name[k],]
	}
	else{
		para_all[k,]<-pre_para_mean
	}
}
##############################################################
#**the code will be used when there is no overlap between a subclass and Prior *********
if (is.na(sum(para_all))){
  for (k in 1:nGene){
    para_all[k,]<-c(rep(1,9),rep(7,9),rep(5,9),7,3)
  }
  para_all[1,]
}

para_file<-paste(date_pre,sample_num,sep="")
para_file_name<-paste(para_file,"para.tmp",sep="_")

write.table(para_all,para_file_name,quote = F,sep = "\t",row.names = F,col.names = F)

#*******************************************************************************************#
LRT <- pValueLRT <- array(NA, c(nGene, nTypeTest))
denominator_help <- array(NA,c(nGene, nTypeTest))
dimnames(denominator_help) <- dimnames(N)
dimnames(LRT) <- dimnames(pValueLRT) <- dimnames(N)
for(k in 1:nGene)
{
  for(j in 1:nTypeTest)
  {
    denominator_help[k,j] <- fun_LRTnew_den(etaEstimate[j, ]) * para_all[k,j] * para_all[k,j]

  }
}
denominator <- rowSums(denominator_help)
denominator_ok <- rep(NA,nGene)
for(k in 1:nGene){
  denominator_ok[k]<-sqrt(denominator[k])
}
rm(denominator,denominator_help)
gc()
for(k in 1:nGene)
{
  for(j in 1:nTypeTest)
  {
    LRT[k, j] <- fun_LRTnew_num(obsData$n[k, j, ],N[k,j],etaEstimate[j, ]) * para_all[k,j] / denominator_ok[k]
  }
}
#----------------------------------------------------------------------------
epsilon<-eps_test*sum(LRT<0)/length(LRT)
print(epsilon)
etaEstimate<-etaEstimate*(1+epsilon)

LRT <- pValueLRT <- array(NA, c(nGene, nTypeTest))
denominator_help <- array(NA,c(nGene, nTypeTest))
dimnames(denominator_help) <- dimnames(N)
dimnames(LRT) <- dimnames(pValueLRT) <- dimnames(N)
for(k in 1:nGene)
{
  for(j in 1:nTypeTest)
  {
    denominator_help[k,j] <- fun_LRTnew_den(etaEstimate[j, ]) * para_all[k,j] * para_all[k,j]

  }
}
denominator <- rowSums(denominator_help)
denominator_ok <- rep(NA,nGene)
for(k in 1:nGene){
  denominator_ok[k]<-sqrt(denominator[k])
}
rm(denominator,denominator_help)
gc()
for(k in 1:nGene)
{
  for(j in 1:nTypeTest)
  {
    LRT[k, j] <- fun_LRTnew_num(obsData$n[k, j, ],N[k,j],etaEstimate[j, ]) * para_all[k,j] / denominator_ok[k]
  }
}
#--------------------------------------------------------------------------
eta_file<-paste(date_pre,sample_num,sep="")
eta_file_name<-paste(eta_file,"eta.tmp",sep="_")

write.table(etaEstimate,eta_file_name,quote = F,sep = "\t",row.names = F,col.names = F)

rm(para_all,denominator_ok,etaEstimate)
gc()
multipleLRT <- rowSums(LRT)

lrt_file<-paste(date_pre,sample_num,sep="")
lrt_file_name<-paste(lrt_file,"lrt.tmp",sep="_")

write.table(multipleLRT,lrt_file_name,quote = F,sep = "\t",row.names = F,col.names = F)



N_file<-paste(date_pre,sample_num,sep="")
N_file_name<-paste(N_file,"N.tmp",sep="_")

write.table(N,N_file_name,quote = F,sep = "\t",row.names = F,col.names = F)

rm(LRT)
gc()

nMutationSilent <- apply(obsData$m, 1, sum)
nMutationMissense <- apply(obsData$n[, 1:9, ], 1, sum)
nMutationNonsense <- apply(obsData$n[, 10:18, ], 1, sum)
nMutationSplicing <- apply(obsData$n[, 19:27, ], 1, sum)
if(ifIndel){
  nMutationFsIndel <- apply(obsData$n[, 28, ], 1, sum)
  nMutationNFsIndel <- apply(obsData$n[, 29, ], 1, sum)
  nMutationTotal <- nMutationSilent + nMutationMissense + nMutationNonsense + nMutationSplicing + nMutationFsIndel + nMutationNFsIndel
  outSummary <- data.frame(gene = rownames(N), total= nMutationTotal, silent = nMutationSilent, 
                           missense = nMutationMissense, nonsense = nMutationNonsense, splicing = nMutationSplicing, Fs_indel = nMutationFsIndel,
 nFs_indel = nMutationNFsIndel, 
                           LRT = multipleLRT)
}else{
  nMutationTotal <- nMutationSilent + nMutationMissense + nMutationNonsense + nMutationSplicing
  outSummary <- data.frame(gene = rownames(N), total= nMutationTotal, silent = nMutationSilent, 
                           missense = nMutationMissense, nonsense = nMutationNonsense, splicing = nMutationSplicing, 
                           LRT = multipleLRT)
}
out_fi <- paste('out_file_',sample_num,sep = "")
out_file <- paste(out_fi,'.tmp',sep="")
out_file_date<-paste(date_pre,out_file,sep='')
write.table(outSummary, out_file_date, row.names = F, quote = F, sep = '\t')
