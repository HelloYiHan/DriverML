#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
infileExp <- args[1]
infileObs <- args[2]
if(is.na(infileExp) | is.na(infileObs))stop('No input files in R!')
geneRatio <- 1
indelRatio <- 0.05
ifHeaderExp <- F
ifHeaderObs <- F
output=paste(unlist(strsplit(infileExp,split = "[.]"))[1],'-eta.tmp',sep='')
#*************************************************************************************#
geneTableExp <- read.table(infileExp, header = ifHeaderExp, row.names = 1)
geneTableObs <- read.table(infileObs, header = ifHeaderObs, row.names = 1)
#****************************************************************#
funLogLikelihood <- function(n, N, eta, alpha)
{
  logL <- sum(-N * (eta + alpha) + log((N * (eta + alpha))^n) - log(factorial(n)))
  return(logL)
}
#----------------------------------------
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
#----------------------------------------------

if(nrow(geneTableExp) != nrow(geneTableObs) | any(rownames(geneTableExp) != rownames(geneTableObs)))stop('Invalid row number or row names of the 2 input files!')
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
  dimnames(arrGeneTableExp) <- list(geneName = rownames(geneTableExp), type = c(colnames(geneTableExp)[-ncol(geneTableExp)], 'Fs_indel', 'nFs_indel'))
}else{
  arrGeneTableExp <- as.matrix(geneTableExp)
}

arrGeneTableObs <- array(as.matrix(geneTableObs), dim = c(nGene, nType, nPeople))
dimnames(arrGeneTableObs) <- list(geneName = rownames(geneTableExp), type = colnames(arrGeneTableExp), peopleIdx = paste('p', 1:nPeople, sep = ''))

for(i in 1:nPeople){
  if(any(temp <- arrGeneTableObs[, , i] > arrGeneTableExp))arrGeneTableObs[, , i][temp] <- arrGeneTableExp[temp]
}

M <- arrGeneTableExp[, 1:9]
N <- arrGeneTableExp[, -(1:9)]
obsData <- list(m = arrGeneTableObs[, 1:9, ], n = arrGeneTableObs[, -(1:9), ])

etaEstimate <- apply(obsData$m, 2:3, sum) / apply(M, 2, sum)
meanEtaPrior <- apply(etaEstimate, 1, mean)
write.table(meanEtaPrior,file=output,row.names=F,col.names=F,quote=F,sep='\t',append=F)
varEtaPrior <- apply(etaEstimate, 1, var)
alphaPrior <- (meanEtaPrior * (1 - meanEtaPrior) / varEtaPrior - 1) * meanEtaPrior
betaPrior <- (meanEtaPrior * (1 - meanEtaPrior) / varEtaPrior - 1) * (1 - meanEtaPrior)
etaEstimate <- (apply(obsData$m, 2:3, sum) + alphaPrior) / (apply(M, 2, sum) + alphaPrior + betaPrior) * geneRatio
etaEstimate[is.na(etaEstimate)] <- min(etaEstimate[!is.na(etaEstimate)])

if(ifIndel){
  etaEstimateIndel <- apply(obsData$m, 3, sum) / sum(M)
  meanEtaPriorIndel <- mean(etaEstimateIndel)
  write.table(meanEtaPriorIndel,file=output,row.names=F,col.names=F,quote=F,sep='\t',append=T)
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
