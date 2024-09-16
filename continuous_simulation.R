
library(Rcpp,lib.loc="/home/slcxding/R_libs")
library(RcppArmadillo,lib.loc="/home/slcxding/R_libs")
library(CompQuadForm,lib.loc="/home/slcxding/R_libs")
library(rareGE,lib.loc="/home/slcxding/R_libs")
library(iSKAT)
library(foreach,lib.loc="/home/slcxding/R_libs")
library(iterators,lib.loc="/home/slcxding/R_libs")
library(doParallel,lib.loc="/home/slcxding/R_libs")
library(gtx,lib.loc="/home/slcxding/R_libs")
library(MASS,lib.loc="/home/slcxding/R_libs")
library(corpcor,lib.loc="/home/slcxding/R_libs")
library(MiSTi,lib.loc="/home/slcxding/R_libs")
library(aGE,lib.loc="/home/slcxding/R_libs")


sourceCpp("/home/slcxding/File/MAGEB.cpp")
source("/home/slcxding/File/ADABFGE.R")
source("/home/slcxding/File/MAGE.R")


################################################################################
###########################        haplotype data generation      #######################
################################################################################

registerDoParallel(15)
set.seed(666)

ns=10000  ## Set the people number 10000
haplotype = read.table("/home/slcxding/File/haplotype.txt",header=F) #choose the 0-1 haplotype file
n.1 = dim(haplotype)[1]
haplotype = as.matrix(haplotype)
sample1=sample(1:n.1,ns,replace=TRUE)
sample2=sample(1:n.1,ns,replace=TRUE)
genotype=haplotype[sample1,]+haplotype[sample2,] #generate the SNP data
nonraredup.common=genotype[,which(colSums(genotype)>2*ns*0.05)] # rare allel frequency >0.05
nonraredup.common=as.matrix(nonraredup.common)
n.2 = dim(nonraredup.common)[2]  ## The dimension is 10000*576, 10000 people, 576SNP

print(n.2)
nonraredup.rare=genotype[,which(colSums(genotype)>2*ns*0.005 & colSums(genotype)<2*ns*0.05)] 
nonraredup.rare=as.matrix(nonraredup.rare)
n.3 = dim(nonraredup.rare)[2]  
print(n.3)

####### Initialization #########

n = 100000 # iteration numbers
ss = 5000 # sample size
sn_common = 10 # Common SNP numbers
sn_rare = 40 # Rare SNP numbers
cn_common = 2 # Common causal SNP numbers
cn_rare = 8 # Rare causal SNP numbers
n_main = 10
n_common.main = 2
n_rare.main = 8
alpha_common_main = 0.09
alpha_rare_main = 0.17

 



pval_GESAT = function(y, x1, x2, e1, G){

    G = as.matrix(G)
    y = as.matrix(y)
    E = as.matrix(e1)
    X = cbind(x1,x2)
    X = as.matrix(X)
    p_gesat = GESAT(G,y,E,X)$pvalue
  
    return(p_gesat)
}


pval_iSKAT = function(y, x1, x2, e1, G){

    G = as.matrix(G)
    y = as.matrix(y)
    E = as.matrix(e1)
    X = cbind(x1,x2)
    X = as.matrix(X)
         p_iSKAT = iSKAT(G,y,E,X)$pvalue

    return(p_iSKAT)
}

pval_rare = function(y, x1, x2, e1, G){

    X_rare = cbind(e1,x1,x2)
    X_rare = as.matrix(X_rare)     
       rare_result = rareGE(y,G,X_rare)
    pval_intfix = rare_result$pINT_FIX
    pval_intran = rare_result$pINT_RAN
    p = cbind(pval_intfix, pval_intran)

    return(p)

}

pval_MiSTi = function(y, x1, x2, e1, G){

n = ncol(G)
data = cbind(y,e1,x1,x2,G)
data = as.data.frame(data)
result = MiSTi(data,p = n,m=1,d=2)$pvalue

return(result)
}


pval_aGE = function(y, x1, x2, e1, G){

cov = cbind(e1,x1,x2)
cov = as.matrix(cov)
result = aGE(y,G,cov,model="gaussian")
p_aGE = result[7:8]

return(p_aGE)
}


pval_ADABF = function(y, x1, x2, e1, G){
  
  cov = cbind(x1,x2)
  p_adabf = ADABFGE(y, G, e1, Y.Type="C", E.Type="D", Cov=cov, Sig=0.00001, Precision.P=1)
  return(p_adabf)
  
}



###### Claculation Function#########

alpha_common = 0.19
alpha_rare = 0.59


method.name = c("p_MAGE_FIX")

cat(c(method.name,"\n"), file="type1.txt",append=TRUE )


#######################################
oper = foreach(i = 1:n, .combine = rbind) %dopar%
  {

    x1 = rnorm(ss, mean=62.4, sd=11.5)											
    x2 = rbinom(ss, prob=0.52, size=1)											
    e1 = rbinom(ss, prob=0.5, size=1)	  
    indicator = sample(1:10000, ss, replace=F)

    indi_common = sample(1:n.2,sn_common,replace = F)
    none_indi_common = sample(1:sn_common,cn_common)
    G_common = rep(10,sn_common*ss)
    G_common = matrix(G_common,ncol=sn_common)
    SNP=NULL
    for (j in 1:sn_common)
    {
      SNP = nonraredup.common[,indi_common[j]]										
      G_common[,j] = SNP[indicator];
    }


    indi_rare = sample(1:n.3,sn_rare,replace = F)
    none_indi_rare = sample(1:sn_rare,cn_rare)
    G_rare = rep(10,sn_rare*ss)
    G_rare = matrix(G_rare,ncol=sn_rare)
    SNP=NULL
    for (j in 1:sn_rare)
    {
      SNP = nonraredup.rare[,indi_rare[j]]										
      G_rare[,j] = SNP[indicator];
    }



    indi.common.main = sample(1:sn_common, n_common.main, replace=F)
    G.common.main = G_common[, indi.common.main]

    indi.rare.main = sample(1:sn_rare, n_rare.main, replace=F)
    G.rare.main = G_rare[, indi.rare.main]




G = cbind(G_common, G_rare)

colnames(G)=rep("gene",dim(G)[2])
for(i in 1:dim(G)[2]){
colnames(G)[i] = paste("gene",i,sep="")
}

														
    epsilon = rnorm(ss, mean=0, sd=1.5)
    SNP_common = G_common[,none_indi_common]
    SNP_rare = G_rare[,none_indi_rare]
    SNP = cbind(SNP_common, SNP_rare)


alpha1 = runif(n_common.main, min=alpha_common_main-0.02, max=alpha_common_main+0.02)
alpha2 = runif(n_rare.main, min=alpha_rare_main-0.02, max=alpha_rare_main+0.02)
beta1 = runif(cn_common, min=alpha_common-0.02, max=alpha_common+0.02)
beta2 = runif(cn_rare, min=alpha_rare-0.02, max=alpha_rare+0.02)

SNP.sign = rep(c(1,-1),(cn_common+cn_rare))

y_common = y_rare = y_common.main = y_rare.main = rep(0,ss)

for (i in 1:cn_common){
y_common = y_common + beta1[i]*SNP[,i]*SNP.sign[i]*e1
}


for (j in (cn_common+1):(cn_common+cn_rare)){
y_rare = y_rare + beta2[j-cn_common]*SNP[,j]*SNP.sign[j]*e1
}


if (n_common.main == 1){
y_common.main = y_common.main + alpha1*G.common.main*SNP.sign[1]
} else {
for (k in 1:n_common.main){
y_common.main = y_common.main + alpha1[k]*G.common.main[,k]*SNP.sign[k]
}
}



if(n_rare.main == 1){
y_rare.main = y_rare.main + alpha2*G.rare.main*SNP.sign[10]
} else{
for (t in 1:n_rare.main){
y_rare.main = y_rare.main + alpha2[t]*G.rare.main[,t]*SNP.sign[t+n_common.main]
}
}



y = 0.05*x1+0.057*x2+0.64*e1+ epsilon + y_common.main + y_rare.main




fix = try(MAGE_FIX.C(y, x1, x2, e1, G),silent=TRUE)
if('try-error' %in% class(fix)){ print("error_fix")
  p_MAGE_FIX = "NA"
}else{
  p_MAGE_FIX = fix
}





tt = cbind(p_MAGE_FIX)
#tt = cbind(p_MAGE_RAN, p_MAGE_FIX, p_MAGE_BUR)
#tt = cbind(p_MAGE_RAN,p_MAGE_FIX,var.1,var.2,var.3,var.4)
#tt = cbind(p_MAGE_RAN, p_MAGE_FIX, p_MAGE_BUR)
cat(c(tt,"\n"), file="type1.txt",append=TRUE )
}








q("no")

















