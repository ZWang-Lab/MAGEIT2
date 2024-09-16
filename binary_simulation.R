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
library(truncnorm,lib.loc="/home/slcxding/R_libs")

sourceCpp("/home/slcxding/File/MAGEB.cpp")
source("/home/slcxding/File/ADABFGE.R")
source("/home/slcxding/File/MAGE.R")

################################################################################
###########################        haplotype data generation      #######################
################################################################################

registerDoParallel(5)
set.seed(666)

ns=80000  ## Set the people number 10000
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

filename = 'type1.txt'
n = 50000 # iteration numbers
ss = 5000 # sample size
sn_common = 10 # Common SNP numbers
sn_rare = 40 # Rare SNP numbers
cn_common = 2 # Common causal SNP numbers
cn_rare = 8 # Rare causal SNP numbers
n_main = 10
n_common.main = 2
n_rare.main = 8
alpha_common_main = 0.1
alpha_rare_main = 0.2
ss.raw = 50000
popu.pre = 0.05

###### Claculation Function#########

alpha_common = 0.3
alpha_rare = 0.88

#method.name = c("p_MAGE_RAN", "p_MAGE_FIX", "p_adabf", "p_GESAT", "p_iSKAT", "p_rarefix", "p_rareran", "p_Mif",  
#"p_Mir",  "p_Mio",  "p_Mia", "p_Miff" )

method.name = c("p_MAGE_RAN", "p_MAGE_FIX")

cat(c(method.name,"\n"), file=filename,append=TRUE )

#######################################
oper = foreach(i = 1:n, .combine = rbind) %dopar%
{

	x1.raw = rnorm(ss.raw, mean=62.4, sd=11.5)											
	x2.raw = rbinom(ss.raw, prob=0.52, size=1)											
	e1.raw = rbinom(ss.raw, prob=0.5, size=1) 
	indicator = sample(1:ns, ss.raw, replace=F)

	indi_common = sample(1:n.2,sn_common,replace = F)
	none_indi_common = sample(1:sn_common,cn_common)
	G_common = rep(10,sn_common*ss.raw)
	G_common = matrix(G_common,ncol=sn_common)
	SNP=NULL
	for (j in 1:sn_common)
	{
		SNP = nonraredup.common[,indi_common[j]]										
		G_common[,j] = SNP[indicator]
	}

	indi_rare = sample(1:n.3,sn_rare,replace = F)
 	none_indi_rare = sample(1:sn_rare,cn_rare)
 	G_rare = rep(10,sn_rare*ss.raw)
	G_rare = matrix(G_rare,ncol=sn_rare)
	SNP=NULL
	for (j in 1:sn_rare)
	{
		SNP = nonraredup.rare[,indi_rare[j]]										
		G_rare[,j] = SNP[indicator]
	}

	G.raw = cbind(G_common, G_rare)

	colnames(G.raw)=rep("gene",dim(G.raw)[2])
	for(i in 1:dim(G.raw)[2]){
		colnames(G.raw)[i] = paste("gene",i,sep="")
	}

	indi.common.main = sample(1:sn_common, n_common.main, replace=F)
	G.common.main = G_common[, indi.common.main]

	indi.rare.main = sample(1:sn_rare, n_rare.main, replace=F)
	G.rare.main = G_rare[, indi.rare.main]
														
	#epsilon = rnorm(ss.raw, mean=0, sd=1.5)
	SNP_common = G_common[,none_indi_common]
	SNP_rare = G_rare[,none_indi_rare]
	SNP = cbind(SNP_common, SNP_rare)

	alpha1 = runif(n_common.main, min=alpha_common_main-0.02, max=alpha_common_main+0.02)
	alpha2 = runif(n_rare.main, min=alpha_rare_main-0.02, max=alpha_rare_main+0.02)
	beta1 = runif(cn_common, min=alpha_common-0.02, max=alpha_common+0.02)
	beta2 = runif(cn_rare, min=alpha_rare-0.02, max=alpha_rare+0.02)

	SNP.sign = rep(c(1,1),(cn_common+cn_rare))

	y_common.main = y_rare.main = rep(0,ss.raw)

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

	linear = -6.2 + 0.05*x1.raw+0.057*x2.raw+0.64*e1.raw + y_common.main + y_rare.main

	pi = exp(linear)/(1+exp(linear))
	y.raw = rbinom(ss.raw, prob=pi, size=1)
	#pre = mean(pi)
	#print(pre)

################# choose case-control samples #############################

	y.case.index = y.raw == 1
      y.control.index = y.raw == 0
      
      y.case = y.raw[y.case.index][1:2500]
      y.control = y.raw[y.control.index][1:2500]
      y = c(y.case, y.control)
          
      x1.case = x1.raw[y.case.index][1:2500]
      x1.control = x1.raw[y.control.index][1:2500]
      x1 = c(x1.case, x1.control)
      
      x2.case = x2.raw[y.case.index][1:2500]
      x2.control = x2.raw[y.control.index][1:2500]
      x2 = c(x2.case, x2.control)
      
      e1.case = e1.raw[y.case.index][1:2500]
      e1.control = e1.raw[y.control.index][1:2500] 
      e1 = c(e1.case, e1.control)
      
      G.case = G.raw[y.case.index,][1:2500,]
      G.control = G.raw[y.control.index,][1:2500,]
      G = rbind(G.case, G.control)


	RAN = try(MAGE_RAN.B(y, x1, x2, e1, G, 2),silent=TRUE)
	if('try-error' %in% class(RAN)){ print("error_RAN")
		p_MAGE_RAN="NA"
	}else{
		p_MAGE_RAN=RAN
	}

	FIX = try(MAGE_FIX.B(y, x1, x2, e1, G),silent=TRUE)
	if('try-error' %in% class(FIX)){ print("error_FIX")
		p_MAGE_FIX="NA"
	}else{
		p_MAGE_FIX=FIX
	}

	#tt = cbind(p_MAGE_RAN, p_MAGE_FIX, p_adabf, p_GESAT, p_iSKAT, p_rarefix, p_rareran, p_Mif,  p_Mir,  p_Mio,  p_Mia, p_Miff )
	#tt = cbind(p_MAGE_RAN, p_MAGE_FIX, p_MAGE_BUR)
	#tt = cbind(p_MAGE_RAN,p_MAGE_FIX,var.1,var.2,var.3,var.4)
	tt = cbind(p_MAGE_RAN, p_MAGE_FIX)
	write.table(tt,filename,quote=F,row.names=F,col.names=F, append=T)
}



q("no")



