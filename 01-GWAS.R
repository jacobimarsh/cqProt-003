#conda create -n r-env r-base r-essentials 
#conda activate r-env 
#conda install -c conda-forge blast=2.5 mkl\
#conda install r-rMVP  

library("rMVP")

MVP.Data(fileVCF="myVCF.vcf",
         filePhe="gwas_popfile.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp"
         )


genotype <- attach.big.matrix("mvp.geno.desc")
phenotype <- read.table("mvp.phe",head=TRUE)
map <- read.table("mvp.geno.map" , head = TRUE)



imMVP <- MVP(
    phe=phenotype,
    geno=genotype,
    map=map,
    #K=Kinship,
    #CV.GLM=Covariates,     ##if you have additional covariates, please keep there open.
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,      ##if you have added PC into covariates, please keep there closed.
    nPC.MLM=3,
    nPC.FarmCPU=3,
    priority="speed",       ##for Kinship construction
    #ncpus=10,
    vc.method="BRENT",      ##only works for MLM
    maxLoop=10,
    method.bin="static",      ## "FaST-LMM", "static" (#only works for FarmCPU)
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU")
)
