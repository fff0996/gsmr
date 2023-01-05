#clumping

G1000_bim <- fread("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_88_MR/33_GSMR/ref/1000G_hm3.bim")
G1000 <- data.frame(G1000_bim)
names(G1000)[2] <- c("SNP")
gwas <- read.csv("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/33_GWAS_result/input_ldsc/X20151",sep="\t")
m <- left_join(G1000,gwas,by='SNP')
m2 <- na.omit(m)
m3 <- m2[,c(7,2,8:14)]
write.table(m3,"/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_88_MR/33_GSMR/T2D/clumpedQTGWAS/X30120.clumpedinput",sep="\t",quote=F,row.names=F)

plink --bfile /BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_88_MR/ref/1000G_hm3 --clump-p1 0.00000005 --clump-p2 1 --clump-r2 0.05 --clump-kb 1000 --clump /BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/33_GWAS_result/input_ldsc/X${a[$i]} --out X${a[$i]}.GWAS.clump

#make inputdata(QT->BT)
qt <- read.table("../clumpedQTGWAS/X20151.GWAS.clump.clumped",header=T)
qt <- qt$SNP
qt <- data.frame(SNP=qt)

gwas <- read.csv("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/33_GWAS_result/input_ldsc/X20151",sep="\t")
m <- left_join(qt,gwas,by="SNP")
m <- na.omit(m)
m <- m[,c(1,4,5,6,7,8,9)]
names(m) <- c("SNP","a1","a2","bzx","bzx_se","bzx_pval","bzx_n")

E11 <- read.table("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_88_MR/11_GWAS_results_inVALIDANCE/E11.GWAS.assoc.logistic",header=T)
E11 <- E11[,c(1,2,3,4,6,7,8,12)]
m2 <- left_join(m,E11,by="SNP")
m3 <- na.omit(m2)
m3 <- m3[,c(1:7,11,12,13,14)]

names(m3) <- c("SNP","a1","a2","bzx","bzx_se","bzx_pval","bzx_n","bzy_n","OR","bzy_se","bzy_pval")
m3["bzy"] <- log(m3$OR)
m3 <- m3[,c(1,2,3,4,5,6,7,12,10,11,8)]
write.table(m3,"X20151.inputgsmr",sep="\t",quote=F,row.names=F)

##          SNP a1 a2    a1_freq     bzx bzx_se  bzx_pval    bzx_n       bzy
## 1 rs10903129  A  G 0.45001947 -0.0328 0.0037 3.030e-17 169920.0  0.008038
## 2 rs12748152  T  C 0.08087758  0.0499 0.0066 3.209e-12 172987.5  0.013671
## 3 rs11206508  A  G 0.14396988  0.0434 0.0055 2.256e-14 172239.0  0.030222
## 4 rs11206510  C  T 0.19128911 -0.0831 0.0050 2.380e-53 172812.0 -0.074519
## 5 rs10788994  T  C 0.18395430  0.0687 0.0049 8.867e-41 172941.9  0.038267
## 6   rs529787  G  C 0.19713099 -0.0553 0.0052 8.746e-24 161969.0  0.001707
##      bzy_se     bzy_pval  bzy_n
## 1 0.0092442 0.3845651000 184305
## 2 0.0185515 0.4611690000 184305
## 3 0.0141781 0.0330400000 184305
## 4 0.0133438 0.0000000234 184305
## 5 0.0118752 0.0012711000 184305
## 6 0.0135491 0.8997431000 184305




library("gsmr")
data("gsmr")


gsmr_data <- read.csv(paste(,"4079",".inputgsmr",sep=""),sep="\t")

freq <- fread("ref/mer_hm3_unr_exc_missnp.freq.frq")
freq <- freq[,c(2,5)]
freq <- data.frame(freq)
gsmr_data <- left_join(gsmr_data,freq,by="SNP")
names(gsmr_data)
gsmr_data <- gsmr_data [,c(1:3,12,4:11)]
names(gsmr_data)[4] <- c("a1_freq")

gsmr_data <- gsmr_data[gsmr_data$bzy_pval > 0.00000005,]

write.table(gsmr_data[,c(1,2)], "gsmr_example_snps.allele", col.names=F, row.names=F, quote=F)
#LDcorrelation
gcta64 --bfile /BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_88_MR/ref/1000G_hm3 
--extract X${a[$i]}.gsmr_example_snps.allele --update-ref-allele gsmr_example_snps.allele --recode --out X${a[$i]}.gsmr_example

# Extract the genotype data from a GWAS dataset using GCTA
gcta64 --bfile gsmr_example --extract gsmr_example_snps.allele --update-ref-allele gsmr_example_snps.allele --recode --out gsmr_example

snp_coeff_id = scan("gsmr_example.xmat.gz", what="", nlines=1)
snp_coeff = read.table("gsmr_example.xmat.gz", header=F, skip=2)
#NA가 있을경우 
mean(snp_coeff$V127,na.rm=TRUE)

snp_id = Reduce(intersect, list(gsmr_data$SNP, snp_coeff_id))
gsmr_data = gsmr_data[match(snp_id, gsmr_data$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]

# Calculate the LD correlation matrix
ldrho = cor(snp_coeff)

# Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
colnames(ldrho) = rownames(ldrho) = snp_coeff_id

snpfreq = gsmr_data$a1_freq             # allele frequencies of the SNPs
bzx = gsmr_data$bzx     # effects of the instruments on risk factor
bzx_se = gsmr_data$bzx_se       # standard errors of bzx
bzx_n = gsmr_data$bzx_n          # GWAS sample size for the risk factor
std_zx = std_effect(snpfreq, bzx, bzx_se, bzx_n)    # perform standardisation
gsmr_data$std_bzx = std_zx$b    # standardized bzx
gsmr_data$std_bzx_se = std_zx$se    # standardized bzx_se
head(gsmr_data)

bzx = gsmr_data$std_bzx    # SNP effects on the risk factor
bzx_se = gsmr_data$std_bzx_se    # standard errors of bzx
bzx_pval = gsmr_data$bzx_pval   # p-values for bzx
bzy = gsmr_data$bzy    # SNP effects on the disease
bzy_se = gsmr_data$bzy_se    # standard errors of bzy
bzy_pval = gsmr_data$bzy_pval   

n_ref = 404 # Sample size of the reference sample
gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snps_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
filtered_index=gsmr_results$used_index

cat("Hypertension -> ",i,"\n",file="output.txt",append=T)
cat("The estimated effect of the exposure on outcome: ",gsmr_results$bxy,"\n",file="output.txt",append=T)
cat("Standard error of bxy: ",gsmr_results$bxy_se,"\n",file="output.txt",append=T)
cat("P-value for bxy: ", gsmr_results$bxy_pval,"\n",file="output.txt",append=T)
cat("Number of pleiotropic outliers: ", length(gsmr_results$pleio_snps),"\n",file="output.txt",append=T)

#cat("The estimated effect of the exposure on outcome: ",gsmr_results$bxy)
#cat("Standard error of bxy: ",gsmr_results$bxy_se)
#cat("P-value for bxy: ", gsmr_results$bxy_pval)
#cat("Number of pleiotropic outliers: ", length(gsmr_results$pleio_snps))
