library("gsmr")
data("gsmr")
head(gsmr_data)

gsmr_data <- read.csv(paste(gwasinput,"4079",".inputgsmr",sep=""),headert=T)

freq <- fread("ref/mer_hm3_unr_exc_missnp.freq.frq")
freq <- freq[,c(2,5)]
freq <- data.frame(freq)
gsmr_data <- left_join(gsmr_data,freq,by="SNP")
names(gsmr_data)
gsmr_data <- gsmr_data [,c(1:3,13,5:12)]
names(gsmr_data)[4] <- c("a1_freq")

gsmr_data <- gsmr_data[gsmr_data$bzy_pval > 0.00000005]

write.table(gsmr_data[,c(1,2)], "gsmr_example_snps.allele", col.names=F, row.names=F, quote=F)
# Extract the genotype data from a GWAS dataset using GCTA
gcta64 --bfile gsmr_example --extract gsmr_example_snps.allele --update-ref-allele gsmr_example_snps.allele --recode --out gsmr_example

snp_coeff_id = scan("gsmr_example.xmat.gz", what="", nlines=1)
snp_coeff = read.table("gsmr_example.xmat.gz", header=F, skip=2)

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
cat("The estimated effect of the exposure on outcome: ",gsmr_results$bxy)
cat("Standard error of bxy: ",gsmr_results$bxy_se)
cat("P-value for bxy: ", gsmr_results$bxy_pval)
cat("Number of pleiotropic outliers: ", length(gsmr_results$pleio_snps))
