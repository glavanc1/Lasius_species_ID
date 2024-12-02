# this script takes four arguments:
# $1 is a vcf file

# $2 is the missing data cutoff, that is the proportion of
# loci with missing genotypes above which an individual is discarded

# $3 is the median ADR threshold value, which defines the percentile of the expected ADR distribution
# that is used as cutoff above which individuals are discarded.
# it takes values between 0 and 1 (default = 0.01, meaning the top 1%)

# $4 is the ADR threshold for genotype correction, which defines the percentile of the expected ADR distribution
# that is used as cutoff above which genotypes are corrected (by calling the genotype homozygous for the read supported by the most reads).
# it takes values between 0 and 1 (default = 0.05, meaning the top 5%)


### Part 1: load and format data

# get arguments from the script
# "trailingonly = T" makes sure only user-specified arguments are used
args <- commandArgs(trailingOnly = T)

# Install or load required packages
# packages <- c("vcfR")
# install.packages(setdiff(packages, rownames(installed.packages())), repos = "http://cran.us.r-project.org")
library(vcfR)

# Function to apply a sum like pmin or pmax
# Written by Sven Hohenstein here:
# https://stackoverflow.com/questions/13123638/there-is-pmin-and-pmax-each-taking-na-rm-why-no-psum
psum <- function(...,na.rm=FALSE) {
    dat <- do.call(cbind,list(...))
    res <- rowSums(dat, na.rm=na.rm)
    idx_na <- !rowSums(!is.na(dat))
    res[idx_na] <- NA
    res
}

# get the name of the vcf from the first given argument
vcf_name <- args[1]

# get the threshold above which to filter individuals with too much contamination
# defined as 1 - third given argument because we want the 1 - xth percentile
# 0.99 by default
threshold.discard = 0.99
threshold.discard = 1 - as.numeric(args[3])

# get the threshold above which to correct putative contaminated genotypes
# defined as 1 - fourth given argument because we want the 1 - xth percentile
# 0.99 by default
threshold.correct = 0.99
threshold.correct = 1 - as.numeric(args[4])


# remove ".vcf" or ".recode.vcf" at the end of the vcf file name
# if there are other "." in the name, only what's in front will be retained
name = unlist(strsplit(x = vcf_name, split = "[.]"))[1]
print(paste("name =", name))

# load the vcf file in R
vcf_all <- read.vcfR(vcf_name)

### Part 2: Discard individuals with too much missing data.

# extract genotypes from the vcf
gt_all <- extract.gt(vcf_all, "GT")

# this function sums the number of missing genotypes for a vector
count.na<-function(x) sum(is.na(x))

# counts the number of missing genotypes
# for each column of the gt_all object (i.e. each individual)
# and divide by the length of these columns
prop_missing <- apply(gt_all, 2, count.na)/dim(gt_all)[[1]]

# get the positions of individuals with fewer missing genotypes than the cutoff
retained_individuals <- which(prop_missing < args[2])

png(paste("missing_data_distribution_", name, ".png", sep=""), width = 600, height = 600)
hist(prop_missing*100, breaks=250, las=1, ylab=NA, main=NA, xlab="Percentage of missing data per individual", col="lightgrey")
abline(v=as.numeric(args[2])*100, lty=2, col="red")
dev.off()

print(paste("Retained", length(retained_individuals), "out of", dim(gt_all)[[2]], "individuals (passing the cutoff of", args[2], "% missing data)."))

vcf <- vcf_all[, c(1, retained_individuals + 1)]

### Part 3: Compute ADR

# extract genotypes from the new vcf
gt <- extract.gt(vcf, "GT")

# extract depth values from the vcf
dp <- extract.gt(vcf, "DP", as.numeric=T)

# get allelic depth for the reference allele at heterozygous sites
ad.ref <<- extract.gt(vcf, "AD", as.numeric = T)
ad.ref[which(dp < 8)] <- NA
ad.ref.het <<- ad.ref
ad.ref.het[which(gt %in% c("0/0", "1/1"))] <- NA

# and the alternative allele
ad.alt <<- dp - ad.ref
ad.alt[which(dp < 8)] <- NA
ad.alt.het <<- ad.alt
ad.alt.het[which(gt %in% c("0/0", "1/1"))] <- NA

# compute allelic depth ratio (ADR), defined as
# depth at the allele supported by most reads divided by total depth at the locus.
# takes two matrices with allelic depths (integer) as input
get.ad.ratio <<- function(ad.ref, ad.alt) return(pmax(ad.ref, ad.alt) / mapply(psum, ad.ref, ad.alt, na.rm=T))

# use the function on the whole dataset
AD.ratio <- get.ad.ratio(as.matrix(ad.ref.het), as.matrix(ad.alt.het))


### Part 3: estimate an ADR cutoff above which to correct genotypes

# get positions of heterozygote genotypes
het_genotypes <- which(gt %in% c("0/1", "1/0"))

# subset depth values at heterozygote genotypes only
dp.het <- dp[het_genotypes]

# this function samples a number of reads for one of the alleles
# following a binomial distribution and equal probabilities for both alleles
# given the depth at the locus
sim.ad <- function(depth) {
  rbinom(n = 1, size = depth, prob = 0.5)
}

# visualize overall observed distribution of ADR at heterozygous loci
# and expected distribution based on depth at these loci
png(paste("ADR_observed_expected", name, ".png", sep=""), width = 600, height = 600)

# Plot a histogram of observed ADR at heterozygous genotypes
hist(AD.ratio, freq=F, las=1, main="", xlab="")

# probably not the most efficient way but it works.
# applies sim.ad.ref on observed depth values at heterozygous loci
sim.ad.ref <- unlist(lapply(X = dp.het, FUN = sim.ad))

# create a dataframe with simulated AD at both alleles
ad.df <- data.frame(A=sim.ad.ref, B=dp.het - sim.ad.ref)

# compute ADR for the simulated AD values
sim.ad.ratio <- get.ad.ratio(apply(ad.df, 1, max), apply(ad.df, 1, min))

# plot it as a density
points(density(sim.ad.ratio, bw = 0.01) , type="l", col="#1D71B8")

# define the cutoff for further filtering
# as the "thresholdth" percentile of the simulated distribution
cutoff.discard <<- sort(sim.ad.ratio)[threshold.discard*length(sim.ad.ratio)]
cutoff.correct <<- sort(sim.ad.ratio)[threshold.correct*length(sim.ad.ratio)]
abline(v=cutoff.discard, col="red", lw=2)
abline(v=cutoff.correct, col="red", lw=, lty=2)
dev.off()

### Part 4: defining contaminated individuals

# We define individuals with the majority of ADR values
# above the cutoff as contaminated (and therefore discarded)

# define the median ADR for each individual
ADR_median <- apply(AD.ratio, 2, median, na.rm=T)

# define discarded individuals as those with
# median ADR above cutoff (i.e. majority of ADR values above cutoff).
retained_ADR <- which(ADR_median < cutoff.discard)
discarded_ADR <- which(ADR_median >= cutoff.discard)

col.discarded <- rep("#00000030", length(ADR_median))
col.discarded[discarded_ADR] <- "#A31B5C54"

print(paste("discarded", length(discarded_ADR), "individuals suspected of contamination"))
print(dimnames(gt)[[2]][discarded_ADR])

### Part 5: visualizing the distribution of ADR at heterozygous loci per individual
pdf(paste("ADR_distribution_per_sample_", name, "_discard", args[3], "_correct", args[4], ".pdf", sep=""), width = 40, height = 40)
par(mfrow=c(1,1), mai=c(1, 0.2, 0.2, 0.2))
plot(NULL, xlim=c(0.45, 1), ylim=c(0,15), xlab="Allelic depth ratio (heterozygous alleles only)", ylab="", yaxt="n", las=1)
for (i in 1:dim(AD.ratio)[[2]]) {
  if (sum(is.na(AD.ratio[,i]) == F) > 5) {
    points(density(AD.ratio[,i][which(AD.ratio[,i] != 1)], na.rm=T) , type="l", col=col.discarded[i])
  }
}
# Add the simulated distribution
points(density(sim.ad.ratio, bw = 0.01) , type="l", lwd=3, col="#00A19A")

# Add lines with the two cutoffs
abline(v=cutoff.discard, col="#A31B5C", lw=2)
abline(v=cutoff.correct, col="#A31B5C", lw=2, lty=2)

dev.off()

print("ADR distribution per sample done.")

### Part 6: filtering the vcf

# copying gt for retained individuals
gt.filtered <- gt[, retained_ADR]

# copying dp for retained individuals
dp.filtered <- dp[, retained_ADR]
dp.filtered[which(is.na(gt.filtered))] <- NA

# copying ad.ref for retained individuals
ad.ref.filtered <- ad.ref[, retained_ADR]
ad.ref.filtered[which(is.na(gt.filtered))] <- NA

# creating ad.alt to be filtered
ad.alt.filtered <- dp[, retained_ADR] - ad.ref[, retained_ADR]
ad.alt.filtered[which(is.na(gt.filtered))] <- NA

# identifying cells of the genotype matrix with ADR > cutoff
# and assigning them to either "conta.alt" or "conta.ref"
# in "conta.alt", the alternative allele has the lower coverage and is assumed to be the result of a contamination
# in "conta.ref", the reference allele has the lower coverage and is assumed to be the result of a contamination

# First subsetting AD.ratio to remove contaminated samples
AD.ratio <- AD.ratio [, retained_ADR]

conta.alt <- which(AD.ratio > cutoff.correct & ad.ref.filtered > ad.alt.filtered)
conta.ref <- which(AD.ratio > cutoff.correct & ad.alt.filtered > ad.ref.filtered)

# correcting genotype values
gt.filtered[conta.alt] <- "0/0"
gt.filtered[conta.ref] <- "1/1"

# correcting depth values
dp.filtered[conta.alt] <- ad.ref.filtered[conta.alt]
dp.filtered[conta.ref] <- ad.alt.filtered[conta.ref]

# correcting allelic depth values
ad.alt.filtered[conta.alt] <- 0
ad.ref.filtered[conta.ref] <- 0

# copying the original vcf, retaining only individuals passing ADR cutoff
vcf.filtered <- vcf[, c(0, retained_ADR) + 1]

# modifying the "FORMAT" field
vcf.filtered@gt[,1] <- "GT:DP:AD"

# creating the other columns of the vcf@gt by pasting GT, DP and AD
# and replacing into vcf.filtered to keep all metadata
vcf.filtered@gt[,2:dim(vcf.filtered@gt)[2]] <- paste(gt.filtered, dp.filtered, paste(ad.ref.filtered, ad.alt.filtered, sep=","), sep=":")
vcf.filtered@gt[,2:dim(vcf.filtered@gt)[2]][which(vcf.filtered@gt[,2:dim(vcf.filtered@gt)[2]] == "NA:NA:NA,NA")] <- "./.:.:.,."

# Writing into a file
write.vcf(vcf.filtered, file = paste(name, "_ADR_filtered_discard", args[3], "_correct", args[4], ".vcf.gz", sep=""))

print(paste("Done. writing file ", name, "_ADR_filtered_discard", args[3], "_correct", args[4], ".vcf.gz", sep=""))
