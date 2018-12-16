#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# v2.0 by Victor Missirian, Elias Scheer and Colette L. Picard
# 11/09/2018
# ------------------------------------------------------------------------------------

# Usage:
# detect_imprinting.R [options] AxB_counts.txt BxA_counts.txt outfile

# Description:
# -------------------------
# Accepts a counts file from cross AxB and a counts file from the reciprocal cross BxA, containing
# allele-specific counts (counts from the A allele and counts from the B allele), with the following 
# format (no header, 3 columns, first = gene ID, counts from A, counts from B):
#		AL16048391	82	814
#		AL481009	440	124
# Counts files for AxB and BxA should contain the same genes in the same order

# Additional options:
# - User can also specify the expected ratio of maternal to paternal reads (default 1, but in
#   the triploid endosperm, which has 2 maternal : 1 paternal genomes, the ratio is 2).
# - Altering --pval_cutoff and --alpha changes the genes for which the IF or CEF is calculated (see below)
#   but is probably unnecessary for the majority of users.

# OVERVIEW:
# (1) Test for imprinting (compute Fisher's exact test p-value, Imprinting Factor)
# (2) Test for cis-effects (compute Fisher's exact test p-value, Cis-effects Factor)
# (3) Output results to indicated output file

# Procedure in more detail:
# -------------------------
# For each gene, performs the following analysis:
# (1) are there enough allelic reads to evaluate imprinting in both AxB and BxA? (see --minallelic option)
# 	- if NO -> STOP (not enough data to evaluate)
#   - ELSE:
# (2) perform Fisher's exact test to evaluate if there is significant parental bias (note that for
# cases where the expected maternal:paternal ratio != 1, counts must be adjusted to simulate the extra
# maternal or paternal chromosomes, since Fisher's can only test whether or not the two crosses have ratio == 1
#
# Let pat_counts = paternal counts, mat_counts = maternal counts, mp_ratio = expected maternal:paternal ratio,
# and tot_counts = mat_counts + pat_counts. Then adjusted counts are computed as:
#
#			adj_pat_counts = round( [ (mp_ratio * pat_counts) / ((mp_ratio * pat_counts) + mat_counts) ] * tot_counts )
#			adj_mat_counts = tot_counts - adj_pat_counts
#
# Counts for AxB and BxA are adjusted separately. The contingency  table for Fisher's exact test is then:
#
#			M = | AxB_A		AxB_B |
#				| BxA_A		BxA_B |
#
# (3) Compute imprinting factor (IF)
# The bias-factor (imprinting factor when evaluating for imprinting, cis-effect factor when evaluating strain
# effects) is effectively the ratio of the lower bound for %A in the direction of the cross with higher %A
# (this will be AxB if gene is maternally biased, BxA if gene is paternally biased) to the higher bound of %
# maternal in the cross with lower %A. Because the expected maternal:paternal ratio may not be 1, this value
# is adjusted for this expected ratio. A brief overview of how this is calculated (using the binconf function):

# CI_AxB <- binconf(AxB_A, AxB_A+AxB_B, alpha=0.05)
# AxB_lower = CI_AxB[,1] 	# lower bound of confidence interval
# AxB_upper = CI_AxB[,2]	# upper bound

# CI_BxA <- binconf(BxA_A, BxA_A+BxA_B, alpha=0.05)
# BxA_lower = CI_BxA[,1] 	# lower bound of confidence interval
# BxA_upper = CI_BxA[,2]	# upper bound

# Next, identify the cross with higher pA = (A counts) / (A+B counts):
# - take the *lower* bound of the CI from the cross with higher pA -> call this p_hi
# - take the *upper* bound of the CI from the cross with lower pA -> call this p_lo
# Convert p_hi and p_lo to r_hi and r_lo, where r = p / (1-p) (ratio of A:B)
# Multiply the ratio from the cross with B as the father by (mp_ratio)^2:
# - if AxB has higher pA than BxA, r_lo was from BxA -> get r_lo* = r_lo * (mp_ratio^2)
# - if BxA has higher pA than AxB, r_hi was from AxB -> get r_hi* = r_hi * (mp_ratio^2)
# Finally, compute bias_factor = r_hi / r_lo

# (4) Imprinting p-values are corrected for multiple hypotheses using the Benjamini-
# Hochberg method; for genes where this adjusted p-value is not significant (>= pval_cutoff)
# the IF is censored

# (5) Compute p-value and bias-factor for cis-effects (Cis-effects factor or CEF)
# This is computed using the same function that computes these values for imprinting,
# by providing the same values to the function (AxB_A, AxB_B, BxA_A, BxA_B) but switching
# the order of BxA_A and BxA_B, so that the contingency table becomes:
#
#			M = | AxB_A		AxB_B |
#				| BxA_B		BxA_A |
#
# Now counts from the same strain, instead of the same parent, are on the diagonals,
# and Fisher's will evaluate whether there is significant bias in favor of a particular
# strain using the same approach. Expected mp_ratio == 1 for all cis-effect tests. Bias
# factor is calculated as above, but with BxA_A and BxA_B flipped, so that CI_BxA is
# the confidence interval around pB = (B counts)/ (A+B counts) instead of pA. Same
# logic applies.


# -------------------------
# Version history:
# v.1.0: initial version - Victor Missirian
# v.1.1: 06/22/2014 - minor changes by Elias Scheer:
#		- change to Fisher's exact test instead of storer-kim method
#		- arguments are passed to script via a separate ARGS file containing the args
# 		- renamed script from detect.imprinting.R to detect.imprinting.FishersOnly.R
# v.1.2: 04/10/2015 - minor changes by CLP:
#		- arguments to the script now provided on the command line instead of in an ARGS file
#		- added some notes to the functions
# 		- renamed script from detect.imprinting.FisherOnly.R to detect_imprinting.R
# v.2.0: 11/09/2018 - major clean up by CLP:
#		- generally condensed and cleaned up code, dropped things we don't use anymore
#		- incorporated bias functions into this script (originally in a separate script called bias.ratio.helper.functions.R)
#		- changed varname syntax to AxB vs BxA instead of RefxOth, OthxRef for consistency with other scripts
# -------------------------

# -------------------------
# Load required R packages
# -------------------------
# First, check that all packages are in fact installed
if (is.element('stats', installed.packages()[,1])==FALSE) { stop("Error: R package stats could not be loaded") }
if (is.element('Hmisc', installed.packages()[,1])==FALSE) { stop("Error: R package Hmisc could not be loaded") }
if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
# Actually load the packages
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(argparse))

# -------------------------
# Functions used in this script
# -------------------------

# simulate_extra_pat - helper function for test_for_imprinting
# -------------------------
# The way we implemented Fisher's exact test (see test_for_imprinting) requires that
# the expected ratio of maternal:paternal be 1, so read counts must be adjusted to 
# "simulate" extra maternal or paternal expression in cases where this ratio is not 1,
# to approximate what would be observed if the ratio was 1. Else if, for example, the
# ratio was 2, the contingency table c(c(2,1),c(1,2)) would appear significant even
# though 2:1 maternal:paternal is expected (and is not imprinting) in this tissue
"simulate_extra_pat" <- function(mat_counts, tot_counts, mp_ratio){
	# in the case of a mat to pat ratio of > 1, we want to increase the proportion of pat alleles
	pat_counts <- tot_counts - mat_counts

	## double number of Y alleles
	pat_frac_adj <- (mp_ratio * pat_counts) / ((mp_ratio * pat_counts) + mat_counts)

	pat_counts_adj <- round(tot_counts * pat_frac_adj)
	mat_counts_adj <- tot_counts - pat_counts_adj

	list(mat_counts_adj = mat_counts_adj, pat_counts_adj = pat_counts_adj)
}

# get_binconf - helper function for get_bias_factor
# -------------------------
# Returns lower and upper bound of Binomial confidence interval
"get_binconf" <- function(x, n, alpha){
	CI <- binconf(x, n, alpha=alpha)
	est_pA <- CI[,1]
	low_pA <- CI[,2]
	high_pA <- CI[,3]

	list(low_pA=low_pA, high_pA=high_pA)
}

# get_bias_factor - helper function for test_for_imprinting
# -------------------------
# Desc - TODO
"get_bias_factor" <- function(AxB_A_counts, AxB_tot_counts, BxA_A_counts, BxA_tot_counts, mp_ratio, favored_direction, alpha){
	
	res <- get_binconf(AxB_A_counts, AxB_tot_counts, alpha)
	AxB_low_pA <- res$low_pA
	AxB_high_pA <- res$high_pA
	
	res <- get_binconf(BxA_A_counts, BxA_tot_counts, alpha)
	BxA_low_pA <- res$low_pA
	BxA_high_pA <- res$high_pA
	
	if (favored_direction == "X") {
		if (BxA_high_pA == 1.0) {
			## since bias ratio is (constant / infinity)
			bias_ratio <- 0
		} else {
			BxA_high_ratioAB <- BxA_high_pA / (1 - BxA_high_pA)
			AxB_low_ratioAB <- AxB_low_pA / (1 - AxB_low_pA)
                        
			## test for higher than expected maternal expression
			bias_ratio <- (AxB_low_ratioAB / ((mp_ratio^2) * BxA_high_ratioAB))
		}
	} else if (favored_direction == "Y") {
		if (AxB_high_pA == 1.0) {
			## since bias ratio is (constant / infinity)
			bias_ratio <- 0
		} else {
			AxB_high_ratioAB <- AxB_high_pA / (1 - AxB_high_pA)
			BxA_low_ratioAB <- BxA_low_pA / (1 - BxA_low_pA)
                        
			## test for higher than expected paternal expression
			bias_ratio <- (((mp_ratio^2) * BxA_low_ratioAB) / AxB_high_ratioAB)
		}
	} else {
		stop("Expected that at this point in the script, favored_direction would be either 'X' or 'Y'.")
	}
	
	if(bias_ratio < 1){
		## set bias_factor to 0, when we cannot call
		## the tested effect with high confidence (determined by alpha)
		bias_factor <- 0
	}else{
		bias_factor <- bias_ratio
	}
	bias_factor
}

# test_for_imprinting - main function for determining imprinting
# -------------------------
## Desc - TODO
"test_for_imprinting" <- function(AxB_A_counts, AxB_B_counts, BxA_A_counts, BxA_B_counts, mp_ratio, alpha) {

	# Determine parent that is expressed more than expected (not necessarily to a significant degree)
	AxB_tot_counts = AxB_A_counts + AxB_B_counts
	BxA_tot_counts = BxA_A_counts + BxA_B_counts
	AxB_pA <- AxB_A_counts / AxB_tot_counts
	BxA_pA <- BxA_A_counts / BxA_tot_counts
							
	if((AxB_pA == 1.0) && (BxA_pA == 1.0)){
		favored_dir <- "none"
	} else if (AxB_pA == 1.0) {
		favored_dir <- "X"
	} else if (BxA_pA == 1.0) {
		favored_dir <- "Y"
	} else {
		AxB_ABratio <- (AxB_pA) / (1 - AxB_pA)
		BxA_ABratio <- (BxA_pA) / (1 - BxA_pA)
		if (AxB_ABratio > (mp_ratio^2) * BxA_ABratio) {
			favored_dir <- "X"
		} else if (AxB_ABratio == (mp_ratio^2) * BxA_ABratio) {
			favored_dir <- "none"
		} else {
			favored_dir <- "Y"
		}
	}
	
	# If expected ratio of maternal to paternal != 1, need to adjust counts for Fisher's exact test
	rr = simulate_extra_pat(AxB_A_counts, AxB_tot_counts, mp_ratio)
	AxB_A_counts_adj = rr$mat_counts_adj
	AxB_B_counts_adj = rr$pat_counts_adj
	
	rr = simulate_extra_pat(BxA_B_counts, BxA_tot_counts, mp_ratio)
	BxA_B_counts_adj = rr$mat_counts_adj
	BxA_A_counts_adj = rr$pat_counts_adj
		
	# Perform Fisher's exact test
	AxB_vals = c(AxB_A_counts_adj, AxB_B_counts_adj)		# AxB maternal, AxB paternal
	BxA_vals = c(BxA_A_counts_adj, BxA_B_counts_adj)		# BxA paternal, BxA maternal
	contingency_mat = rbind(AxB_vals, BxA_vals)
	
	fisher_res <- fisher.test(contingency_mat)
	pvalue <- fisher_res$p.value
	
	# calculate bias factor (imprinting factor if testing imprinting, cis-effect factor if testing
	# for cis-effects)	
	if (favored_dir != "none") {
		bias_factor <- get_bias_factor(AxB_A_counts, AxB_tot_counts, BxA_A_counts, BxA_tot_counts, mp_ratio, favored_dir, alpha)
	} else {
		bias_factor = 0
	}
	
	list(favored_dir = favored_dir, pval = pvalue, bias_factor = bias_factor)	
}


# -------------------------
# MAIN
# -------------------------

# read in command-line arguments
parser = ArgumentParser()
parser$add_argument("AxB_file", nargs=1, help = "tab-delimited text file of raw allelic counts for cross AxB (gene ID, A counts, B counts)")
parser$add_argument("BxA_file", nargs=1, help = "tab-delimited text file of raw allelic counts for cross BxA (gene ID, A counts, B counts)")
parser$add_argument("outfile", nargs=1, help = "name for output file")
parser$add_argument("--mp_ratio", type="double", default = 1, help = "expected ratio of maternal to paternal reads (1 for most tissues, 2 for endosperm)")
parser$add_argument("--pval_cutoff", type="double", default = 0.05, help = "p-value cutoff to consider imprinting or cis-effects")
parser$add_argument("--alpha", type="double", default = 0.05, help = "alpha value used to compute Binomial confidence intervals used in IF/CEF calculations - default 0.05")
parser$add_argument("--minallelic", type="integer", default = 1, help = "Min. number of allelic reads in each direction of cross in order to evaluate imprinting")

args <- commandArgs(trailingOnly = TRUE)

# the -0 option is for checking that all dependencies are installed; exit gracefully here if so
if (length(args) == 1 & args[1] == "-0") {
	q(save="no", status=0, runLast=FALSE)
}

if (length(args) <= 2) {
	cat("Usage: detect_imprinting.R [options] AxB_counts.txt BxA_counts.txt outfile\n")
	cat("----------------------\n")
	cat("Required arguments:\n")
	cat("AxB_counts.txt : tab-delimited text file of raw allelic counts for cross AxB (gene ID, A counts, B counts)\n")
	cat("BxA_counts.txt : tab-delimited text file of raw allelic counts for cross BxA (gene ID, A counts, B counts)\n")
	cat("outfile.txt : name to use for output file\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--mp_ratio <int> : ratio of expected maternal to paternal expression (e.g. 1 for embryo, 2 for endosperm) - default 1\n")
	cat("--pval_cutoff <float> : p value cutoff for imprinted locus (IF is not reported for genes w/ padj > pval_cutoff) - default 0.05\n")
	cat("--alpha <float> : alpha value used to compute Binomial confidence intervals used in IF/CEF calculations - default 0.05\n")
	cat("--minallelic <int> : min number of allelic reads in both AxB and BxA to evaluate imprinting - default 1\n")
	cat("----------------------\n")
	stop("At least one required argument not provided. See usage above.")
}

opt <- parser$parse_args()

AxB_file = opt$AxB_file
BxA_file = opt$BxA_file
outfile = opt$outfile
mp_ratio = opt$mp_ratio
pval_cutoff = opt$pval_cutoff
alpha = opt$alpha
minallelic = opt$minallelic

cat("\nRunning detect_imprinting.R v.2.0 (11/09/2018)\n")
cat("----------------------\n")
cat("RefxAlt infile: ", AxB_file,"\n")
cat("AltxRef infile: ", BxA_file,"\n")
cat("Output file: ", outfile,"\n")
cat("----------------------\n")
cat("Additional parameters:\n")
cat("Expected ratio of maternal to paternal reads: ", mp_ratio,"\n")
cat("p-value threshold: ", pval_cutoff,"\n")
cat("alpha cutoff: ", alpha,"\n")
cat("Min. number of allelic reads in each dir of cross: ", minallelic,"\n")
cat("----------------------\n")
cat("\n")

# check input files exist
if (! file.exists(AxB_file)) {
	stop("Error: input file ",AxB_file," could not be opened")
}
if (! file.exists(BxA_file)) {
	stop("Error: input file ",BxA_file," could not be opened")
}

# read in input files
AxB_data <- read.table(AxB_file, header=FALSE, sep="\t", stringsAsFactors = FALSE)
BxA_data <- read.table(BxA_file, header=FALSE, sep="\t", stringsAsFactors = FALSE)

# check that columns 2 and 3 of both files can be coerced to numeric
if (suppressWarnings(all(!is.na(as.numeric(AxB_data[,2])))) == TRUE) { AxB_data[,2] = as.numeric(AxB_data[,2]) } else { stop("Error: column 2 of AxB input file contains non-numeric values") }		
if (suppressWarnings(all(!is.na(as.numeric(AxB_data[,3])))) == TRUE) { AxB_data[,3] = as.numeric(AxB_data[,3]) } else { stop("Error: column 3 of AxB input file contains non-numeric values") }	
if (suppressWarnings(all(!is.na(as.numeric(BxA_data[,2])))) == TRUE) { BxA_data[,2] = as.numeric(BxA_data[,2]) } else { stop("Error: column 2 of BxA input file contains non-numeric values") }	
if (suppressWarnings(all(!is.na(as.numeric(BxA_data[,3])))) == TRUE) { BxA_data[,3] = as.numeric(BxA_data[,3]) } else { stop("Error: column 3 of BxA input file contains non-numeric values") }	

# check that same genes in same order provided in both input files
genelist_AxB <- AxB_data[,1]
genelist_BxA <- BxA_data[,1]

if(length(genelist_AxB) != length(genelist_BxA)){
	stop(sprintf("AxB and BxA lists have different lengths: %d and %d\n", length(genelist_AxB), length(genelist_BxA)))
}
for(i in 1:length(genelist_AxB)){
	if(genelist_AxB[i] != genelist_BxA[i]){
		stop(paste("AxB and BxA should contain loci in same order; encountered different locus at same index: ", i, sep=""))
	}
}
genelist <- genelist_AxB
totgenes <- length(genelist)

# combine lists of counts for both crosses
colnames(AxB_data) = c("locus_name","AxB_A_counts","AxB_B_counts")
colnames(BxA_data) = c("locus_name","BxA_A_counts","BxA_B_counts")
AxB_data$AxB_tot_counts = AxB_data$AxB_A_counts + AxB_data$AxB_B_counts
BxA_data$BxA_tot_counts = BxA_data$BxA_A_counts + BxA_data$BxA_B_counts
alldata = cbind(AxB_data,BxA_data[,-1])

# mask sites with fewer than minallelic counts in either direction of the cross
alldata$mincounts_mask <- ifelse(alldata$AxB_tot_counts >= minallelic & alldata$BxA_tot_counts >= minallelic,TRUE, FALSE)


# -------------------------
# Test for imprinting
# -------------------------
# calculate initial p-value using Fisher's exact test
cat(" - Testing genes for significant parental bias\n")
pvals = array(data=NA, dim=totgenes)
fav_parent = array(data=NA, dim=totgenes)
bias_factor = array(data=NA, dim=totgenes)

for (i in 1:totgenes) {
	if (alldata$mincounts_mask[i] == TRUE) {
		res = test_for_imprinting(alldata$AxB_A_counts[i], alldata$AxB_B_counts[i], alldata$BxA_A_counts[i], alldata$BxA_B_counts[i], mp_ratio, alpha) 
		pvals[i] = res$pval
		bias_factor[i] = as.numeric(sprintf("%.3f", res$bias_factor))

		# translate favored dir to favored parent
		if (res$favored_dir == 'X') {
			fav_parent[i] <- 'mother'
		} else if (res$favored_dir == 'Y') {
			fav_parent[i] <- 'father'
		} else {
			fav_parent[i] = 'none'
		}
	} else {
		fav_parent[i] = 'none'
		pvals[i] = 1.0
		bias_factor[i] = NA
	}
}

alldata$imprinting_pval = pvals
alldata$favored_parent = fav_parent
alldata$IFs = bias_factor

# adjust p-values for multiple testing
alldata$imprinting_padj = array(NA, dim=totgenes)
alldata$imprinting_padj[alldata$mincounts_mask] = p.adjust(alldata$imprinting_pval[alldata$mincounts_mask], method="BH")
alldata$imprinting_padj[!alldata$mincounts_mask] = 1.0

# mask imprinting factor if p-value >= cutoff
alldata$IFs = ifelse(alldata$imprinting_padj < pval_cutoff,alldata$IFs,'NA')	


# -------------------------
# Test for cis-effects
# -------------------------
# essentially the same protocol as testing for imprinting, except we flip the order
# of BxA_A and BxA_B provided to test_for_imprinting, and expected ratio is always 1
cat(" - Testing genes for significant strain bias\n")
pvals = array(data=NA, dim=totgenes)
fav_strain = array(data=NA, dim=totgenes)
bias_factor = array(data=NA, dim=totgenes)

for (i in 1:totgenes) {
	if (alldata$mincounts_mask[i] == TRUE) {
		res = test_for_imprinting(alldata$AxB_A_counts[i], alldata$AxB_B_counts[i], alldata$BxA_B_counts[i], alldata$BxA_A_counts[i], 1, alpha) 
		pvals[i] = res$pval
		bias_factor[i] = as.numeric(sprintf("%.3f", res$bias_factor))

		# translate favored dir to favored strain
		if (res$favored_dir == 'X') {
			fav_strain[i] <- 'A'
		} else if (res$favored_dir == 'Y') {
			fav_strain[i] <- 'B'
		} else {
			fav_strain[i] = 'none'
		}
	} else {
		fav_strain[i] = 'none'
		pvals[i] = 1.0
		bias_factor[i] = NA
	}
}

alldata$cis_effect_pval = pvals
alldata$favored_strain = fav_strain
alldata$CEF = bias_factor

# adjust p-values for multiple testing
alldata$cis_effect_padj = array(NA, dim=totgenes)
alldata$cis_effect_padj[alldata$mincounts_mask] = p.adjust(alldata$cis_effect_pval[alldata$mincounts_mask], method="BH")
alldata$cis_effect_padj[!alldata$mincounts_mask] = 1.0

# mask imprinting factor if p-value >= cutoff
alldata$CEF = ifelse(alldata$cis_effect_padj < pval_cutoff,alldata$CEF,'NA')	


# -------------------------
# Output results
# -------------------------
# get subset of fields to output
alldata$mp_ratio = array(data=mp_ratio, dim=totgenes)
output = alldata[ , c("locus_name","AxB_A_counts","AxB_B_counts","BxA_A_counts","BxA_B_counts","mp_ratio","favored_parent","imprinting_padj","IFs","favored_strain","cis_effect_padj","CEF")]
write.table(output, outfile, sep = "\t", quote = FALSE, row.names = FALSE)









