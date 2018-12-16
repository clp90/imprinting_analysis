#!/usr/bin/env Rscript
if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
if (is.element('plyr', installed.packages()[,1])==FALSE) { stop("Error: R package plyr could not be loaded") }
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(plyr))

# ------------------------------------------------------------------------------------
# v2.0 by Colette L. Picard
# 11/08/2018
# ------------------------------------------------------------------------------------

# Usage:
# filter_imprinted.R [options] infile.txt outprefix
#

# Description:
# -------------------------
# Given the output of detect_imprinted.R, further filters the obtained genes
# with additional thresholds. Also creates plot of % maternal in AxB vs BxA.

# infile.txt = the output of detect_imprinting.R
# outprefix = a prefix (incl. path) for output files (whole dataset, nonzero subset, impr lists, plot)

# OVERVIEW:
# (1) Filter provided gene list according to all filters provided
# (2) Summarize results (genes passing each successive filter)
# (3) Make scatterplots of MEGs and PEGs
# (4) Output results for all genes, genes passing minallelic cutoffs only, and output list of MEG, list of PEGs

# Additional options:
# - User can provide a list of genes to filter out of the final list of MEGs/PEGs (e.g. genes whose
# imprinting status cannot be determined because they are highly expressed in seed coat)
# - By default, no filtering on CEF occurs; use --CEF to set a max cis-effects factor value for MEGs/PEGs

# -------------------------
# Version history:
# v.1.0: initial build - 04/10/2015
# v.1.1: (05/18/2015) added option to filter out a list of genes from the analysis (e.g. seedcoat)
# and added ability to add additional information to the output files (see --addinfo option)
# v.2.0: (11/11/2018) somewhat major edit, particularly to varnames to bring syntax in line with rest of
# imprinting suite (e.g. AxB instead of REFxOTH), minor clean-up, dropped --addinfo option
# and added minallelic option. Switched to argparse from optparse for parsing command-line options.
# -------------------------
parser = ArgumentParser()
parser$add_argument("infile", nargs=1, help = "an output file from detect_imprinted.R, tab-delimited txt")
parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")
parser$add_argument("--ratio", type="double", default = 1, help = "expected ratio of maternal to paternal reads (1 for most tissues, 2 for endosperm)")
parser$add_argument("--minallelic", type="integer", default = 1, help = "Min. number of allelic reads in each direction of cross in order to evaluate imprinting")
parser$add_argument("--pval", type="double", default = 0.01, help = "Maximum (adjusted) p-value for imprinted locus, default 0.01")
parser$add_argument("--IF", type="double", default = 2, help = "Minimum imprinting factor (IF) for imprinted locus, default 2")
parser$add_argument("--CEF", type="double", default = -1, help = "Maximum cis-effects factor (CEF) for imprinted locus, default no filtering on CEF")
parser$add_argument("--MEG", type="double", default = 70.0, help = "Minimum % maternal reads for a MEG, default 70")
parser$add_argument("--PEG", type="double", default = 30.0, help = "Minimum % maternal reads for a MEG, default 30")
parser$add_argument("--nameA", type="character", default = "strainA", help = "Name of strain A")
parser$add_argument("--nameB", type="character", default = "strainB", help = "Name of strain B")
parser$add_argument("--filter", type="character", default = "", help = "List of genes to filter out of analysis (e.g. high expressed seedcoat)")
parser$add_argument("--pointsize", type="double", default = 0.25, help = "Size of points (adjust to taste), default 0.25")

args <- commandArgs(trailingOnly = TRUE)

# the -0 option is for checking that all dependencies are installed; exit gracefully here if so
if (length(args) == 1 & args[1] == "-0") {
	q(save="no", status=0, runLast=FALSE)
}

if (length(args) <= 1) {
	cat("Usage: filter_imprinted.R [options] infile.txt outprefix\n")
	cat("----------------------\n")
	cat("infile.txt : an output file from detect_imprinted.R\n")
	cat("outprefix : prefix to use for all output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--ratio <int> : ratio of expected maternal to paternal expression (e.g. 1 for embryo, 2 for endosperm) - default 1\n")
	cat("--minallelic <int> : min number of allelic reads in both AxB and BxA to evaluate - default 1\n")
	cat("--pval <float> : p value cutoff for imprinted locus - default 0.01\n")
	cat("--IF <float> : Minimum imprinting factor (IF) for imprinted locus - default 2\n")
	cat("--CEF <float> : Maximum cis-effects factor (CEF) for imprinted locus - default none\n")
	cat("--MEG <float> : Minimum % maternal reads for a MEG - default 70\n")
	cat("--PEG <float> : Maximum % maternal reads for a PEG - default 30\n")
	cat("--nameA <str> : Name of reference strain - default strainA\n")
	cat("--nameB <str> : Name of alternate strain - default strainB\n")
	cat("--filter <str> : Path to file containing list of genes to filter out of analysis (e.g. seedcoat) - default none\n")
	cat("--pointsize <str> : Size of points (adjust to taste) - default 0.25\n")
	cat("----------------------\n")
	stop("At least one required argument not provided. See usage above.")
}

opt <- parser$parse_args()
infile = opt$infile
outprefix = opt$outprefix

# filters used in this script
cat("----------------------\n")
cat("Input file: ",infile,"\n")
cat("Prefix for output files: ",outprefix,"\n")
cat("Reference strain is:",opt$nameA,"\n")
cat("Alternate strain is:",opt$nameB,"\n")
cat("----------------------\n")
cat("Applying the following filters to the imprinted genes:\n")
cat("Maximum p-value cutoff: ",opt$pval,"\n")
cat("Minimum imprinting factor (IF): ",opt$IF,"\n")
if (opt$CEF == -1) { 
	cat("Not filtering on cis-effects (CEF)\n")
} else {
	cat("Maximum cis-effects factor (CEF): ",opt$CEF,"\n")
}
cat("Minimum % maternal reads for MEGs: ",opt$MEG,"\n")
cat("Maximum % maternal reads for PEGs: ",opt$PEG,"\n")
if (! opt$filter == "") { cat("Filtering out all genes in: ",opt$filter,"\n") }
cat("----------------------\n")

# check input file exists
if (! file.exists(infile)) {
	stop("Error: input file ",infile," could not be opened")
}

alldata = read.table(infile, header=TRUE, sep="\t", stringsAsFactors = FALSE)
frac_mat_meg = opt$MEG / 100
frac_mat_peg = opt$PEG / 100

# ----------------------
# get total counts and percent maternal for all genes
# ----------------------
alldata$tot_AxB = alldata$AxB_A_counts + alldata$AxB_B_counts		# total counts in the RefxAlt cross
alldata$tot_BxA = alldata$BxA_A_counts + alldata$BxA_B_counts		# total counts in the AltxRef cross

alldata$pmat_AxB = alldata$AxB_A_counts / alldata$tot_AxB			# percent maternal in RefxAlt
alldata$pmat_BxA = alldata$BxA_B_counts / alldata$tot_BxA			# percent maternal in AltxRef

# ----------------------
# create masks for the different filters above
# ----------------------
# minallelic filter
alldata$mincountspass = ifelse(alldata$tot_AxB >= opt$minallelic & alldata$tot_BxA > opt$minallelic,TRUE, FALSE)	
alldata$status = ifelse(alldata$mincountspass == TRUE, "TBD", "low_counts")

# if provided, censor all genes in --filter
if (! opt$filter == "") {
	to_filt = data.frame(read.table(opt$filter, header=FALSE, sep = '\t'))
	names(to_filt)[names(to_filt)=="V1"] <- "locus_name"
	to_filt$censored = as.logical(rep(1,length(to_filt$locus_name)))	

	# merge the gene IDs to be dropped with the dataset
	alldata = merge(alldata,to_filt,by="locus_name", all.x = TRUE)
	alldata$censored[is.na(alldata$censored)] = FALSE	
} else {
	alldata$censored = FALSE
}
alldata$status = ifelse(alldata$status == "TBD" & alldata$censored == TRUE, "censored", alldata$status)

# p-value filter
alldata$pval_pass = ifelse(alldata$imprinting_padj < opt$pval,TRUE, FALSE)	
alldata$status = ifelse(alldata$status == "TBD" & alldata$pval_pass == FALSE, "fail_pval_cutoff", alldata$status)

# imprinting factor cutoff
alldata$IF_pass = ifelse(alldata$IFs >= opt$IF & !is.na(alldata$IFs),TRUE, FALSE)
alldata$status = ifelse(alldata$status == "TBD" & alldata$IF_pass == FALSE, "fail_IF_cutoff", alldata$status)

# cis-effects factor cutoff
if (opt$CEF == -1) { 
	alldata$CEF_pass = as.logical(rep(1,length(opt$CEF)))
} else {
	alldata$CEF_pass = ifelse(alldata$CEF < opt$CEF | is.na(alldata$CEF),TRUE, FALSE)
	alldata$status = ifelse(alldata$status == "TBD" & alldata$CEF_pass == FALSE, "fail_CEF_cutoff", alldata$status)
}

# percent maternal cutoff - MEGs
alldata$pmat_pass = ifelse(alldata$pmat_AxB >= frac_mat_meg & alldata$pmat_BxA >= frac_mat_meg & !is.na(alldata$pmat_BxA) & !is.na(alldata$pmat_AxB) & alldata$favored_parent == "mother", TRUE, FALSE)
alldata$pmat_pass = ifelse(alldata$pmat_AxB <= frac_mat_peg & alldata$pmat_BxA <= frac_mat_peg & !is.na(alldata$pmat_BxA) & !is.na(alldata$pmat_AxB) & alldata$favored_parent == "father", TRUE, alldata$pmat_pass)
alldata$status = ifelse(alldata$status == "TBD" & alldata$pmat_pass == FALSE, "fail_pmat_cutoff", alldata$status)

# imprinted genes that pass all the cutoffs above
# (note - all genes that still have status "TBD" passed all the cutoffs above)
alldata$isMEG = ifelse(alldata$favored_parent == "mother" & alldata$status == "TBD", TRUE, FALSE)
alldata$isPEG = ifelse(alldata$favored_parent == "father" & alldata$status == "TBD", TRUE, FALSE)
alldata$status = ifelse(alldata$isMEG == TRUE, "MEG", alldata$status)
alldata$status = ifelse(alldata$isPEG == TRUE, "PEG", alldata$status)

# move status column to the end
alldata = alldata[c(setdiff(names(alldata), "status"), "status")]


# ----------------------
# Summarize results
# ----------------------
cat("Summarizing filtering results:\n")
cat("Total number of loci in file:",length(alldata$locus_name),"\n")
num_eval = nrow(alldata[alldata$mincountspass == TRUE,])
cat(" - ",num_eval," genes had allelic counts >= ",opt$minallelic," in both crosses and could be evaluated\n",sep="")

if (! opt$filter == "") {
	num_censored = nrow(alldata[alldata$status == "censored",]); f_censored = num_censored / num_eval * 100
	cat(sprintf(" - Of these, %d (%.1f%%) were in the --filter file and censored\n", num_censored, f_censored))
	rem = num_eval - num_censored; num_eval = rem
	cat(" - Of the remaining",num_eval,"loci evaluated:\n")
} else {
	cat(" - Of these ",num_eval,":\n",sep="")
}

num_pval_fail = nrow(alldata[alldata$status == "fail_pval_cutoff",]); f_pval_fail = num_pval_fail / num_eval * 100
cat(sprintf("   - %d (%.1f%%) failed the p-value cutoff (adjusted imprinting p-value >= %s)\n", num_pval_fail, f_pval_fail, as.character(opt$pval)))

num_IF_fail = nrow(alldata[alldata$status == "fail_IF_cutoff",]); f_IF_fail = num_IF_fail / num_eval * 100
cat(sprintf("   - %d (%.1f%%) passed the p-value cutoff but failed the IF cutoff (IF < %s)\n", num_IF_fail, f_IF_fail, as.character(opt$IF)))

if (opt$CEF != -1) {
	num_CEF_fail = nrow(alldata[alldata$status == "fail_CEF_cutoff",]); f_CEF_fail = num_CEF_fail / num_eval * 100
	cat(sprintf("   - %d (%.1f%%) passed the p-value and IF cutoffs but failed the CEF cutoff (CEF >= %s)\n", num_CEF_fail, f_CEF_fail, as.character(opt$CEF)))
}

num_pmat_fail = nrow(alldata[alldata$status == "fail_pmat_cutoff",]); f_pmat_fail = num_pmat_fail / num_eval * 100
cat(sprintf("   - %d (%.1f%%) passed the p-value and IF cutoffs but failed the pmat cutoff (pmat < %s for MEGs, or pmat > %s for PEGs)\n", num_pmat_fail, f_pmat_fail, as.character(opt$MEG), as.character(opt$PEG)))

num_imprinted = nrow(alldata[(alldata$status == "MEG" | alldata$status == "PEG"),]); f_imprinted = num_imprinted / num_eval * 100
cat(sprintf("   - %d (%.1f%%) passed all filters for imprinting\n", num_imprinted, f_imprinted))
cat('\n')
cat("----------------------\n")
cat("Number of MEGs passing all filters:",nrow(alldata[alldata$status == "MEG", ]),"\n")
cat("Number of PEGs passing all filters:",nrow(alldata[alldata$status == "PEG", ]),"\n")
cat("----------------------\n")
cat('\n')


# ----------------------
# output whole dataset with filters, as well as the "minallelic" subset, and list of MEGs and PEGs
# ----------------------
minallelic = alldata[alldata$mincountspass == TRUE, ]
MEGs = alldata[alldata$isMEG == TRUE,1]
PEGs = alldata[alldata$isPEG == TRUE,1]

write.table(alldata, paste(outprefix,"_all.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(minallelic, paste(outprefix,"_mincounts.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(MEGs, paste(outprefix,"_MEGs.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(PEGs, paste(outprefix,"_PEGs.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ----------------------
# make plot of log(MAT/PAT) for AxB vs BxA
# ----------------------
# add pseudocounts to allow plotting loci where counts for denom = 0
alldata$log_mp_AxB = log2((alldata$AxB_A_counts+1) / (alldata$AxB_B_counts+1))
alldata$log_mp_BxA = log2((alldata$BxA_B_counts+1) / (alldata$BxA_A_counts+1))

# color points based on if they're MEGs (red) or PEGs (blue) else gray
alldata$color = "ivory4"
alldata$color[alldata$isMEG==TRUE] = "indianred3"
alldata$color[alldata$isPEG==TRUE] = "dodgerblue2"

# only plot loci with nonzero totals for both crosses
nomiss = subset(alldata, alldata$mincountspass == TRUE)
nomiss = nomiss[with(nomiss, order(isMEG, isPEG)), ]		# plot imprinted genes on top

# scatterplot
png(paste(outprefix,'_plot_all_minallelic.png',sep=""), width = 6, height = 6.5, units = 'in', res = 300)
plot(nomiss$log_mp_AxB, nomiss$log_mp_BxA,
	xlab = "", sub = paste("log2 (m/p reads)\n",opt$nameA,"x",opt$nameB,"\n"), 
	mtext(paste("\nlog2 (m/p reads)\n",opt$nameB,"x",opt$nameA,"\n"), 2, line=1), ylab = "",
	pch = 20, cex=opt$pointsize, ylim=c(-10,10), xlim=c(-10,10), col=nomiss$color)
abline(h = log2(opt$ratio), col = "gray")
abline(v = log2(opt$ratio), col = "gray")
abline(a = 0, b = 1, col = "gray", lty=2)
graphics.off()

if (! opt$filter == "") { 
	# repeat with another color showing genes dropped b/c of filt list, and also make plot
	nomiss$color[nomiss$censored==TRUE] = "yellowgreen" 
	
	png(paste(outprefix,'_plot_censored_in_green.png',sep=""), width = 6, height = 6.5, units = 'in', res = 300)
	plot(nomiss$log_mp_AxB, nomiss$log_mp_BxA,
		xlab = "", sub = paste("log2 (m/p reads)\n",opt$nameA,"x",opt$nameB,"\n"), 
		mtext(paste("\nlog2 (m/p reads)\n",opt$nameB,"x",opt$nameA,"\n"), 2, line=1), ylab = "",
		pch = 20, cex=opt$pointsize, ylim=c(-10,10), xlim=c(-10,10), col=nomiss$color)
	abline(h = log2(opt$ratio), col = "gray")
	abline(v = log2(opt$ratio), col = "gray")
	abline(a = 0, b = 1, col = "gray", lty=2)
	graphics.off()
	
	nocensored = subset(nomiss, nomiss$censored == FALSE)
	
	png(paste(outprefix,'_plot_no_censored.png',sep=""), width = 6, height = 6.5, units = 'in', res = 300)
	plot(nocensored$log_mp_AxB, nocensored$log_mp_BxA,
		xlab = "", sub = paste("log2 (m/p reads)\n",opt$nameA,"x",opt$nameB,"\n"), 
		mtext(paste("\nlog2 (m/p reads)\n",opt$nameB,"x",opt$nameA,"\n"), 2, line=1), ylab = "",
		pch = 20, cex=opt$pointsize, ylim=c(-10,10), xlim=c(-10,10), col=nocensored$color)
	abline(h = log2(opt$ratio), col = "gray")
	abline(v = log2(opt$ratio), col = "gray")
	abline(a = 0, b = 1, col = "gray", lty=2)
	graphics.off()	
}























