#!/bin/bash

# ------------------------------------------------------------------------------------
# v2.1 by Colette L. Picard
# 12/17/2018
# ------------------------------------------------------------------------------------

# Usage:
# call_imprinted.sh [options] -1 cross1.bam -2 cross2.bam -s snp_file.bed -o outdir

# -------------------------
# Version history:
# v.1.0: initial build - 04/30/2015
# v.2.0: major rebuild - 11/08/2018
# v.2.1: added -N option (htseq-count option --nonunique all); changed default to --nonunique none - 12/17/2018
# -------------------------

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
v.2.1 by Colette L Picard, 12/17/2018
-------------------------------
Usage:
call_imprinted.sh [options] -1 cross1.bam -2 cross2.bam -s snp_file.bed -R ratio -o outdir
or
call_imprinted.sh [options] -x AxB_A_counts.txt -y AxB_B_counts.txt -X BxA_A_counts.txt -Y BxA_B_counts.txt -R ratio -o outdir
-------------------------------

The goal of this script is to identify imprinted genes, starting either from (1) read
alignments (BAM) or (2) per-gene allelic read counts. Imprinted genes are defined as
genes that are preferentially expressed from either the maternal or paternal allele, and
are usually identified using RNA-seq data from reciprocal crosses between two accessions
or species, so that polymorphisms can be used to distinguish parent-of-origin of sequencing
reads. 

OVERVIEW:
-----if starting from BAM-----
(1) assign mapped reads to parent of origin using SNPs
(2) count allelic reads at each gene using htseq-count
(3) perform preliminary imprinting assessment for each gene using counts
(4) filter preliminary list of imprinted genes to get final list
-----if starting from counts-----
(1) perform preliminary imprinting assessment for each gene using counts
(2) filter preliminary list of imprinted genes to get final list

-------------------------------
ADDITIONAL BACKGROUND:
This script assumes that the data derive from a pair of reciprocal crosses between
strain A and strain B:

	A (mother) x B (father)
	B (mother) x A (father)
	
For each gene, the number of reads that derive from strains A and B (determined using SNPs)
is counted, producing four counts per gene. Genes are imprinted if counts from the maternal
parent are significantly higher in both directions of the cross. In contrast, non-imprinted
genes are expressed approximately equally from the maternal and paternal alleles. Reciprocal
data are also important to be able to distinguish imprinting from strain effects - that is,
if the strain A allele of a gene is more highly expressed due to a genetic or epigenetic
difference, then counts the A allele will be higher regardless of the direction of the cross.
	
	A allele, AxB (M)	B allele, AxB (F)	A allele, BxA (F)	B allele, BxA (M)
	----------			----------			----------			----------			
	high				high				high				high				-> not imprinted
	low					low					low					low					-> not imprinted
	----------			----------			----------			----------			
	high				low					high				low					-> strain-specific (A allele >> B allele)
	low					high				low					high				-> strain-specific (B allele >> A allele)
	----------			----------			----------			----------			
	high				low					low					high				-> imprinted (maternally expressed - MEG)
	low					high				high				low					-> imprinted (paternally expressed - PEG)
	
-------------------------------
ADDITIONAL CONSIDERATIONS:
- The required argument -R is the expected ratio of maternal to paternal reads
for non-imprinted genes in the tissue assayed. Although in most tissues this will be 1,
imprinting in flowering plants usually occurs in the endosperm, which is a triploid
tissue with two maternal genomes and one paternal genome. That means that the expected ratio 
of maternal to paternal reads for most endosperms is 2, not 1.
- If starting from BAM files:
	(1) BAM files can also be (uncompressed) SAM files; either is fine
	(2) you must also provide a file containing a list of SNPs in BED format, with column 1 
		= chromosome, column 2 = SNP position (BED files are 0-based, meaning the first position 
		in the genome is position 1), column 3 = SNP position + 1, and 
		column 4 = (strain A allele)>(strain B allele). For example:
			scaffold_1	9630	9631	T>C
			scaffold_1	11142	11143	G>A
			scaffold_1	13391	13392	A>G
	(3) you must also provide a GTF file containing annotations for known transcripts; this 
		is used by htseq-count to count allelic reads in each gene.
	(4) if alignments are from a paired-end run, you may want to filter out discordant read pairs
	(e.g. keep only proper pairs) since it is unclear how to treat discordant pairs when assigning to 
	parent-of-origin. Normally for a proper pair, if one mate overlaps a SNP, the entire pair is assigned
	to that parent-of-origin even if the other mate has no SNPs, but whether that behavior should
	extent to discordant pairs (e.g. very far apart, wrong orientation, etc.) is unclear.
	By default, discordant pairs are treated as two singletons, but this might lead to double counting
	(since each singleton counts as 1). Alternately, discordant pairs can be treated like proper pairs 
	if the -l flag is used.
	(5) strongly recommend stranded RNA-seq data if possible. By default, htseq-count is in --nonunique none
	mode, so reads overlapping more than one annotation (on the same strand, if RNA-seq data is stranded)
	will not be counted. This is commonly an issue for non-stranded data, where a read may map to a locus
	with both a sense and antisense annotation - since it is unclear which transcript generated the read,
	it is censored. However, this behavior can be turned off by using -N with this script (equivalent to
	the htseq-count option --nonunique all) - with this option, reads overlapping multiple annotations
	will be counted for all of them. Note that this is generally not recommended.
- The location for output files is provided with -o and must be a non-existant or empty
folder, else the script will (by default) refuse to overwrite its contents. Note this can be overridden
with the -r option, which you use at your own risk.
-------------------------------

Required and optional arguments are listed below. Default values, when applicable, are given
at the end of the line in square brackets. Required programs and scripts also listed.

User-specified options:
Required arguments:
	-o outdir : name of directory for all output files
Required arguments for starting from BAM:
	-1 AxB : bam or sam file with tophat-aligned reads for the AxB cross (where mother is A)
	-2 BxA : bam or sam file with tophat-aligned reads for the BxA cross (where mother is B)
	-S SNPs : bed file containing SNPs for allele calling
	-G annot : (htseq-count) GTF/GFF annotation file of known transcripts
Required arguments for starting from counts:
	-x AxB_A : counts file for A reads, from AxB
	-y AxB_B : counts file for B reads, from AxB
	-X BxA_A : counts file for A reads, from BxA
	-Y BxA_B : counts file for B reads, from BxA
Additional options:
	-n name : prefix, used for some output files ["expt"]
	-R ratio : expected ratio of maternal to paternal expression (e.g. 2 for endosperm, 1 for embryo) [1]
	-A strainA : name of strain A ["strainA"]
	-B strainB : name of strain B ["strainB"]
	-m htseq_mode : htseq-count mode (see https://htseq.readthedocs.io/en/master/count.html) - must be "union" "intersection_strict" or "intersection_nonempty" (note that --nonunique all is used) ["union"]
	-c minallelic : (filter_imprinted) minimum number of allelic reads required in each direction of the cross [5]
	-p pval : (filter_imprinted) p-value cutoff for gene to be called imprinted (Fishers Exact Test, BH-corrected) [0.01]
	-I minIF : (filter_imprinted) minimum imprinting factor (IF) for gene to be called imprinted (see Gehring et al. 2011) [2]
	-C maxCEF : (filter_imprinted) maximum cis-effects factor (CEF) for gene to be called imprinted []
	-M minpmMEG : (filter_imprinted) minimum % maternal reads for maternally expressed imprinted genes (MEGs) [70]
	-P maxpmPEG : (filter_imprinted) maximum % maternal reads for paternally expressed imprinted genes (PEGs) [30]
	-f filter : (filter_imprinted) list of genes to filter out of analysis (e.g. high seedcoat expression) [""]
	-i pointsize : (filter_imprinted) size to plot points in scatterplots, adjust up or down to taste (or just re-run filter_imprinted directly) [0.25]
	-s path_to_scripts : path to folder containing all required scripts (note - $scriptDir in default is the location of this script) ["$scriptDir/scripts"]
	-a memalloc : memory (in Gb) to allocate to java when running MarkDuplicates [10]
Flag options:
	-l : (assign_to_allele) treat all read pairs as pairs for assign_to_allele and counting purposes (default proper pairs only, non-proper pairs treated as singletons) [relaxed=false]
	-D : (assign_to_allele) remove PCR duplicates after assigning reads to parent-of-origin [dedup=false]
	-F : (htseq-count) only run htseq-count on the reads that could be assigned to a parent of origin (default also runs htseq-count on all reads; use if time is an issue) [skipall=false]
	-N : (htseq-count) for reads overlapping multiple annotations, count for all overlapping annotations (by default, these reads are censored; see htseq-count --nonunique all) [alluniq=false]
	-w : (htseq-count) AxB library is stranded with first read in pair on same strand as feature [AxBfirststrand=false]
	-W : (htseq-count) BxA library is stranded with first read in pair on same strand as feature [BxAfirststrand=false]
	-v : (htseq-count) AxB library is stranded with first read in pair on opposite strand as feature [AxBsecstrand=false]
	-V : (htseq-count) BxA library is stranded with first read in pair on opposite strand as feature [BxAsecstrand=false]
	-r : allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!) [overwrite=false]
	-0 : checks that all required programs installed on PATH and all required helper scripts can be located, then exits without running
	-h : prints this version and usage information

Required installed on user PATH:
	- python (tested on v2.7.6)
	- R (tested on v3.4.4)
	- samtools (by Heng Li, tested on v1.7)
	- java (tested on v.1.7.0_181, OpenJDK Runtime Environment (IcedTea 2.6.14) (7u181-2.6.14-0ubuntu0.2), OpenJDK 64-Bit Server VM (build 24.181-b01, mixed mode))
	
Required python modules:
	- htseq (by Simon Anders, tested on v0.9.1)

Must be in path_to_scripts (if not on default path, specify location with -s):
	- assign_to_allele.py (by Colette Picard, v.1.4)
	- detect_imprinting.R (by Victor Missirian, edits by Elias Scheer and Colette Picard, v.2.0)
	- filter_imprinted.R (by Colette Picard, v.1.2)
	- merge_by_column.R (by Colette Picard, v.1.2)
	- MarkDuplicates.jar (from the picard-tools suite, Broad Institute, downloaded Oct. 2, 2015)
		
------------------------------------------------------------------------------------
EOF

[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# ----------------------
# Get user-specified arguments
# ----------------------

# Initiate environment
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 
workdir=$( pwd )											# working directory

# Required arguments:
# ----------------------
outdir=""							# name of directory for all output files

# Required arguments for starting from BAM:
# ----------------------
AxB=""							# bam or sam file with tophat-aligned reads for the AxB cross (where mother is A)
BxA=""							# bam or sam file with tophat-aligned reads for the BxA cross (where mother is B)
SNPs=""							# bed file containing SNPs for allele calling
annot=""							# (htseq-count) GTF/GFF annotation file of known transcripts

# Additional options:
# ----------------------
name="expt"								# prefix, used for some output files
ratio="1"							# expected ratio of maternal to paternal expression (e.g. 2 for endosperm, 1 for embryo)
strainA="strainA"							# name of strain A
strainB="strainB"							# name of strain B
htseq_mode="union"						# htseq-count mode (see https://htseq.readthedocs.io/en/master/count.html) - must be "union" "intersection_strict" or "intersection_nonempty" (note that --nonunique all is used)
minallelic=5							# (filter_imprinted) minimum number of allelic reads required in each direction of the cross
pval=0.01								# (filter_imprinted) p-value cutoff for gene to be called imprinted (Fishers Exact Test, BH-corrected)
minIF=2								# (filter_imprinted) minimum imprinting factor (IF) for gene to be called imprinted (see Gehring et al. 2011)
maxCEF=""							# (filter_imprinted) maximum cis-effects factor (CEF) for gene to be called imprinted
minpmMEG=70								# (filter_imprinted) minimum % maternal reads for maternally expressed imprinted genes (MEGs)
maxpmPEG=30								# (filter_imprinted) maximum % maternal reads for paternally expressed imprinted genes (PEGs)
filter=""							# (filter_imprinted) list of genes to filter out of analysis (e.g. high seedcoat expression)
pointsize=0.25							# (filter_imprinted) size to plot points in scatterplots, adjust up or down to taste (or just re-run filter_imprinted directly) 
path_to_scripts="$scriptDir/helper_scripts"							# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
memalloc=10							# memory (in Gb) to allocate to java when running MarkDuplicates

# Required arguments for starting from counts:
# ----------------------
AxB_A=""							# counts file for A reads, from AxB
AxB_B=""							# counts file for B reads, from AxB
BxA_A=""							# counts file for A reads, from BxA
BxA_B=""							# counts file for B reads, from BxA

# Flag options:
# ----------------------
relaxed=false							# (assign_to_allele) treat all read pairs as pairs for assign_to_allele and counting purposes (default proper pairs only, non-proper pairs treated as singletons)
dedup=false							# (assign_to_allele) remove PCR duplicates after assigning reads to parent-of-origin
skipall=false							# (htseq-count) only run htseq-count on the reads that could be assigned to a parent of origin (default also runs htseq-count on all reads; use if time is an issue)
alluniq=false							# (htseq-count) for reads overlapping multiple annotations, count for all overlapping annotations (by default, these reads are censored; see htseq-count --nonunique all)
AxBfirststrand=false							# (htseq-count) AxB library is stranded with first read in pair on same strand as feature
BxAfirststrand=false							# (htseq-count) BxA library is stranded with first read in pair on same strand as feature
AxBsecstrand=false							# (htseq-count) AxB library is stranded with first read in pair on opposite strand as feature
BxAsecstrand=false							# (htseq-count) BxA library is stranded with first read in pair on opposite strand as feature
overwrite=false							# allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!)

checkdep=false

# ----------------------
while getopts "R:o:1:2:S:G:x:y:X:Y:n:A:B:m:c:p:I:C:M:P:f:i:s:a:lDFNwWvVr0h" opt; do
	case $opt in
		R)	# expected ratio of maternal to paternal expression (e.g. 2 for endosperm, 1 for embryo)
			ratio="$OPTARG"
			;;
		o)	# name of directory for all output files
			outdir="$OPTARG"
			;;
		1)	# bam or sam file with tophat-aligned reads for the AxB cross (where mother is A)
			AxB="$OPTARG"
			;;
		2)	# bam or sam file with tophat-aligned reads for the BxA cross (where mother is B)
			BxA="$OPTARG"
			;;
		S)	# bed file containing SNPs for allele calling
			SNPs="$OPTARG"
			;;
		G)	# (htseq-count) GTF/GFF annotation file of known transcripts
			annot="$OPTARG"
			;;
		x)	# counts file for A reads, from AxB
			AxB_A="$OPTARG"
			;;
		y)	# counts file for B reads, from AxB
			AxB_B="$OPTARG"
			;;
		X)	# counts file for A reads, from BxA
			BxA_A="$OPTARG"
			;;
		Y)	# counts file for B reads, from BxA
			BxA_B="$OPTARG"
			;;
		n)	# prefix, used for some output files
			name="$OPTARG"
			;;
		A)	# name of strain A
			strainA="$OPTARG"
			;;
		B)	# name of strain B
			strainB="$OPTARG"
			;;
		m)	# htseq-count mode (see https://htseq.readthedocs.io/en/master/count.html) - must be "union" "intersection_strict" or "intersection_nonempty" (note that --nonunique all is used)
			htseq_mode="$OPTARG"
			;;
		c)	# (filter_imprinted) minimum number of allelic reads required in each direction of the cross
			minallelic="$OPTARG"
			;;
		p)	# (filter_imprinted) p-value cutoff for gene to be called imprinted (Fishers Exact Test, BH-corrected)
			pval="$OPTARG"
			;;
		I)	# (filter_imprinted) minimum imprinting factor (IF) for gene to be called imprinted (see Gehring et al. 2011)
			minIF="$OPTARG"
			;;
		C)	# (filter_imprinted) maximum cis-effects factor (CEF) for gene to be called imprinted
			maxCEF="$OPTARG"
			;;
		M)	# (filter_imprinted) minimum % maternal reads for maternally expressed imprinted genes (MEGs)
			minpmMEG="$OPTARG"
			;;
		P)	# (filter_imprinted) maximum % maternal reads for paternally expressed imprinted genes (PEGs)
			maxpmPEG="$OPTARG"
			;;
		f)	# (filter_imprinted) list of genes to filter out of analysis (e.g. high seedcoat expression)
			filter="$OPTARG"
			;;
		i)	# (filter_imprinted) size to plot points in scatterplots, adjust up or down to taste (or just re-run filter_imprinted directly) 
			pointsize="$OPTARG"
			;;
		s)	# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
			path_to_scripts="$OPTARG"
			;;
		a)	# memory (in Gb) to allocate to java when running MarkDuplicates
			memalloc="$OPTARG"
			;;
		l)	# (assign_to_allele) treat all read pairs as pairs for assign_to_allele and counting purposes (default proper pairs only, non-proper pairs treated as singletons)
			relaxed=true
			;;
		D)	# (assign_to_allele) remove PCR duplicates after assigning reads to parent-of-origin
			dedup=true
			;;
		F)	# (htseq-count) only run htseq-count on the reads that could be assigned to a parent of origin (default also runs htseq-count on all reads; use if time is an issue)
			skipall=true
			;;
		N)	# (htseq-count) for reads overlapping multiple annotations, count for all overlapping annotations (by default, these reads are censored; see htseq-count --nonunique all)
			alluniq=true
			;;
		w)	# (htseq-count) AxB library is stranded with first read in pair on same strand as feature
			AxBfirststrand=true
			;;
		W)	# (htseq-count) BxA library is stranded with first read in pair on same strand as feature
			BxAfirststrand=true
			;;
		v)	# (htseq-count) AxB library is stranded with first read in pair on opposite strand as feature
			AxBsecstrand=true
			;;
		V)	# (htseq-count) BxA library is stranded with first read in pair on opposite strand as feature
			BxAsecstrand=true
			;;
		r)	# allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!)
			overwrite=true
			;;
		0)	# check dependencies ok then exit
			checkdep=true
			;;
		h)	# print usage and version information to stdout and exit
			echo "$usage"
			exit 0
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

# ----------------------
# Helper functions for this script
# ----------------------
err_msg ()
# prints an error message both to stdout and to the log file, then exits
# usage: err_msg msg logfile
{
	printf "Error: $1 \n"
	printf "Exited due to error: $1 \n" >> "$2"
	exit 1	
}

displaytime () {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  [[ $D > 0 ]] && printf '%d days ' $D
  [[ $H > 0 ]] && printf '%d hours ' $H
  [[ $M > 0 ]] && printf '%d minutes ' $M
  [[ $D > 0 || $H > 0 || $M > 0 ]] && printf 'and '
  printf '%d seconds\n' $S
}

compress_sam ()
# sorts and compresses a SAM file to BAM and then indexes it; removes the original SAM file
# Usage: compress_sam infile scriptlog
{
	infile="$1"; logfile="$2"
	# grab filename minus extension
	ff="${infile%.*}"
	if [[ -f "$infile" && -s "$infile" ]]; then
		samtools sort "$infile" -o "${ff}.bam" -T "${ff}"
		[ $? != 0 ] && err_msg "samtools sort failed for input file $infile" "$logfile"
		samtools index "${ff}.bam"
		[ $? != 0 ] && err_msg "samtools index failed for input file ${ff}.bam" "$logfile"
		rm "$infile"
	fi
}

# ----------------------
# Main code
# ----------------------

# Check that all programs required on PATH are installed
# ----------------------
command -v python >/dev/null 2>&1 || { echo "Error: python is required on PATH but was not found"; exit 1; }
command -v R >/dev/null 2>&1 || { echo "Error: R is required on PATH but was not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools is required on PATH but was not found"; exit 1; }
command -v java >/dev/null 2>&1 || { echo "Error: java is required on PATH but was not found"; exit 1; }

# Check that required python modules are installed
# ----------------------
python -c "import HTSeq" 2>/dev/null; [ $? != 0 ] && { echo "Error: required python module htseq not installed"; exit 1; }

# Check that all other required helper scripts can be located
# ----------------------
[ -f "$path_to_scripts/assign_to_allele.py" ] || { echo "Error: required helper script assign_to_allele.py was not found in $path_to_scripts (see -s option)"; exit 1; }
[ -f "$path_to_scripts/detect_imprinting.R" ] || { echo "Error: required helper script detect_imprinting.R was not found in $path_to_scripts (see -s option)"; exit 1; }
[ -f "$path_to_scripts/filter_imprinted.R" ] || { echo "Error: required helper script filter_imprinted.R was not found in $path_to_scripts (see -s option)"; exit 1; }
[ -f "$path_to_scripts/merge_by_column.R" ] || { echo "Error: required helper script merge_by_column.R was not found in $path_to_scripts (see -s option)"; exit 1; }
[ -f "$path_to_scripts/MarkDuplicates.jar" ] || { echo "Error: required helper script MarkDuplicates.jar was not found in expected location $path_to_scripts"; exit 1; }

# Check that all R scripts have all their dependencies installed
# ----------------------
$path_to_scripts/detect_imprinting.R -0; [ $? != 0 ] && exit 1
$path_to_scripts/filter_imprinted.R -0; [ $? != 0 ] && exit 1
$path_to_scripts/filter_imprinted.R -0; [ $? != 0 ] && exit 1

# Done checking all requirements. Stop here if -0 flagged.
# ----------------------
"$checkdep" && exit 0

# Check all required inputs are provided
# ----------------------
[ -z "$ratio" ] && { echo "Error: -R ratio is a required argument (expected ratio of maternal to paternal expression (e.g. 2 for endosperm, 1 for embryo))"; exit 1; }
[ -z "$outdir" ] && { echo "Error: -o outdir is a required argument (name of directory for all output files)"; exit 1; }

# If any inputs provided for 'starting from counts', check all inputs provided
# ----------------------
if [[ ! -z "$AxB_A" || ! -z "$AxB_B" || ! -z "$BxA_A" || ! -z "$BxA_B" ]]; then
	mode="counts"
	[ -z "$AxB_A" ] && { echo "Error: if starting from counts, -x AxB_A is a required argument (counts file for A reads, from AxB)"; exit 1; }
	[ -z "$AxB_B" ] && { echo "Error: if starting from counts, -y AxB_B is a required argument (counts file for B reads, from AxB)"; exit 1; }
	[ -z "$BxA_A" ] && { echo "Error: if starting from counts, -X BxA_A is a required argument (counts file for A reads, from BxA)"; exit 1; }
	[ -z "$BxA_B" ] && { echo "Error: if starting from counts, -Y BxA_B is a required argument (counts file for B reads, from BxA)"; exit 1; }
	# check that provided input files can be opened
	[ -f "$AxB_A" ] || { echo "Error: could not open AxB, A counts file $AxB_A"; exit 1; }
	[ -f "$AxB_B" ] || { echo "Error: could not open AxB, B counts file $AxB_B"; exit 1; }
	[ -f "$BxA_A" ] || { echo "Error: could not open BxA, A counts file $BxA_A"; exit 1; }
	[ -f "$BxA_B" ] || { echo "Error: could not open BxA, B counts file $BxA_B"; exit 1; }
	# check that SAM/BAM, SNP and annotation files NOT provided
	[ ! -z "$AxB" ] && { echo "Error: if starting from counts, -1 AxB should not be used"; exit 1; }
	[ ! -z "$BxA" ] && { echo "Error: if starting from counts, -2 BxA should not be used"; exit 1; }
	[ ! -z "$SNPs" ] && { echo "Error: if starting from counts, -S SNPs should not be used"; exit 1; }
	[ ! -z "$annot" ] && { echo "Error: if starting from counts, -G annot should not be used"; exit 1; }	
fi

# If any inputs provided for 'starting from BAM', check all inputs provided
# ----------------------
if [[ ! -z "$AxB" || ! -z "$BxA" || ! -z "$SNPs" || ! -z "$annot" ]]; then
	mode="BAM"
	[ -z "$AxB" ] && { echo "Error: if starting from SAM/BAM, -1 AxB is a required argument (bam or sam file with tophat-aligned reads for the AxB cross (where mother is A))"; exit 1; }
	[ -z "$BxA" ] && { echo "Error: if starting from SAM/BAM, -2 BxA is a required argument (bam or sam file with tophat-aligned reads for the BxA cross (where mother is B))"; exit 1; }
	[ -z "$SNPs" ] && { echo "Error: if starting from SAM/BAM, -S SNPs is a required argument (bed file containing SNPs for allele calling)"; exit 1; }
	[ -z "$annot" ] && { echo "Error: if starting from SAM/BAM, -G annot is a required argument ((htseq-count) GTF/GFF annotation file of known transcripts)"; exit 1; }
	# check that provided input files can be opened
	[ -f "$AxB" ] || { echo "Error: could not open AxB input file $AxB"; exit 1; }
	[ -f "$BxA" ] || { echo "Error: could not open BxA input file $BxA"; exit 1; }
	[ -f "$SNPs" ] || { echo "Error: Error: could not open SNP file $SNPs"; exit 1; }
	[ -f "$annot" ] || { echo "Error: Error: could not open GTF annotation file $annot"; exit 1; }
	# check that count files NOT provided
	[ ! -z "$AxB_A" ] && { echo "Error: if starting from SAM/BAM, -x AxB_A should not be used"; exit 1; }
	[ ! -z "$AxB_B" ] && { echo "Error: if starting from SAM/BAM, -y AxB_B should not be used"; exit 1; }
	[ ! -z "$BxA_A" ] && { echo "Error: if starting from SAM/BAM, -X BxA_A should not be used"; exit 1; }
	[ ! -z "$BxA_B" ] && { echo "Error: if starting from SAM/BAM, -Y BxA_B should not be used"; exit 1; }
	# get library type (ignored if counts)
	AxBtype="no"
	[[ "$AxBfirststrand" = "true" && "$AxBsecstrand" = "true" ]] && { echo "Error: cannot use both -w and -v at the same time (pick one)"; exit 1; }
	[ "$AxBfirststrand" = "true" ] && AxBtype="yes"
	[ "$AxBsecstrand" = "true" ] && AxBtype="reverse"
	BxAtype="no"
	[[ "$BxAfirststrand" = "true" && "$BxAsecstrand" = "true" ]] && { echo "Error: cannot use both -w and -v at the same time (pick one)"; exit 1; }
	[ "$BxAfirststrand" = "true" ] && BxAtype="yes"
	[ "$BxAsecstrand" = "true" ] && BxAtype="reverse"
fi

# Check if neither SAM/BAM nor counts files were provided
# ----------------------
[[ -z "$AxB_A" && -z "$AxB" ]] && { echo "Error: please provide either SAM/BAM files or count files"; exit 1; }

# Check if both SAM/BAM and counts files were provided
# ----------------------
[[ ! -z "$AxB_A" && ! -z "$AxB" ]] && { echo "Error: please provide either SAM/BAM files or count files"; exit 1; }

# Check that htseq-count mode provided is allowed
# ----------------------
if [[ "$htseq_mode" != "union" && "$htseq_mode" != "intersection_strict" && "$htseq_mode" != "intersection_nonempty" ]]; then
	echo "Error: value for -m must be 'union', 'intersection_strict' or 'intersection_nonempty'"
	exit 1
fi

# If output folder doesn't exist, make it, else overwrite if user used -r, else error
# ----------------------
if [ ! -d "$outdir" ]; then
	mkdir "$outdir"
else
	if [ "$(ls -A $outdir)" ]; then
		if [ "$overwrite" = "true" ]; then
			echo "Overwriting previous contents of output dir $outdir"
			rm -rf "$outdir"/*	
		else
			echo "Error: provided output directory is not empty. To allow overwrite of non-empty dir, use -r flag. WARNING: all existing files in -outdir- will be deleted. Seriously."
			exit 1
		fi
	fi
fi


log="$outdir/${name}_log.txt" 	# create log file
time_start=$(date)				# time run was started
time_ss=$(date +%s)				# time run was started (in seconds)

# Output user-derived options to stdout and to log file
# ----------------------
echo "Running call_imprinted.sh v2.0 (11/08/2018):" | tee "$log"
echo "Run start on: $time_start" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Working directory: $( pwd )" | tee -a "$log"
echo "Output directory: $outdir" | tee -a "$log"
echo "Log file: $log" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
if [ ! -z "$AxB_A" ]; then
	echo "Starting from counts files:" | tee -a "$log"
	echo "AxB, A counts file: $AxB_A" | tee -a "$log"
	echo "AxB, B counts file: $AxB_B" | tee -a "$log"
	echo "BxA, A counts file: $BxA_A" | tee -a "$log"
	echo "BxA, B counts file: $BxA_B" | tee -a "$log"
else
	echo "Starting from SAM/BAM alignment files:" | tee -a "$log"
	echo "AxB alignments file: $AxB" | tee -a "$log"
	[ "$AxBtype" = "yes" ] && { echo "AxB library is stranded with first read in pair on same strand as feature" | tee -a "$log"; }
	[ "$AxBtype" = "reverse" ] && { echo "AxB library is stranded with first read in pair on opposite strand as feature" | tee -a "$log"; }
	[ "$AxBtype" = "no" ] && { echo "AxB library is not stranded" | tee -a "$log"; }
	echo "BxA alignments file: $BxA" | tee -a "$log"
	[ "$BxAtype" = "yes" ] && { echo "BxA library is stranded with first read in pair on same strand as feature" | tee -a "$log"; }
	[ "$BxAtype" = "reverse" ] && { echo "BxA library is stranded with first read in pair on opposite strand as feature" | tee -a "$log"; }
	[ "$BxAtype" = "no" ] && { echo "BxA library is not stranded" | tee -a "$log"; }
	echo "List of A>B SNPs: $SNPs" | tee -a "$log"
	echo "GTF or GFF annotation file: $annot" | tee -a "$log"
	[ "$relaxed" = "true" ] || { echo "When assigning reads to parent-of-origin, counting discordant pairs as pairs instead of singletons" | tee -a "$log"; }
	[ "$relaxed" = "false" ] || { echo "When assigning reads to parent-of-origin, counting discordant pairs as singletons instead of pairs" | tee -a "$log"; }
	[ "$dedup" = "true" ] || { echo "PCR duplicates will be removed with MarkDuplicates before getting allele-specific counts with htseq-count" | tee -a "$log"; }
	[ "$skipall" = "true" ] || { echo "Only performing htseq-count on allelic reads" | tee -a "$log"; }
fi
echo "-------------------------" | tee -a "$log"
echo "Prefix for output files: $name" | tee -a "$log"
echo "Strain A name: $strainA" | tee -a "$log"
echo "Strain B name: $strainB" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Imprinted genes must pass the following filters:" | tee -a "$log"
echo "Minimum number of allelic reads required to assess imprinting (in both directions of cross): $minallelic" | tee -a "$log"
echo "p-value (Fisher's exact test, BH-corrected) < $pval" | tee -a "$log"
echo "Imprinting Factor (IF) >= $minIF" | tee -a "$log"
[ -z "$maxCEF" ] && { echo "No cis-effects factor cutoff" | tee -a "$log"; } || { echo "Cis-effects Factor (CEF) <= $maxCEF" | tee -a "$log"; }
echo "For MEGs, % maternal must be >= $minpmMEG in both AxB and BxA" | tee -a "$log"
echo "For PEGs, % maternal must be <= $maxpmPEG in both AxB and BxA" | tee -a "$log"
[ -z "$filter" ] || { echo "Additionally, genes in this list will be excluded from the analysis: $filter" | tee -a "$log"; }
echo "-------------------------" | tee -a "$log"
echo "Additional parameters:" | tee -a "$log"
[ -z "$addinfo" ] || { echo "Additionally, genes in this list will be excluded from the analysis: $filter" | tee -a "$log"; }
echo "Htseq-count is being run in mode: $htseq_mode" | tee -a "$log"
[ "$alluniq" = "true" ] && { echo "Htseq-count will count reads overlapping multiple annotations once for each annotation (WARNING: not recommended)" | tee -a "$log"; }
echo "Location of helper scripts: $path_to_scripts" | tee -a "$log"
echo "-------------------------" | tee -a "$log"


# ----------------------
# Step 0: check input files ok
# ----------------------
printf "\nStep 0: Checking that all inputs are ok\n" | tee -a "$log"

if [ "$mode" = "BAM" ]; then
	# SNP file should have four columns, no header, etc.
	res=$( awk -F$'\t' '{OFS=FS} {if (NF != 4) {print "COLERR"} else if ($2 != $3-1) {print "POSERR"} else if ($4 !~ /^[ATGC]>[ATGC]$/) {print "SNPERR"}}' "$SNPs" )
	if [ ! -z "$res" ]; then
		echo "Error(s) detected in SNP file:"
		rres=$( head -1 $res )
		[ "$rres" = "COLERR" ] && echo "SNP file must have exactly four columns (chr, 0-based position, pos+1, SNP)"
		[ "$rres" = "POSERR" ] && echo "SNP column 3 must be equal to 1 + value of column 2 (0-based position of SNP)"
		[ "$rres" = "SNPERR" ] && echo "4th column (SNP) must be in the form strain1allele>strain2allele (e.g. A>T)"
		echo "First few lines of file:"
		head -3 "$SNPs"
		exit 1
	fi
fi


# ----------------------
# Step 1: convert both input files to SAM, sort and assign to parent-of-origin
# ----------------------
if [ "$mode" = "BAM" ]; then
	printf "\nStep 1: assigning mapped reads to parent of origin using SNPs\n" | tee -a "$log"; ts=$(date +%s)

	mkdir "$outdir/assign_to_allele"
	
	# sort both input SAM/BAM files and convert to SAM if needed
	echo " - Sorting input files by position" | tee -a "$log"
		
	samtools sort "$AxB" -o "$outdir/${name}_AxB.sam" -T "$outdir/AxB_tmp" > /dev/null 2>&1
	samtools sort "$BxA" -o "$outdir/${name}_BxA.sam" -T "$outdir/BxA_tmp" > /dev/null 2>&1

# -----------ASSIGN_TO_ALLELE-----------
	echo " - Assigning reads to parent of origin using SNPs" | tee -a "$log"
	[ "$relaxed" = "true" ] && relstr=" --relaxed" || relstr=""
	$path_to_scripts/assign_to_allele.py "$SNPs" "$outdir/${name}_AxB.sam" "$outdir/assign_to_allele/${name}_${strainA}x${strainB}" --refname "$strainA" --altname "$strainB"${relstr} > "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_log.txt"
	[ $? != 0 ] && err_msg "error running assign_to_allele, see $outdir/assign_to_allele/${name}_${strainA}x${strainB}_log.txt" "$log"

	$path_to_scripts/assign_to_allele.py "$SNPs" "$outdir/${name}_BxA.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}" --refname "$strainA" --altname "$strainB"${relstr} > "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_log.txt"
	[ $? != 0 ] && err_msg "error running assign_to_allele, see $outdir/assign_to_allele/${name}_${strainB}x${strainA}_log.txt" "$log"
	
	# summarize results
	echo "" | tee -a "$log"
	echo "Summary of results (indicated percentages are of pairs/singletons assessed):" | tee -a "$log"
	echo "  For reads from ${strainA}x${strainB}:" | tee -a "$log"
	awk '$0 ~ /assigned/ {print "    "$0}' "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_log.txt" | tee -a "$log"
	echo "  For reads from ${strainB}x${strainA}:" | tee -a "$log"
	awk '$0 ~ /assigned/ {print "    "$0}' "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_log.txt" | tee -a "$log"

# -----------REMOVE PCR DUPLICATES-----------
	if [ "$dedup" = "true" ]; then
		echo "" | tee -a "$log"
		echo " - Removing presumed PCR duplicates from allele-specific alignments" | tee -a "$log"
		if [ "$skipall" = "false" ]; then
			samfiles=( "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainB}.sam" "$outdir/${name}_AxB.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainB}.sam" "$outdir/${name}_BxA.sam" )
		else
			samfiles=( "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainB}.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainB}.sam" )
		fi
		
		mkdir -p "$outdir/assign_to_allele/MarkDuplicates_logs" "$outdir/assign_to_allele/BAM_files_with_PCR_duplicates_included" "$outdir/assign_to_allele/flagstat_summaries"
		
		for ((i=0;i<${#samfiles[@]};++i)); do
			ff="${samfiles[i]}"
			bb="${ff%.*}"
			base=$(basename "$ff")	# filename without path
			bbase="${base%.*}"		# filename without path or extension
			
			# sort and compress SAM files
			compress_sam "$ff"
			mv "${bb}.bam" "$outdir/assign_to_allele/BAM_files_with_PCR_duplicates_included/${bbase}.bam"
			mv "${bb}.bam.bai" "$outdir/assign_to_allele/BAM_files_with_PCR_duplicates_included/${bbase}.bam.bai"
			
			# run MarkDuplicates
			java -Xmx${memalloc}g -jar $path_to_scripts/MarkDuplicates.jar I="$outdir/assign_to_allele/BAM_files_with_PCR_duplicates_included/${bbase}.bam" O="${bb}.bam" METRICS_FILE="$outdir/assign_to_allele/MarkDuplicates_logs/MarkDuplicates_metrics_${bbase}.txt" REMOVE_DUPLICATES=true > "$outdir/assign_to_allele/MarkDuplicates_logs/MarkDuplicates_log_${bbase}.txt" 2>&1
			[ $? != 0 ] && err_msg "MarkDuplicates failed, see $outdir/assign_to_allele/MarkDuplicates_logs/MarkDuplicates_log_${bbase}.txt" "$log"
			
			# summarize result
			samtools flagstat "${bb}.bam" > "$outdir/assign_to_allele/flagstat_summaries/${bbase}_summary.txt"
			
			# convert to SAM and delete BAM
			samtools view -ho "${bb}.sam" "${bb}.bam"
			rm "${bb}.bam"			
		done	
		
		echo "" | tee -a "$log"
		echo "Summary of results after removing PCR duplicates (indicated percentages are of pairs/singletons assessed):" | tee -a "$log"
		lloc="$outdir/assign_to_allele/flagstat_summaries"
		AxB_A_p=$( awk '$0 ~ /itself and mate mapped/ {print $1 / 2} ' "$lloc/expt_${strainA}x${strainB}_${strainA}_summary.txt" )
		AxB_B_p=$( awk '$0 ~ /itself and mate mapped/ {print $1 / 2} ' "$lloc/expt_${strainA}x${strainB}_${strainB}_summary.txt" )
		BxA_A_p=$( awk '$0 ~ /itself and mate mapped/ {print $1 / 2} ' "$lloc/expt_${strainB}x${strainA}_${strainA}_summary.txt" )
		BxA_B_p=$( awk '$0 ~ /itself and mate mapped/ {print $1 / 2} ' "$lloc/expt_${strainB}x${strainA}_${strainB}_summary.txt" )

		AxB_A_s=$( awk '$0 ~ /singletons/ {print $1 / 2} ' "$lloc/expt_${strainA}x${strainB}_${strainA}_summary.txt" )
		AxB_B_s=$( awk '$0 ~ /singletons/ {print $1 / 2} ' "$lloc/expt_${strainA}x${strainB}_${strainB}_summary.txt" )
		BxA_A_s=$( awk '$0 ~ /singletons/ {print $1 / 2} ' "$lloc/expt_${strainB}x${strainA}_${strainA}_summary.txt" )
		BxA_B_s=$( awk '$0 ~ /singletons/ {print $1 / 2} ' "$lloc/expt_${strainB}x${strainA}_${strainB}_summary.txt" )
	
		echo "	For reads from ${strainA}x${strainB}:" | tee -a "$log"
		echo "		# pairs assigned to ${strainA}: $AxB_A_p" | tee -a "$log"
		echo "		# pairs assigned to ${strainB}: $AxB_B_p" | tee -a "$log"
		echo "		# singletons assigned to ${strainA}: $AxB_A_s" | tee -a "$log"
		echo "		# singletons assigned to ${strainB}: $AxB_B_s" | tee -a "$log"
	
		echo "	For reads from ${strainB}x${strainA}:" | tee -a "$log"
		echo "		# pairs assigned to ${strainA}: $BxA_A_p" | tee -a "$log"
		echo "		# pairs assigned to ${strainB}: $BxA_B_p" | tee -a "$log"
		echo "		# singletons assigned to ${strainA}: $BxA_A_s" | tee -a "$log"
		echo "		# singletons assigned to ${strainB}: $BxA_B_s" | tee -a "$log"
	
	fi

	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"
	

# ----------------------
# Step 2: count allelic reads at each gene using htseq-count
# ----------------------
	printf "\nStep 2: counting allelic reads at each gene using htseq-count\n" | tee -a "$log"; ts=$(date +%s)

	mkdir "$outdir/counts_per_gene"
	
	# note - htseq-count assumes BAM files are sorted by name so that entries of paired reads are adjacent;
	# output of assign_to_allele should have all pairs adjacent already so no need for additional sort
	if [ "$skipall" = "false" ]; then
		samfiles=( "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainB}.sam" "$outdir/${name}_AxB.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainB}.sam" "$outdir/${name}_BxA.sam" )
		libtypes=( "$AxBtype" "$AxBtype" "$AxBtype" "$BxAtype" "$BxAtype" "$BxAtype" )
		cross=( "${strainA}x${strainB}" "${strainA}x${strainB}" "${strainA}x${strainB}" "${strainB}x${strainA}" "${strainB}x${strainA}" "${strainB}x${strainA}" )
		strain=( "$strainA" "$strainB" "all" "$strainA" "$strainB" "all" )
	else
		samfiles=( "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainB}.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainA}.sam" "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainB}.sam" )
		libtypes=( "$AxBtype" "$AxBtype" "$BxAtype" "$BxAtype" )
		cross=( "${strainA}x${strainB}" "${strainA}x${strainB}" "${strainB}x${strainA}" "${strainB}x${strainA}" )
		strain=( "$strainA" "$strainB" "$strainA" "$strainB" )
	fi
	
# -----------HTSEQ-COUNT-----------
	[ "$dedup" = "true" ] && sortstr=" --order pos" || sortstr="" 
	[ "$alluniq" = "true" ] && uniqstr="all" || uniqstr="none"
	for ((i=0;i<${#samfiles[@]};++i)); do
		echo " - Counting ${strain[i]} reads in ${cross[i]}" | tee -a "$log"
		if [ -f "${samfiles[i]}" ]; then		
			python -m HTSeq.scripts.count --nonunique "$uniqstr"${sortstr} -m "$htseq_mode" -s "${libtypes[i]}" "${samfiles[i]}" "$annot" > "$outdir/counts_per_gene/${name}_${cross[i]}_${strain[i]}_counts.txt" 2> "$outdir/counts_per_gene/htseq_count_${cross[i]}_${strain[i]}_log.txt"
			[ $? != 0 ] && err_msg "error running htseq-count, see $outdir/counts_per_gene/htseq_count_${cross[i]}_${strain[i]}_log.txt" "$log"	
		else
			err_msg "no ${strain[i]} reads detected in ${cross[i]}" "$log"
		fi
	done
	
	# compress or remove SAM files that are no longer needed
	echo " - Cleaning up and compressing SAM files" | tee -a "$log"
	compress_sam "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainA}.sam"
	compress_sam "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_${strainB}.sam"
	compress_sam "$outdir/assign_to_allele/${name}_${strainA}x${strainB}_none.sam"

	compress_sam "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainA}.sam"
	compress_sam "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_${strainB}.sam"
	compress_sam "$outdir/assign_to_allele/${name}_${strainB}x${strainA}_none.sam"
	
	# assign count files to variables to match set-up in "counts" mode
	AxB_A="$outdir/counts_per_gene/${name}_${strainA}x${strainB}_${strainA}_counts.txt"							# counts file for A reads, from AxB
	AxB_B="$outdir/counts_per_gene/${name}_${strainA}x${strainB}_${strainB}_counts.txt"							# counts file for B reads, from AxB
	BxA_A="$outdir/counts_per_gene/${name}_${strainB}x${strainA}_${strainA}_counts.txt"							# counts file for A reads, from BxA
	BxA_B="$outdir/counts_per_gene/${name}_${strainB}x${strainA}_${strainB}_counts.txt"							# counts file for B reads, from BxA
	
	# done with these, can safely remove
	rm -f "$outdir/${name}_AxB.sam" "$outdir/${name}_BxA.sam"

	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"
	
fi


# ----------------------
# Step 3: perform preliminary imprinting assessment for each gene using counts
# ----------------------
echo "" | tee -a "$log"
[ "$mode" = "BAM" ] && echo "Step 3: performing preliminary imprinting assessment for each gene using counts" | tee -a "$log"
[ "$mode" = "counts" ] && echo "Step 1: performing preliminary imprinting assessment for each gene using counts" | tee -a "$log"
ts=$(date +%s)

mkdir -p "$outdir/counts_per_gene"

# combine the two counts files for each cross into one file - geneID countsA countsB
echo " - Combining counts files for each cross"
echo "geneID	countsA" > "$outdir/counts_per_gene/${strainA}x${strainB}_${strainA}_tmp.txt"
awk -F$'\t' '{if ($1 !~ /^__/) {print $0}}' "${AxB_A}" >> "$outdir/counts_per_gene/${strainA}x${strainB}_${strainA}_tmp.txt"		# htseq-count outputs lines starting with __ that are summary, not gene counts; delete
echo "geneID	countsB" > "$outdir/counts_per_gene/${strainA}x${strainB}_${strainB}_tmp.txt"
awk -F$'\t' '{if ($1 !~ /^__/) {print $0}}' "${AxB_B}" >> "$outdir/counts_per_gene/${strainA}x${strainB}_${strainB}_tmp.txt"
echo "geneID	countsA" > "$outdir/counts_per_gene/${strainB}x${strainA}_${strainA}_tmp.txt"
awk -F$'\t' '{if ($1 !~ /^__/) {print $0}}' "${BxA_A}" >> "$outdir/counts_per_gene/${strainB}x${strainA}_${strainA}_tmp.txt"
echo "geneID	countsB" > "$outdir/counts_per_gene/${strainB}x${strainA}_${strainB}_tmp.txt"
awk -F$'\t' '{if ($1 !~ /^__/) {print $0}}' "${BxA_B}" >> "$outdir/counts_per_gene/${strainB}x${strainA}_${strainB}_tmp.txt"

$path_to_scripts/merge_by_column.R "$outdir/counts_per_gene/${strainA}x${strainB}_${strainA}_tmp.txt" "$outdir/counts_per_gene/${strainA}x${strainB}_${strainB}_tmp.txt" geneID "$outdir/counts_per_gene/${name}_${strainA}x${strainB}_tmp.txt" > "$outdir/counts_per_gene/${name}_${strainA}x${strainB}_merge_log.txt"
$path_to_scripts/merge_by_column.R "$outdir/counts_per_gene/${strainB}x${strainA}_${strainA}_tmp.txt" "$outdir/counts_per_gene/${strainB}x${strainA}_${strainB}_tmp.txt" geneID "$outdir/counts_per_gene/${name}_${strainB}x${strainA}_tmp.txt" > "$outdir/counts_per_gene/${name}_${strainB}x${strainA}_merge_log.txt"

tail -n+2 "$outdir/counts_per_gene/${name}_${strainA}x${strainB}_tmp.txt" > "$outdir/counts_per_gene/${name}_${strainA}x${strainB}_counts_merged.txt"	# remove header
tail -n+2 "$outdir/counts_per_gene/${name}_${strainB}x${strainA}_tmp.txt" > "$outdir/counts_per_gene/${name}_${strainB}x${strainA}_counts_merged.txt"	# remove header

rm $outdir/counts_per_gene/*_tmp.txt

# -----------DETECT_IMPRINTING-----------
echo " - Performing preliminary imprinting analysis"
mkdir "$outdir/imprinting"
$path_to_scripts/detect_imprinting.R "$outdir/counts_per_gene/${name}_${strainA}x${strainB}_counts_merged.txt" "$outdir/counts_per_gene/${name}_${strainB}x${strainA}_counts_merged.txt" "$outdir/imprinting/${name}_imprinting_calls.txt" --mp_ratio 2 --minallelic "$minallelic" > "$outdir/imprinting/${name}_imprinting_calls_log.txt"
[ $? != 0 ] && err_msg "error running detect_imprinting.R, see $outdir/imprinting/${name}_imprinting_calls_log.txt" "$log"

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"


# ----------------------
# Step 4: filter preliminary imprinting calls to obtain final set of imprinted genes
# ----------------------
echo "" | tee -a "$log"
[ "$mode" = "BAM" ] && echo "Step 4: filtering preliminary imprinting calls to obtain final set of imprinted genes" | tee -a "$log"
[ "$mode" = "counts" ] && echo "Step 2: filtering preliminary imprinting calls to obtain final set of imprinted genes" | tee -a "$log"
ts=$(date +%s)

[ -z "$maxCEF" ] && cefstr="" || cefstr=" --CEF $maxCEF"
[ -z "$filter" ] && filtstr="" || filtstr=" --filter $filter"

# -----------FILTER_IMPRINTED-----------
$path_to_scripts/filter_imprinted.R "$outdir/imprinting/${name}_imprinting_calls.txt" "$outdir/imprinting/${name}_imprinting_filtered" --ratio "$ratio" --nameA "$strainA" --nameB "$strainB" --minallelic "$minallelic" --pval "$pval" --IF "$minIF" --MEG "$minpmMEG" --PEG "$maxpmPEG" --pointsize "$pointsize"${cefstr}${filtstr} > "$outdir/imprinting/${name}_imprinting_filtering_log.txt"
[ $? != 0 ] && err_msg "error running filter_imprinted.R, see $outdir/imprinting/${name}_imprinting_filtering_log.txt" "$log"

echo "" | tee -a "$log"
echo "Summary of results:" | tee -a "$log"

[ -z "$filter" ] && loglines=13 || loglines=14
tail -${loglines} "$outdir/imprinting/${name}_imprinting_filtering_log.txt" | tee -a "$log"

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"

echo "" | tee -a "$log"
time_end=$(date)	# time run ended
time_es=$(date +%s)	# time run ended
echo "Run ended $time_end" | tee -a "$log"
echo "Total time elapsed: $( displaytime $(($time_es - $time_ss)) )" | tee -a "$log"





