#!/bin/bash

# ------------------------------------------------------------------------------------
# v1.0 by Colette L. Picard
# 11/13/2018
# ------------------------------------------------------------------------------------
# This script runs a few basic tests to ensure that rna_seq_map.sh, call_imprinting.sh,
# get_homology.py and make_metagenome.py can run successfully and that all their required 
# dependencies are installed.

start=$(date)
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )		# location of this script
echo "--------------------"
echo "Testing installation of all scripts from the github repository imprinting_analysis (2018)"
echo "Run date and time: $start"
echo "Location of currently running script: $scriptDir"
echo "--------------------"

# ---------------------
# Check all dependencies
# ---------------------
echo "Testing installation of all scripts and dependencies"

command -v python >/dev/null 2>&1 || { echo "Error: python is required on PATH but was not found"; exit 1; }
pyversion=$( python -c 'import sys; print(sys.version_info[0])' )
[ "$pyversion" -ne 2 ] && { echo "Error: the default python version (given by python --version) must be 2.x.x"; exit 1; }

[ -f "$scriptDir/comp_imprinting.py" ] || { echo "Error: expected comp_imprinting.py in the same directory as this script"; exit 1; }
echo "comp_imprinting.py - OK"

[ -f "$scriptDir/get_homologs.py" ] || { echo "Error: expected get_homologs.py in the same directory as this script"; exit 1; }
python -c "import intermine" 2>/dev/null; [ $? != 0 ] && { echo "Error: required python module intermine not installed (required for get_homologs.py)"; exit 1; }
echo "get_homologs.py - OK"

[ -f "$scriptDir/make_metagenome.py" ] || { echo "Error: expected make_metagenome.py in the same directory as this script"; exit 1; }
python -c "import Bio" 2>/dev/null; [ $? != 0 ] && { echo "Error: required python module Bio not installed (required for make_metagenome.py)"; exit 1; }
echo "make_metagenome.py - OK"

[ -f "$scriptDir/rna_seq_map.sh" ] || { echo "Error: expected rna_seq_map.sh in the same directory as this script"; exit 1; }
$scriptDir/rna_seq_map.sh -0
[ $? = 0 ] && echo "rna_seq_map.sh - OK" || exit 1

[ -f "$scriptDir/call_imprinting.sh" ] || { echo "Error: expected call_imprinting.sh in the same directory as this script"; exit 1; }
$scriptDir/call_imprinting.sh -0
[ $? = 0 ] && echo "call_imprinting.sh - OK" || exit 1

echo ""
echo "All dependencies OK"
echo ""

# ---------------------
# Run tests on each main script to make sure it's working correctly
# ---------------------

# check all directories and files exist (they should have been copied verbatim from the github repository)
testdir="$scriptDir/tests"
testinputs="$testdir/test_inputs"
[ -d "$testdir" ] || { echo "Error: could not find directory ${testdir}"; exit 1; }
[ -d "$testinputs" ] || { echo "Error: could not find directory ${testinputs}"; exit 1; }
expected="$testdir/expected"; mkdir -p "$expected"

# Run tests on get_homologs
# ---------------------
echo "Trying a test run of get_homologs.py to make sure everything is working correctly..."
testoutdir="$scriptDir/tests/get_homologs"; mkdir -p "$testoutdir"
exp="$expected/get_homologs"

[ -f "$testinputs/genelist.txt" ] || { echo "Error: could not find file $testinputs/genelist.txt"; exit 1; } 

$scriptDir/get_homologs.py "A. thaliana" "A. lyrata" "$testinputs/genelist.txt" "$testoutdir/homolog_list.txt" > "$testoutdir/log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with get_homologs.py, see $testoutdir/log.txt"; exit 1; }

# check output file matches
cmp --silent "$testoutdir/homolog_list.txt" "$exp/homolog_list.txt" || { echo "Error: output file $testoutdir/homolog_list.txt doesn't match expected output $exp/homolog_list.txt"; exit 1; }
echo "get_homologs.py - all OK"


# Run tests on make_metagenome
# ---------------------
echo "Trying a few test runs of make_metagenome.py to make sure everything is working correctly..."
testoutdir="$scriptDir/tests/make_metagenome"; mkdir -p "$testoutdir"
exp="$expected/make_metagenome"

[ -f "$testinputs/chr1_subset.fa" ] || { echo "Error: could not find file $testinputs/chr1_subset.fa"; exit 1; } 
[ -f "$testinputs/GTF_subset.gtf" ] || { echo "Error: could not find file $testinputs/GTF_subset.gtf"; exit 1; } 
[ -f "$testinputs/snps_subset.bed" ] || { echo "Error: could not find file $testinputs/snps_subset.bed"; exit 1; } 
[ -f "$testinputs/snps_subset_w.bed" ] || { echo "Error: could not find file $testinputs/snps_subset.bed"; exit 1; } 

[ -f "$exp/out.fa" ] || { echo "Error: could not find file $exp/out.fa"; exit 1; } 
[ -f "$exp/out_metachrom.txt" ] || { echo "Error: could not find file $exp/out_metachrom.txt"; exit 1; } 
[ -f "$exp/out_metagtf.gtf" ] || { echo "Error: could not find file $exp/out_metagtf.gtf"; exit 1; } 

# test with GTF option not used
$scriptDir/make_metagenome.py "$testinputs/snps_subset.bed" "$testinputs/chr1_subset.fa" "$testoutdir/out" > "$testoutdir/out_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with make_metagenome.py, see $testoutdir/out_log.txt"; exit 1; }
# check all output files match
cmp --silent "$testoutdir/out.fa" "$exp/out.fa" || { echo "Error: output file $testoutdir/out.fa doesn't match expected output $exp/out.fa"; exit 1; }
cmp --silent "$testoutdir/out_metachrom.txt" "$exp/out_metachrom.txt" || { echo "Error: output file $testoutdir/out_metachrom.txt doesn't match expected output $exp/out_metachrom.txt"; exit 1; }

# test with GTF option enabled
$scriptDir/make_metagenome.py "$testinputs/snps_subset.bed" "$testinputs/chr1_subset.fa" "$testoutdir/out" --GTF "$testinputs/GTF_subset.gtf" > "$testoutdir/out_log2.txt"
[ $? != 0 ] && { echo "Error: something went wrong with make_metagenome.py, see $testoutdir/out_log2.txt"; exit 1; }
cmp --silent "$testoutdir/out_metagtf.gtf" "$exp/out_metagtf.gtf" || { echo "Error: output file $testoutdir/out_metagtf.gtf doesn't match expected output $exp/out_metagtf.gtf"; exit 1; }
cmp --silent "$testoutdir/out.fa" "$exp/out.fa" || { echo "Error: output file $testoutdir/out.fa doesn't match expected output $exp/out.fa"; exit 1; }
cmp --silent "$testoutdir/out_metachrom.txt" "$exp/out_metachrom.txt" || { echo "Error: output file $testoutdir/out_metachrom.txt doesn't match expected output $exp/out_metachrom.txt"; exit 1; }
echo "make_metagenome.py - all OK"

# test with SNPs to both A and B
$scriptDir/make_metagenome.py "$testinputs/snps_subset.bed" "$testinputs/chr1_subset.fa" "$testoutdir/out2" --ASNPs "$testinputs/snps_subset_w.bed" > "$testoutdir/out_log3.txt" 
[ $? != 0 ] && { echo "Error: something went wrong with make_metagenome.py, see $testoutdir/out_log3.txt"; exit 1; }
cmp --silent "$testoutdir/out2.fa" "$exp/out2.fa" || { echo "Error: output file $testoutdir/out2.fa doesn't match expected output $exp/out2.fa"; exit 1; }


# ---------------------
# Run tests on rna_seq_map
# ---------------------
echo ""
echo "Trying a few test runs of rna_seq_map.sh to make sure everything is working correctly..."
echo "Please be patient, this will take a minute or two"
testoutdir="$scriptDir/tests/rna_seq_map"; mkdir -p "$testoutdir"
exp="$expected/rna_seq_map"

[ -f "$testinputs/chr1_subset.fa" ] || { echo "Error: could not find FASTA genome $testinputs/chr1_subset.fa"; exit 1; } 
[ -f "$testinputs/GTF_subset.gtf" ] || { echo "Error: could not find GTF annotation file $testinputs/GTF_subset.gtf"; exit 1; } 
[ -f "$testinputs/subset_f.fq" ] || { echo "Error: could not find file $testinputs/subset_f.fq"; exit 1; } 
[ -f "$testinputs/subset_r.fq" ] || { echo "Error: could not find file $testinputs/subset_r.fq"; exit 1; } 
[ -f "$testinputs/subset_CxV_f.fq" ] || { echo "Error: could not find file $testinputs/subset_CxV_f.fq"; exit 1; } 
[ -f "$testinputs/subset_CxV_r.fq" ] || { echo "Error: could not find file $testinputs/subset_CxV_r.fq"; exit 1; } 
[ -f "$testinputs/subset_VxC_f.fq" ] || { echo "Error: could not find file $testinputs/subset_VxC_f.fq"; exit 1; } 
[ -f "$testinputs/subset_VxC_r.fq" ] || { echo "Error: could not find file $testinputs/subset_VxC_r.fq"; exit 1; } 

# make STAR-converted genomes using STAR --runMode genomeGenerate from both the regular FASTA genome
# and the metagenome generated above
mkdir -p "$scriptDir/tests/STAR_regular_genome"
mkdir -p "$scriptDir/tests/STAR_metagenome"

# regular genome
STAR --runMode genomeGenerate --outFileNamePrefix "$scriptDir/tests/STAR_regular_genome/build" --genomeDir "$scriptDir/tests/STAR_regular_genome" --genomeFastaFiles "$scriptDir/tests/test_inputs/chr1_subset.fa" --sjdbGTFfile "$scriptDir/tests/test_inputs/GTF_subset.gtf" --sjdbOverhang 30 --genomeSAindexNbases 6 > "$scriptDir/tests/STAR_regular_genome/build_log.txt" 2>&1
[ $? != 0 ] && { echo "Error: something went wrong with building the STAR index, see $scriptDir/tests/STAR_regular_genome/build_log.txt"; exit 1; }

# metagenome genome
STAR --runMode genomeGenerate --outFileNamePrefix "$scriptDir/tests/STAR_regular_genome/build" --genomeDir "$scriptDir/tests/STAR_metagenome" --genomeFastaFiles "$scriptDir/tests/make_metagenome/out.fa" --sjdbGTFfile "$scriptDir/tests/make_metagenome/out_metagtf.gtf" --sjdbOverhang 30 --genomeSAindexNbases 7 > "$scriptDir/tests/STAR_metagenome/build_log.txt" 2>&1
[ $? != 0 ] && { echo "Error: something went wrong with building the STAR index, see $scriptDir/tests/STAR_metagenome/build_log.txt"; exit 1; }

# single end, regular mode
$scriptDir/rna_seq_map.sh -1 "$testinputs/subset_f.fq" -g "$scriptDir/tests/STAR_regular_genome" -o "$testoutdir/test1" -3r > "$testoutdir/test1_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with rna_seq_map.sh, see $testoutdir/test1_log.txt"; exit 1; }

# paired end, regular mode
$scriptDir/rna_seq_map.sh -1 "$testinputs/subset_f.fq" -2 "$testinputs/subset_r.fq" -g "$scriptDir/tests/STAR_regular_genome" -o "$testoutdir/test2" -3r > "$testoutdir/test2_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with rna_seq_map.sh, see $testoutdir/test2_log.txt"; exit 1; }

# single end, metagenome mode
$scriptDir/rna_seq_map.sh -1 "$testinputs/subset_f.fq" -g "$scriptDir/tests/STAR_metagenome" -C "$scriptDir/tests/make_metagenome/out_metachrom.txt" -o "$testoutdir/test3" -3r > "$testoutdir/test3_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with rna_seq_map.sh, see $testoutdir/test3_log.txt"; exit 1; }

# paired end, metagenome mode
$scriptDir/rna_seq_map.sh -1 "$testinputs/subset_CxV_f.fq" -2 "$testinputs/subset_CxV_r.fq" -g "$scriptDir/tests/STAR_metagenome" -C "$scriptDir/tests/make_metagenome/out_metachrom.txt" -o "$testoutdir/test4" -3r > "$testoutdir/test4_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with rna_seq_map.sh, see $testoutdir/test4_log.txt"; exit 1; }

# paired end, metagenome mode
$scriptDir/rna_seq_map.sh -1 "$testinputs/subset_VxC_f.fq" -2 "$testinputs/subset_VxC_r.fq" -g "$scriptDir/tests/STAR_metagenome" -C "$scriptDir/tests/make_metagenome/out_metachrom.txt" -o "$testoutdir/test5" -3r > "$testoutdir/test5_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with rna_seq_map.sh, see $testoutdir/test5_log.txt"; exit 1; }

# check outputs match expected
for run in test1 test2 test3 test4 test5; do
	samtools view -ho - "$testoutdir/$run/STAR/expt_unique_alignments.bam" | awk '$1 !~ /^@/' > "$testoutdir/$run/STAR/expt_unique_alignments.sam"
	cmp --silent "$testoutdir/$run/STAR/expt_unique_alignments.sam" "$exp/$run/expt_unique_alignments.sam" || { echo "Error: output file $testoutdir/$run/STAR/expt_unique_alignments.sam doesn't match expected output $exp/$run/expt_unique_alignments.sam"; exit 1; }
	samtools view -ho - "$testoutdir/$run/STAR/expt_unique_alignments_dedup.bam" | awk '$1 !~ /^@/' > "$testoutdir/$run/STAR/expt_unique_alignments_dedup.sam"
	cmp --silent "$testoutdir/$run/STAR/expt_unique_alignments_dedup.sam" "$exp/$run/expt_unique_alignments_dedup.sam" || { echo "Error: output file $testoutdir/$run/STAR/expt_unique_alignments_dedup.sam doesn't match expected output $exp/$run/expt_unique_alignments_dedup.sam"; exit 1; }
done
echo "rna_seq_map.sh - all OK"

# ---------------------
# Run tests on call_imprinting
# ---------------------
echo ""
echo "Trying a few test runs of call_imprinting.sh to make sure everything is working correctly..."
echo "Please be patient, this will take a minute or two"
testoutdir="$scriptDir/tests/call_imprinting"; mkdir -p "$testoutdir"
exp="$expected/call_imprinting"

# will use outputs from the rna_seq_map runs above
AxB_bam="$scriptDir/tests/rna_seq_map/test4/STAR/expt_unique_alignments_dedup.bam"
BxA_bam="$scriptDir/tests/rna_seq_map/test5/STAR/expt_unique_alignments_dedup.bam"

AxB_A="$scriptDir/tests/test_inputs/MN47xKar_MN47_counts.txt"; [ -f "$AxB_A" ] || { echo "Error: couldn't find test input file $AxB_A"; exit 1; }
AxB_B="$scriptDir/tests/test_inputs/MN47xKar_Kar_counts.txt"; [ -f "$AxB_B" ] || { echo "Error: couldn't find test input file $AxB_B"; exit 1; }
BxA_A="$scriptDir/tests/test_inputs/KarxMN47_MN47_counts.txt"; [ -f "$BxA_A" ] || { echo "Error: couldn't find test input file $BxA_A"; exit 1; }
BxA_B="$scriptDir/tests/test_inputs/KarxMN47_Kar_counts.txt"; [ -f "$BxA_B" ] || { echo "Error: couldn't find test input file $BxA_B"; exit 1; }
filt="$scriptDir/tests/test_inputs/filter_out_test.txt"; [ -f "$filt" ] || { echo "Error: couldn't find test input file $filt"; exit 1; }

# starting from BAM
$scriptDir/call_imprinting.sh -R 2 -o "$testoutdir/test1" -1 "$AxB_bam" -2 "$BxA_bam" -S "$testinputs/snps_subset.bed" -G "$testinputs/GTF_subset.gtf" -r > "$testoutdir/test1_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with call_imprinting.sh, see $testoutdir/test1_log.txt"; exit 1; }

# starting from counts
$scriptDir/call_imprinting.sh -R 2 -o "$testoutdir/test2" -x "$AxB_A" -y "$AxB_B" -X "$BxA_A" -Y "$BxA_B" -A "MN47" -B "Kar" -M 85 -P 50 -f "$filt" -i 2 -r > "$testoutdir/test2_log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with call_imprinting.sh, see $testoutdir/test2_log.txt"; exit 1; }

# check outputs
cmp --silent "$testoutdir/test1/imprinting/expt_imprinting_filtered_all.txt" "$exp/test1/expt_imprinting_filtered_all.txt" || { echo "Error: output file $testoutdir/test1/imprinting/expt_imprinting_filtered_all.txt doesn't match expected output $exp/test1/expt_imprinting_filtered_all.txt"; exit 1; }
cmp --silent "$testoutdir/test2/imprinting/expt_imprinting_filtered_all.txt" "$exp/test2/expt_imprinting_filtered_all.txt" || { echo "Error: output file $testoutdir/test2/imprinting/expt_imprinting_filtered_all.txt doesn't match expected output $exp/test2/expt_imprinting_filtered_all.txt"; exit 1; }

echo "call_imprinting.sh - all OK"

echo ""
echo "All tests successful"

# Run tests on comp_imprinting.py
# ---------------------
echo "Trying a test run of comp_imprinting.py to make sure everything is working correctly..."
testoutdir="$scriptDir/tests/comp_imprinting"; mkdir -p "$testoutdir"
exp="$expected/comp_imprinting"

[ -f "$testinputs/species1_status.txt" ] || { echo "Error: could not find file $testinputs/species1_status.txt"; exit 1; } 
[ -f "$testinputs/species2_status.txt" ] || { echo "Error: could not find file $testinputs/species2_status.txt"; exit 1; } 
[ -f "$testinputs/homologs.txt" ] || { echo "Error: could not find file $testinputs/homologs.txt"; exit 1; } 
[ -f "$testinputs/pathway_info.txt" ] || { echo "Error: could not find file $testinputs/pathway_info.txt"; exit 1; } 

$scriptDir/comp_imprinting.py "$testinputs/species1_status.txt" "$testinputs/species2_status.txt" "$testinputs/homologs.txt" "$testoutdir/comp_imprinting" --pathway "$testinputs/pathway_info.txt" > "$testoutdir/log.txt"
[ $? != 0 ] && { echo "Error: something went wrong with comp_imprinting.py, see $testoutdir/log.txt"; exit 1; }

# check output file matches
cmp --silent "$testoutdir/comp_imprinting.txt" "$exp/comp_imprinting.txt" || { echo "Error: output file $testoutdir/comp_imprinting.txt doesn't match expected output $exp/comp_imprinting.txt"; exit 1; }
cmp --silent "$testoutdir/comp_imprinting_pathways.txt" "$exp/comp_imprinting_pathways.txt" || { echo "Error: output file $testoutdir/comp_imprinting_pathways.txt doesn't match expected output $exp/comp_imprinting_pathways.txt"; exit 1; }
echo "comp_imprinting.py - all OK"





























