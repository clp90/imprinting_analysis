# Imprinting Analysis v.1.0
This is a small suite of scripts designed to provide a consistent pipeline for imprinting analyses starting from raw RNA-seq data or count data. The ultimate goal of this project is to enable comparisons between imprinting calls made in different datasets or species. See [TBD] for the accompanying publication, which has a more in-depth description of the pipeline. The approach is based on the analysis in Gehring *et al.* 2011.

Imprinting is usually assessed using data from a pair of reciprocal crosses between two different strains or species. Throughout this tutorial, the two strains/species will be referred to as A and B, and the reciprocal crosses as AxB and BxA.

This pipeline should be suitable for assessing imprinting in any species, including plants and mammals. This pipeline is also suitable for tissues where the maternal : paternal dosage is not 1 (e.g. the triploid endosperm of most flowering plants). See the usage information for `call_imprinting.sh` for more details.
	
# Overview of pipeline
1. Convert genome into a "metagenome" using `make_metagenome.py` (optional), and index genome
2. Align sequencing reads from AxB and BxA to the (meta)genome using `rna_seq_map.sh`
3. Count allele-specific reads in both crosses and assess imprinting using `call_imprinting.sh`

Alternatively, `call_imprinting.sh` can be used directly with allele-specific count data (see below).
- `get_homologs.py` and `comp_imprinting.py` are also included in this repository to facilitate imprinting comparisons across species, where gene homology must be taken into account. 
	- `get_homologs.py` queries phytozome (Goodstein *et al.* 2011) to obtain homology information in a target species for a set of genes in a query species. `get_homologs.py` is plant-specific, and will not produce homolog information for species not in the phytozome database.
	- `comp_imprinting.py` performs a detailed comparison of imprinting between two different species, given imprinting data (ideally obtained from `call_imprinting.sh`, with both analyses using the same parameters).

### Required inputs if starting from raw RNA-seq data:
- two FASTQ files, one from AxB and one from BxA
```
@NB501288_421_HMLG5BGX7:1:11103:8457:7032#ACTCGCTATACTCCTT/1
GTAATATCCTCACAAAAAATTTCTTTCTAAGACATAGT
+NB501288_421_HMLG5BGX7:1:11103:8457:7032#ACTCGCTATACTCCTT/1
A6AAAEE/EEEEEEEEEEEEEAEAEEEAEEEEEEEEEE
@NB501288_421_HMLG5BGX7:1:11105:9575:9978#ACTCGCTATACTCCTT/1
CCTCATAGTTTTCTAGCACCTGAGACTCTGAAAAACTC
+NB501288_421_HMLG5BGX7:1:11105:9575:9978#ACTCGCTATACTCCTT/1
AAAAAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEE
```
- a BED file containing SNPs between AxB and BxA, where the SNP is given as [A allele]>[B allele]
	- note that BED files use 0-based indexing
```
Chr1	770	771	T>C
Chr1	964	965	C>T
Chr1	1112	1113	A>G
Chr1	1428	1429	C>A
```
- a FASTA file of the genome for read mapping
```
>Chr1
AGATTTTCACAACCACCACACAATTTATAACATTTAACAACTCATCATTTCAAGATAACAAGGAATTTAAACAGTGGGACATTATTTCAAGCTTGCAGTG
TATGAGAAATAACTGAAAAATATTTGTGGTCATCATAAATGAAATTTGTACATTTTAGCTCATTTAAGTTGTAGATAGACCACAAAAAGAAAACGGCTCA
```
- a GTF file of genome annotations
```
Chr1  Araport11  5UTR  959  1201  .  +  .  transcript_id "AT1G06530.1"; gene_id "AT1G06530"
Chr1  Araport11  exon  959  2379  .  +  .  transcript_id "AT1G06530.1"; gene_id "AT1G06530"
Chr1  Araport11  start_codon  1202  1204  .  +  .  transcript_id "AT1G06530.1"; gene_id "AT1G06530"
```

### Required inputs if starting from allele-specific count data:
- per-gene read counts from the A allele and B allele, in both AxB and BxA (4 files total)
```
AT1G04483	0
AT1G06530	24
AT1G06540	11
AT1G06550	4
```

# Installation
To install the scripts included in this project, simply download the repository, which includes the main scripts `make_metagenome.py`, `rna_seq_map.sh`, `call_imprinting.sh`, `get_homologs.py` and `comp_imprinting.py`, as well as several helper scripts in their own subfolder `/helper_scripts` and a testing script `run_tests.sh` that can verify the installation. 

Several other programs are required to run these scripts. These are all freely available, and are listed below along with installation instructions. Version indicated in square brackets was used in testing; if something is not working, consider installing the exact version in brackets below.

## Required additional software
All programs listed below must be on your `$PATH` - see installation instructions below for an example setup.

- `fastqc` (Andrews 2010) [0.11.8]
- `trim_galore` (Krueger 2012) [0.5.0]
- `samtools` (Li *et al.* 2009) [1.9]
- `STAR` (Dobin *et al.* 2012) [2.6.1d]
- `Java 1.8+` [1.8.0_191]
- `Python 2.7.11+` [2.7.15]
	- `HTSeq` (Anders *et al.* 2014) [0.11.0]
	- `biopython` (Cock *et al.* 2009) [1.72]
	- `cutadapt` (Martin 2011) [1.18]
	- `intermine` (Smith *et al.* 2012) [1.11.0]
- `R` [3.5.1]
	- `argparse` (Davis 2018) [1.1.1-1]
	- `Hmisc` (Harrell 2018) [4.1-1]

`MarkDuplicates.jar` from the picard-tools suite (Copyright (c) 2017 Broad Institute) is also required, but is provided in this repository (in the `helper_scripts` subdirectory; version 1.121 downloaded Oct. 2, 2015) for convenience. The picard-tools suite is licensed under the MIT license.

## Installation Instructions
Instructions were tested on Mac OS Sierra (v.10.12.6), in a bash shell.

### Basic Setup
To begin, open a new terminal window and make a new folder (for example, in `/Documents`) to store downloaded programs, and navigate to that directory. Let `installdir` be the name of that directory. For example, in terminal, type:
```
installdir="/path/to/where/programs/will/be/installed"
mkdir -p "$installdir"
cd "$installdir"
```
Replace `/path/to/where/programs/will/be/installed` with the actual path to the new folder you want to create. All instructions provided below assume that your current directory is the `$installdir`, and all programs (except `Java`, `python`, and `R`, which have their own installers) will be installed in this directory.

### Installing `Java`, `python` and `R`
Install the latest versions of `Java` (at least v.1.8), `python` and `R`. Installation of all three is generally straightforward - use the links provided below and follow the installation instructions for your setup (PC/Mac/etc.).
- `java` can be downloaded from https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html 
- `python` can be downloaded from https://www.python.org/downloads/release
	- python scripts in this repository require python 2.x (e.g. 2.7.15)
- `R` can be downloaded from https://cran.r-project.org/mirrors.html.
	- note that X11/XQuartz is required

### Installing required `python` modules
Python modules can be installed using `pip install`, which is included by default in python 2.7.11+. In your terminal, type:
```
pip install --upgrade pip
pip install --user --upgrade HTSeq
pip install --user --upgrade biopython
pip install --user --upgrade cutadapt
pip install --user --upgrade intermine
```
- **Note to Mac users**: you may need to run the `Install Certificates.command` script included in your python installation. It should be in your `Applications` folder in a folder called e.g. `Python 2.7`. Double click the file to run it; this will ensure that intermine can access the URL for phytozome in `get_homologs.py`. Otherwise you may get an error like `SSL error - certificate verify failed` when running the script.

From your terminal, run the following to make `cutadapt` visible to `trim_galore`:
```
iloc=$( pip show cutadapt 2>/dev/null | awk '$1~/Location:/{print $2}' )
nn=$( pip show cutadapt -f 2>/dev/null | awk '$1~/cutadapt$/{print $1}' )
ln -s "$iloc/$nn"
./cutadapt --help
```
This was successful if the last line produces usage information for `cutadapt`, instead of an error.

### Installing required `R` modules
The R modules `argparse` and `Hmisc` can be installed from inside `R`. In your terminal, type:
```
R
install.packages("devtools")
```
The first line launches R, the second installs `devtools`, a handy utility for installing some packages. You will need to choose a mirror for the installation before proceeding. Ensure the installation process completes, then install the two required modules by typing:
```
devtools::install_github("trevorld/argparse")
install.packages("Hmisc")
q()
```
The last line exits the R session.

### Installing `FASTQC`
Installation instructions and documentation are available at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Download the code from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip. Unzip the file and drop it into your `$installdir`. Assuming that the directory of the expanded zip archive is named `FastQC`, run the following in terminal:
```
mv FastQC fastqc_dir
ln -s fastqc_dir/fastqc
./fastqc --help
```
This should produce usage info for `FastQC` if it has been properly installed. The `ln -s` creates a symbolic link in your `$installdir` to the actual `fastqc` executable.

### Installing `trim_galore`
Documentation is available at http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. Download the zip file from https://github.com/FelixKrueger/TrimGalore/archive/0.5.0.zip, unzip it and copy the resulting directory to your `$installdir`. Assuming the expanded zip archive is named `TrimGalore-0.5.0`, run in terminal:
```
ln -s TrimGalore-0.5.0/trim_galore
./trim_galore --help
```
This should produce usage info for `trim_galore` if it has been properly installed.

### Installing `samtools`
Download the latest release from github (for example, from https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2), unzip the archive and copy the directory to your `$installdir`. Assuming the expanded zip archive is named `samtools-1.9`, run in terminal:
```
mkdir samtools_dir
cd samtools-1.9
./configure --prefix $(pwd)/../samtools_dir
make
make install
cd ..
rm -rf samtools-1.9 
ln -s samtools_dir/bin/samtools
./samtools
```
This should produce usage info for `samtools` if it has been properly installed.

### Installing `STAR`
Download the latest release from github (for example, https://github.com/alexdobin/STAR/archive/master.zip), unzip the archive and copy the directory to your `$installdir`. Assuming the expanded zip archive is named `STAR-master`, run in terminal:
- On Mac OSX:
```
ln -s STAR-master/bin/MacOSX_x86_64/STAR
./STAR
```
- On linux:
```
ln -s STAR-master/bin/Linux_x86_64/STAR
./STAR
```
This should produce usage info for `STAR` if it has been properly installed.

### Finishing Installation
The programs above must be on your `$PATH`  in order to be used by the scripts provided in this repository. Complete the installation by adding the `$installdir` to your `$PATH`:
```
export PATH="$PATH:$installdir"
```
**Note: this command will need to be repeated every time a new terminal is opened.**
By adding a symbolic link in the `$installdir` to the `fastqc`, `trim_galore`, `samtools` and `STAR` executables above, and then adding the `$installdir` to the `$PATH`, this effectively adds all four executables to the `$PATH` as well.

### Testing Installation
The script `run_tests.sh` included in this repository checks that all required dependencies are installed, and performs mini test runs to check that everything is working properly. Change directory to the location where this repository's folder was installed (which should contain `run_tests.sh`, `make_metagenome.py`, `rna_seq_map.sh`, `call_imprinting.sh` and `get_homology.py`), and run:

```
./run_tests.sh
```
If this run completes without any errors, all scripts should be ready for use. See below for some examples.

## Examples
### Starting from raw RNA-seq data:
Let `AxB.fq` and `BxA.fq` be a pair of raw sequencing files (FASTQ format), with genome `genome.fa` (FASTA format), annotations `annot.gtf` (GTF format), and SNP file `snps.bed` (BED format). Also let `$scriptDir` be the folder where this repository was installed (containing the scripts `run_tests.sh`, `make_metagenome.py`, etc.).

- (optional) Make a "metagenome" out of `genome.fa` using `snps.bed`:
```
$scriptDir/make_metagenome.py "snps.bed" "genome.fa" "metagenome" --GTF "annot.gtf"
```
- Index genome (TOP) or metagenome (BOTTOM) for use with STAR:
```
mkdir "STAR_index"
STAR --runMode genomeGenerate --outFileNamePrefix "STAR_genome_index" --genomeDir "STAR_genome_index" --genomeFastaFiles "genome.fa" --sjdbGTFfile "annot.gtf" --sjdbOverhang 30
```
```
mkdir "STAR_index"
STAR --runMode genomeGenerate --outFileNamePrefix "STAR_metagenome_index" --genomeDir "STAR_metagenome_index" --genomeFastaFiles "metagenome.fa" --sjdbGTFfile "metagenome_metagtf.gtf" --sjdbOverhang 30
```
- Filter and align reads to the STAR-indexed genome (TOP) or metagenome (BOTTOM):
```
$scriptDir/rna_seq_map.sh -1 "AxB.fq" -g "STAR_genome_index" -o "AxB_mapping_dir"
$scriptDir/rna_seq_map.sh -1 "BxA.fq" -g "STAR_genome_index" -o "BxA_mapping_dir"
```
```
$scriptDir/rna_seq_map.sh -1 "AxB.fq" -g "STAR_metagenome_index" -C "metagenome_metachrom.txt" -o "AxB_mapping_dir"
$scriptDir/rna_seq_map.sh -1 "BxA.fq" -g "STAR_metagenome_index" -C "metagenome_metachrom.txt" -o "BxA_mapping_dir"
```
- Get allele-specific counts and evaluate imprinting:
```
AxB_bam="AxB_mapping_dir/STAR/expt_unique_alignments_dedup.bam"
BxA_bam="BxA_mapping_dir/STAR/expt_unique_alignments_dedup.bam"
$scriptDir/call_imprinting.sh -o "imprinting_dir" -1 "$AxB_bam" -2 "$BxA_bam" -S "snps.bed" -G "annot.gtf" -A "A" -B "B"
```
- Note: we recommend mapping to a metagenome when possible to avoid mapping bias in favor of the strain with greater sequence similarity to the sequenced genome.

### Starting from count data:
Let `AxB_A.txt` and `AxB_B.txt` be per-gene allele-specific read counts from A and B respectively in the AxB cross, and `BxA_A.txt` and `BxA_B.txt` be per-gene allele-specific read counts from A and B respectively in the BxA cross. Again, let `$scriptDir` be the folder where this repository was installed.

- Evaluate imprinting using allele-specific counts:
```
$scriptDir/call_imprinting.sh -o "imprinting_dir" -x "AxB_A.txt" -y "AxB_B.txt" -X "BxA_A.txt" -Y "BxA_B.txt" -A "A" -B "B"
```

### Other examples:

Additionally, we provide two scripts, `get_homologs.py` and `comp_imprinting.py`, that can be useful for various imprinting analyses. Example commands for both scripts are given below.

#### Example `get_homologs.py` usage:
`get_homologs.py` takes a list of gene IDs, a query_species, and a target_species, and outputs the homologs of each gene ID from the query_species in the target_species. Given `genelist.txt`:
```
AT5G49160
AT1G69770
AT5G14620
AT2G32370
AT5G10140
```
This file contains *Arabidopsis thaliana* gene IDs; to obtain their *Arabidopsis lyrata* homologs, run:
```
$scriptDir/get_homologs.py "A. thaliana" "A. lyrata" "genelist.txt" "homolog_list.txt"
cat "homolog_list.txt"
```
which produces:
```
A. thaliana  A. lyrata
AT5G10140  AL6G20600
AT5G10140  AL6G20630
AT1G69770  AL2G28960
AT5G49160  AL8G22690
AT5G14620  AL6G25400
```
Note that species names must be in double quotes, in the form "[first letter of genus]. [species name]".

#### Example `comp_imprinting.py` usage:

`comp_imprinting.py` is designed to compare imprinting between two different species. Assuming data from both species has been processed according to the pipeline above to identify imprinted genes, this script takes:
- file containing the imprinting "status" and parental bias of each gene (columns 1, 7, and 25 of the [prefix]_filtered_all.txt file output by `call_imprinting.sh`) for species 1 and species 2 (2 files total):
```
AT1G02580	none	low_counts
AT1G15215	none	low_counts
AT1G63020	father	fail_pmat_cutoff
AT2G23380	mother	fail_pval_cutoff
AT2G27040	mother	fail_IF_cutoff
```
- A file containing species 1 -> species 2 homologs (e.g. `get_homologs.py` output):
```
A. thaliana	Z. mays
AT3G20740	GRMZM2G118205
AT3G20740	GRMZM2G148924
AT2G23380	GRMZM2G157820
```
The script will then assess the degree to which imprinted genes in species 1 are also imprinted in species 2, based on the homology information provided. The script will output a summary to `stdout` and a more detailed table in an output file with the specified prefix. An example command is given below:

```
$scriptDir/comp_imprinting.py species1_status.txt pecies2_status.txt homologs.txt comp_1vs2 --pathway pathway_info.txt
```

Note: optionally, users can also provide information on genes in pathways, protein complexes, or other groups of interest using the `--pathway` option. The script will then run an additional analysis to examine whether imprinting is conserved within these gene groups. User must provide a list of genes and the pathway/group for each gene, in the following format:
```
AT1G15215	RdDM
AT1G63020	RdDM
AT3G20740	PRC2
AT2G35670	PRC2
```
In this example, there are two groups, RdDM and PRC2, with two members each. Imprinting is considered conserved on the pathway level if either of those genes or their homologs are imprinted in both species.

## Additional Notes:
- Run `make_metagenome.py`, `rna_seq_map.sh`, `call_imprinting.sh` and `get_homology.py` without arguments to get usage info and additional information about each script.
- `rna_seq_map.sh` supports paired-end reads, just provide the forward reads file with `-1` and reverse reads with `-2`

## FAQ
- **I ran `rna_seq_map.sh` and few or no reads passed the quality filtering step (step 2), what's wrong?**

Be sure to check the quality encoding on your reads. `rna_seq_map.sh` assumes a PHRED+64 encoding by default; the user 			must use the `-3` option if they are encoded as PHRED+33 instead. The wikipedia page for the FASTQ file format (https://en.wikipedia.org/wiki/FASTQ_format) has a great explanation of the difference between the different quality encodings.
- **Running `get_homologs.py` produces an error like `SSL error - certificate verify failed`?**

See **Note to Mac users** above - this is a known bug/feature of the Python installer (affects both 2.x and 3.x from what I can tell) documented here: https://bugs.python.org/issue29480. The solution is to run the `Install Certificates.command` script included in your python installation. It should be in your `Applications` folder in a folder called e.g. `Python 2.7`. Double click the file to run it.

- **Running `get_homologs.py` produces a `WebserviceError` or `Internal Server Error`?**

This seems to be what happens when the phytozome server is down. The only real solution is to wait for the server to come back up.


## Citations
Anders S, Pyl PT, Huber W. HTSeq - A Python framework to work with high-throughput sequencing data. Bioinformatics 2014; doi: https://doi.org/10.1093/bioinformatics/btu638.

Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc.

Broad Institute. Picard Tools. 2018. Available online at: http://broadinstitute.github.io/picard/. Accessed: Oct. 2, 2015.

Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 2009; 25(11):1422-1423.

Davis TL. argparse: Command line optional and positional argument parser. 2018. R package version 1.1.1-1. Available online at: https://CRAN.R-project.org/package=argparse.

Dobin A, Davis CA, Schlesinger F, Drenkow J, *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 2012; 29(1):15-21.

Gehring M, Missirian V, Henikoff S. Genomic analysis of parent-of-origin allelic expression in *Arabidopsis thaliana* seeds. *PLoS ONE* 2011; doi: https://doi.org/10.1371/journal.pone.0023687.

Goodstein DM, Shu S, Howson R, *et al.* Phytozome: a comparative platform for green plant genomics. *Nucleic Acids Res.* 2011;40:D1178-86.

Harrell FE Jr, with contributions from Dupont C and many others. Hmisc: Harrell Miscellaneous. 2018. R package version 4.1-1. Available online at: https://CRAN.R-project.org/package=Hmisc.

Krueger F. Trim Galore. 2012. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/.

*Li H, Handsaker B*, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. *Bioinformatics* 2009; 25: 2078-9. *equal contribution*

Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 2011; 17(1):10-12. doi: https://doi.org/10.14806/ej.17.1.200.

Smith RN, Aleksic J, Butano D, Carr A, Contrino S, Hu F, Lyne M, Lyne R, Kalderimis A, Rutherford K, Stepan R, Sullivan J, Wakeling M, Watkins X, Micklem G. InterMine: a flexible data warehouse system for the integration and analysis of heterogeneous biological data. *Bioinformatics* 2012; 28(23):3163-5.
