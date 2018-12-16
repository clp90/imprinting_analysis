#!/usr/bin/env python

''' 
-------------------------
Usage: compare_imprinting.py status1.txt status2.txt homologs.txt outprefix [options]

This is a short script for comparing imprinting between two different species (species 1 and species 2). 
Uses output from call_imprinting.sh ("locus_name", "status" and "favored_parent" columns of *_filtered_all.txt file).

v.1.0	11/29/2018
by Colette Picard

Version history:
v.1.0 - initial build 11/29/2018

-------------------------
'''
 
import sys, os, re, argparse
plots = True
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
except ImportError:
	print "Warning: could not load module matplotlib. Plots will not be generated. To install matplotlib, run:"
	print "pip install matplotlib"
	plots = False

if len(sys.argv) == 1:
	print "-------------------------"
	print "compare_imprinting v1.0		by Colette L. Picard, 11/29/2018"
	print "-------------------------"
	print """
This is a short script for comparing imprinting between two different species (species 1 and species 2). 
Uses output from call_imprinting.sh ("locus_name", "status" and "favored_parent" columns of *_filtered_all.txt file).

All five arguments are required (and all files should NOT have headers)

status1.txt 	- list of imprinting "statuses" in species 1 (3 columns: gene ID, direction of parental bias, status; no header)
status2.txt		- list of imprinting "statuses" in species 2 (3 columns: gene ID, direction of parental bias, status; no header)
homolog.txt		- list of homologs (2 columns: species 1 gene ID, species 2 gene ID; no header)
outprefix 		- name for output file, without extension
type			- type of imprinting being compared, "MEG" or "PEG"

Example inputs:

(1) status1.txt (per-gene imprinting status in species 1)
AT1G01030	none	low_counts
AT1G01040	father	fail_pval_cutoff
AT1G01050	none	low_counts

(2) status2.txt (per-gene imprinting status in species 2)
AC148167.6_FG001	father	PEG
AC149475.2_FG002	father	fail_pval_cutoff
AC149475.2_FG003	father	fail_pval_cutoff

(3) homolog.txt (column order doesn't matter)
AT5G23110	GRMZM2G084819
AT4G35800	GRMZM2G044306
AT5G65930	GRMZM2G070273

Outputs, for each gene in implist, the "most imprinted" homolog in species 2 (if tied, chosen at random),
and the imprinting status of that homolog. Also outputs a summary to stdout.

*********************
NOTE: to also perform an analysis over known interactors, pathways, or complexes, provide 
optional file pathway_info.txt (5th input file):

Example pathway_info.txt for also doing a pathway/interaction/complex analysis:
AT3G43920	RdDM
AT2G27040	RdDM
AT2G33830	RdDM
AT3G20740	PRC2
AT2G35670	PRC2
AT1G02580	PRC2
GRMZM2G157820	PRC2
GRMZM5G875502	PRC2
GRMZM2G043484	PRC2
etc.

Gene ID from either species 1 or species 2 in first column, pathway ID in second column. All 
genes with same pathway ID (here "RdDM" or "PRC2") will be considered together. Note that
either species 1 or 2 ID can be provided, and scripts assumes that any homologs of that gene
in the other species -also- belong to the same pathway.

** Note: providing --pathway does not override the homologs-level analysis; both will be performed **

"""
	print "Usage: compare_imprinting.py status1.txt status2.txt homologs.txt outprefix [options]"
	print "-------------------------"
	print " --pathway : file containing groups of genes (e.g. belonging in same pathway) - 2 columns: gene ID, gene group"
	print " --species1 : name of species from first input file (default \"species1\""
	print " --species2 : name of species from second input file (default \"species2\""
	sys.exit(1)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('status1', help = 'list of imprinting "statuses" in species 1 (3 columns: gene ID, favored_parent, status)')
parser.add_argument('status2', help = 'list of imprinting "statuses" in species 2 (3 columns: gene ID, favored_parent, status)')
parser.add_argument('homolog', help = 'list of homologs (2 columns: species 1 gene ID, species 2 gene ID)')
parser.add_argument('outprefix', help = 'Prefix for output files {outprefix}.txt, {outprefix}_pie.png and {outprefix}_pie_all.png')
parser.add_argument('--pathway', default = "", help = 'Information about known pathways/complexes/etc. for some genes (2 columns: gene ID, pathway ID). See help for example.')
parser.add_argument('--species1', default = "species1", help = 'Name of species from first input file')
parser.add_argument('--species2', default = "species2", help = 'Name of species from second input file')

args = parser.parse_args()

status1 = args.status1
status2 = args.status2
homolog = args.homolog
outprefix = args.outprefix
if args.pathway != "":
	pathway = args.pathway
	use_pathway = True
else:
	use_pathway = False
species1 = args.species1
species2 = args.species2

#-------------------------------------------------------------
values = {"MEG":5, "PEG":5, "fail_pmat_cutoff":4, "fail_CEF_cutoff":3, "fail_IF_cutoff":2, "fail_pval_cutoff":1, "low_counts":0}
values_rev_MEG = {5:"MEG",4:"matbias_fail_pmat_cutoff", 3:"matbias_fail_CEF_cutoff", 2:"matbias_fail_IF_cutoff", 1:"fail_pval_cutoff", 0:"low_counts",-2:"patbias_fail_IF_cutoff",-3:"patbias_fail_CEF_cutoff",-4:"patbias_fail_pmat_cutoff",-5:"PEG"}
values_rev_PEG = {-5:"MEG",-4:"matbias_fail_pmat_cutoff", -3:"matbias_fail_CEF_cutoff", -2:"matbias_fail_IF_cutoff", 1:"fail_pval_cutoff", 0:"low_counts",2:"patbias_fail_IF_cutoff",3:"patbias_fail_CEF_cutoff",4:"patbias_fail_pmat_cutoff",5:"PEG"}
adj_MEG = {"mother":1, "father":-1, "none":1}
adj_PEG = {"mother":-1, "father":1, "none":1}

# check all inputs can be opened
try:
	st1 = open(status1, 'r') 
except IOError, e:
	print e
	print 'Could not open',species1,'status file',status1
	sys.exit(2)
	
try:
	st2 = open(status2, 'r') 
except IOError, e:
	print e
	print 'Could not open',species2,'status file',status2
	sys.exit(2)

try:
	hom = open(homolog, 'r') 
except IOError, e:
	print e
	print 'Could not open homolog list',homolog
	sys.exit(2)

try:
	out = open(outprefix+'.txt', 'w') 
except IOError, e:
	print e
	print 'Could not create output file ',outprefix+'.txt'
	sys.exit(2)
	
if use_pathway == True:
	try:
		pwy = open(pathway, 'r') 
	except IOError, e:
		print e
		print 'Could not open pathway file ',pathway
		sys.exit(2)

	try:
		out2 = open(outprefix+'_pathways.txt', 'w') 
	except IOError, e:
		print e
		print 'Could not create output file ',outprefix+'_pathways.txt'
		sys.exit(2)
	

# read in list of genes assayed in species 1 and their statuses
s1_status_MEG = {}			# store ranking relative to MEG (higher value is closer to being MEG)
s1_status_PEG = {}			# store ranking relative to PEG (higher value is closer to being PEG)
# values 1->5 correspond to fail_pval -> imprinted, -2->-5 = fail_IF -> imprinted in other direction
# -6 == censored
line = st1.readline()
while line:
	ll = line.strip().split('\t')
	if len(ll) != 3:
		print "Error: three columns expected in",species1,"status file,",len(ll),"detected"
		sys.exit(1)
	if ll[0] in s1_status_MEG:
		print "Error: gene",ll[0],"appears more than once in",species1,"status file"
		sys.exit(1)
	if ll[0] in s1_status_PEG:
		print "Error: gene",ll[0],"appears more than once in",species1,"status file"
		sys.exit(1)
	else:
		# convert each gene's status to numeric ranking
		if ll[2] == "censored":
			s1_status_MEG[ll[0]] = -6
			s1_status_PEG[ll[0]] = -6
		else:
			# ranking vs. MEGs
			adj = adj_MEG[ll[1]]
			val = values[ll[2]] * adj
			if val == -1:
				val = 1		
			s1_status_MEG[ll[0]] = val
		
			# ranking vs. PEGs
			adj = adj_PEG[ll[1]]
			val = values[ll[2]] * adj
			if val == -1:
				val = 1		
			s1_status_PEG[ll[0]] = val
	line = st1.readline()
st1.close()

genelist_s1 = set(s1_status_MEG.keys()+s1_status_PEG.keys())


# repeat in species 2
s2_status_MEG = {}			# store ranking relative to MEG (higher value is closer to being MEG)
s2_status_PEG = {}			# store ranking relative to PEG (higher value is closer to being PEG)
line = st2.readline()
while line:
	ll = line.strip().split('\t')
	if len(ll) != 3:
		print "Error: three columns expected in",species2,"status file,",len(ll),"detected"
		sys.exit(1)
	if ll[0] in s2_status_MEG:
		print "Error: gene",ll[0],"appears more than once in",species2,"status file"
		sys.exit(1)
	if ll[0] in s2_status_PEG:
		print "Error: gene",ll[0],"appears more than once in",species2,"status file"
		sys.exit(1)
	else:
		# convert each gene's status to numeric ranking
		if ll[2] == "censored":
			s2_status_MEG[ll[0]] = -6
			s2_status_PEG[ll[0]] = -6
		else:
			# ranking vs. MEGs
			adj = adj_MEG[ll[1]]
			val = values[ll[2]] * adj
			if val == -1:
				val = 1		
			s2_status_MEG[ll[0]] = val
		
			# ranking vs. PEGs
			adj = adj_PEG[ll[1]]
			val = values[ll[2]] * adj
			if val == -1:
				val = 1		
			s2_status_PEG[ll[0]] = val
	line = st2.readline()
st2.close()


# read in list of species 1 -> 2 homologs; detect pathway data if included in col 3
line = hom.readline()
s1_to_s2 = {}
s2_to_s1 = {}
	
while line:
	ll = line.strip().split('\t')
	if len(ll) != 2:
		print "Error: two columns expected in list of homologs,",len(ll),"detected"
		print line
		sys.exit(1)
	
	# read in so that order of genelist file (species 1>2 or 2>1) doesn't matter
	if ll[0] in genelist_s1:
		if ll[0] not in s1_to_s2:
			s1_to_s2[ll[0]] = [ll[1]]
		else:
			s1_to_s2[ll[0]].append(ll[1])			
		if ll[1] not in s2_to_s1:
			s2_to_s1[ll[1]] = [ll[0]]
		else:
			s2_to_s1[ll[1]].append(ll[0])			
	elif ll[1] in genelist_s1:
		if ll[1] not in s1_to_s2:
			s1_to_s2[ll[1]] = [ll[0]]
		else:
			s1_to_s2[ll[1]].append(ll[0])			
		if ll[0] not in s2_to_s1:
			s2_to_s1[ll[0]] = [ll[1]]
		else:
			s2_to_s1[ll[0]].append(ll[1])			
	line = hom.readline()
hom.close()


# if pathway data provided, read that in
if use_pathway == True:
	line = pwy.readline()
	pathways = {}
	dropped = 0
	
	while line:
		ll = line.strip().split('\t')
		if len(ll) != 2:
			print "Error: two columns expected in pathway data,",len(ll),"detected"
			print line
			sys.exit(1)
			
		if ll[1] not in pathways:
			pathways[ll[1]] = {"species1":[],"species2":[]}	
						
		if ll[0] in s1_status_MEG:
			if ll[0] not in pathways[ll[1]]["species1"]:
				pathways[ll[1]]["species1"].append(ll[0])
			if ll[0] in s1_to_s2:
				for s2_gene in s1_to_s2[ll[0]]:
					if s2_gene not in pathways[ll[1]]["species2"]:
						pathways[ll[1]]["species2"].append(s2_gene)
		elif ll[0] in s2_status_MEG:
			if ll[0] not in pathways[ll[1]]["species2"]:
				pathways[ll[1]]["species2"].append(ll[0])
			if ll[0] in s2_to_s1:
				for s1_gene in s2_to_s1[ll[0]]:
					if s1_gene not in pathways[ll[1]]["species1"]:
						pathways[ll[1]]["species1"].append(s1_gene)
		else:
			dropped += 1
		line = pwy.readline()
	
	if dropped != 0:		
		print "Pathway information was provided for",dropped,"genes that didn't have imprinting data; these were censored"
	
	pwy.close()


# put it all together; first on the homolog level
print "Analysing conservation of imprinted expression between",species1,"and",species2,"homologs"
out.write('geneID_species1\tstatus_species1\tbest_homolog_species2\tbest_status_species2\tdiff\n')

# save summaries of statuses of the s2 homologs of s1 MEGs and PEGs for summary
counts_s1_MEG = {"no_homolog":0, "censored":0, "no_data":0, "MEG":0,"matbias_fail_pmat_cutoff":0, "matbias_fail_CEF_cutoff":0, "matbias_fail_IF_cutoff":0, "fail_pval_cutoff":0, "low_counts":0,"patbias_fail_IF_cutoff":0,"patbias_fail_CEF_cutoff":0,"patbias_fail_pmat_cutoff":0,"PEG":0}
counts_s1_PEG = {"no_homolog":0, "censored":0, "no_data":0, "MEG":0,"matbias_fail_pmat_cutoff":0, "matbias_fail_CEF_cutoff":0, "matbias_fail_IF_cutoff":0, "fail_pval_cutoff":0, "low_counts":0,"patbias_fail_IF_cutoff":0,"patbias_fail_CEF_cutoff":0,"patbias_fail_pmat_cutoff":0,"PEG":0}

values_rev_MEG = {-7:"no_data", -6:"censored", 5:"MEG",4:"matbias_fail_pmat_cutoff", 3:"matbias_fail_CEF_cutoff", 2:"matbias_fail_IF_cutoff", 1:"fail_pval_cutoff", 0:"low_counts",-2:"patbias_fail_IF_cutoff",-3:"patbias_fail_CEF_cutoff",-4:"patbias_fail_pmat_cutoff",-5:"PEG"}
values_rev_PEG = {-7:"no_data", -6:"censored", -5:"MEG",-4:"matbias_fail_pmat_cutoff", -3:"matbias_fail_CEF_cutoff", -2:"matbias_fail_IF_cutoff", 1:"fail_pval_cutoff", 0:"low_counts",2:"patbias_fail_IF_cutoff",3:"patbias_fail_CEF_cutoff",4:"patbias_fail_pmat_cutoff",5:"PEG"}

x = 0
for s1 in genelist_s1:
	if s1 not in s1_to_s2:
		best_str = "no_homolog"
		best_homolog = "none"
		diff = ""

		if s1_status_MEG[s1] == 5:
			counts_s1_MEG[best_str]+=1
		if s1_status_PEG[s1] == 5:
			counts_s1_PEG[best_str]+=1	
	else:
		best = -7					# default, "no data"
		best_homolog = "none"	
		diff = ""
				
		# gene is maternally biased in species 1, find most maternally-biased homolog
		if s1_status_MEG[s1] >= 2:
			for s2 in s1_to_s2[s1]:
				if s2 in s2_status_MEG and s2_status_MEG[s2] > best:
					best = s2_status_MEG[s2]
					best_homolog = s2
			best_str = values_rev_MEG[best]
			if best_homolog != "none" and abs(s2_status_MEG[best_homolog]) <= 5 and s2_status_MEG[best_homolog] != 0:
				diff = s1_status_MEG[s1] - s2_status_MEG[best_homolog]
			
			if s1_status_MEG[s1] == 5:
				counts_s1_MEG[best_str]+=1

		# gene is paternally biased in species 1, find most paternally-biased homolog
		elif s1_status_MEG[s1] <= -2 and s1_status_MEG[s1] >= -5:
			for s2 in s1_to_s2[s1]:
				if s2 in s2_status_PEG and s2_status_PEG[s2] > best:
					best = s2_status_PEG[s2]
					best_homolog = s2
			best_str = values_rev_PEG[best]
			if best_homolog != "none" and abs(s2_status_PEG[best_homolog]) <= 5 and s2_status_PEG[best_homolog] != 0:
				diff = s1_status_PEG[s1] - s2_status_PEG[best_homolog]

			if s1_status_PEG[s1] == 5:
				counts_s1_PEG[best_str]+=1
	
		# gene is not parentally biased in species 1; find most parentally-biased homolog
		elif s1_status_MEG[s1] >= -5:
			bestpar = "mat"
			for s2 in s1_to_s2[s1]:
				if s2 in s2_status_MEG and s2_status_MEG[s2] > best:
					best = s2_status_MEG[s2]
					best_homolog = s2
					bestpar = "mat"
				if s2 in s2_status_PEG and s2_status_PEG[s2] > best:
					best = s2_status_PEG[s2]
					best_homolog = s2
					bestpar = "pat"
			if bestpar == "mat":
				best_str = values_rev_MEG[best]
				if best_homolog != "none" and abs(s2_status_MEG[best_homolog]) <= 5 and s2_status_MEG[best_homolog] != 0:
					diff = s1_status_MEG[s1] - s2_status_MEG[best_homolog]
			else:
				best_str = values_rev_PEG[best]
				if abs(s2_status_PEG[best_homolog]) <= 5 and s2_status_MEG[best_homolog] != 0:
					diff = s1_status_MEG[s1] - s2_status_PEG[best_homolog]
					
		# else, gene was censored or had no data in species 1			
	
	if s1_status_MEG[s1] > -2:
		cur_str = values_rev_MEG[s1_status_MEG[s1]]
	else:
		cur_str = values_rev_PEG[s1_status_PEG[s1]]
	
	out.write(s1+'\t'+cur_str+'\t'+best_homolog+'\t'+best_str+'\t'+str(diff)+'\n')


# output summary of results of homologs-level analysis; MEGs and PEGs in species 1 only
print "Done. Results saved to",outprefix+".txt"
cts = [counts_s1_MEG, counts_s1_PEG]
types = ["MEGs","PEGs"]
types2 = ["MEG","PEG"]
bias = ["maternally","paternally"]
prefix = ["matbias_","patbias_"]
print ""
print "Summary of results for MEGs and PEGs in",species1+":"
print "----------------------"

for i in range(0,2):
	j = abs(i-1)
	tot = sum(cts[i].values())
	print "A total of",tot,types[i],"were detected in",species1
	if tot > 0:
		print " -",cts[i]["no_homolog"],"("+str(round(cts[i]["no_homolog"]/float(tot)*100,1))+"%)",types[i],"had no homolog in",species2
		print " -",cts[i]["no_data"],"("+str(round(cts[i]["no_data"]/float(tot)*100,1))+"%)",types[i],"had homolog in",species2+", but no data was available for any homolog"
		print " -",cts[i]["low_counts"],"("+str(round(cts[i]["low_counts"]/float(tot)*100,1))+"%)",types[i],"had homolog in",species2+", but too few allele-specific reads were available for any homolog to evaluate"
		print " -",cts[i]["censored"],"("+str(round(cts[i]["censored"]/float(tot)*100,1))+"%)",types[i],"had homolog in",species2+", but all homolog(s) with data were censored in species 2 analysis"

		could_eval = tot - cts[i]["no_homolog"] - cts[i]["no_data"] - cts[i]["low_counts"] - cts[i]["censored"]
		print ""
		print "The remaining",could_eval,types[i],"("+str(round(could_eval/float(tot)*100,1))+"%)","had at least one homolog that could be evaluated for imprinting in",species2

		if could_eval > 0:
			print "Of these:"
			print " -",cts[i][types2[i]],"("+str(round(cts[i][types2[i]]/float(could_eval)*100,1))+"%) were also",types[i],"in",species2
			print " -",cts[i][prefix[i]+"fail_pmat_cutoff"],"("+str(round(cts[i][prefix[i]+"fail_pmat_cutoff"]/float(could_eval)*100,1))+"%) were substantially",bias[i],"biased in",species2,"but failed the % maternal cutoff (failed pmat cutoff)"
			print " -",cts[i][prefix[i]+"fail_CEF_cutoff"],"("+str(round(cts[i][prefix[i]+"fail_CEF_cutoff"]/float(could_eval)*100,1))+"%) were substantially",bias[i],"biased in",species2,"but also exhibited significant strain bias (failed CEF cutoff)"
			print " -",cts[i][prefix[i]+"fail_IF_cutoff"],"("+str(round(cts[i][prefix[i]+"fail_IF_cutoff"]/float(could_eval)*100,1))+"%) were significantly but not substantially",bias[i],"biased in",species2,"(failed IF cutoff)"
			print ""
			print " -",cts[i]["fail_pval_cutoff"],"("+str(round(cts[i]["fail_pval_cutoff"]/float(could_eval)*100,1))+"%) were not significantly parentally biased in",species2,"(failed pval cutoff)"
			print ""
			print " -",cts[i][prefix[j]+"fail_IF_cutoff"],"("+str(round(cts[i][prefix[j]+"fail_IF_cutoff"]/float(could_eval)*100,1))+"%) were significantly but not substantially",bias[j],"biased in",species2,"(failed IF cutoff)"
			print " -",cts[i][prefix[j]+"fail_CEF_cutoff"],"("+str(round(cts[i][prefix[j]+"fail_CEF_cutoff"]/float(could_eval)*100,1))+"%) were substantially",bias[j],"biased in",species2,"but also exhibited significant strain bias (failed CEF cutoff)"
			print " -",cts[i][prefix[j]+"fail_pmat_cutoff"],"("+str(round(cts[i][prefix[j]+"fail_pmat_cutoff"]/float(could_eval)*100,1))+"%) were substantially",bias[j],"biased in",species2,"but failed the % maternal cutoff (failed pmat cutoff)"
			print " -",cts[i][types2[j]],"("+str(round(cts[i][types2[j]]/float(could_eval)*100,1))+"%) were also",types[j]+"s in",species2

	print "----------------------"
	print ""

# make pie-chart summaries
if plots == True and tot > 0 and could_eval > 0:
	print ""
	print "Making summary pie charts"
	
	for i in range(0,2):
		y1 = []

		if types[i] == "MEGs":
			x1 = ["no homolog", "no data", "low counts", "censored", "conserved MEG", "matbias fail pmat", "matbias fail CEF", "matbias fail IF", "fail pval", "patbias fail IF", "patbias fail CEF", "patbias fail pmat", "PEG" ]
			for vv in ["no_homolog", "no_data", "low_counts", "censored", "MEG", "matbias_fail_pmat_cutoff", "matbias_fail_CEF_cutoff", "matbias_fail_IF_cutoff", "fail_pval_cutoff", "patbias_fail_IF_cutoff", "patbias_fail_CEF_cutoff", "patbias_fail_pmat_cutoff", "PEG"]:
				y1.append(counts_s1_MEG[vv])
			c1 = ["#DF8344", "#D9D9D9", "#A6A6A6", "#7EAB55", "#B02418", "#D28589", "#F4B8B9", "#F9DBDF", "#FDF2D0", "#E0EBF6", "#A4C2E3", "#6A99D0", "#4E71BE" ]
		else:
			x1 = ["no homolog", "no data", "low counts", "censored", "conserved PEG", "patbias fail pmat", "patbias fail CEF", "patbias fail IF", "fail pval", "matbias fail IF", "matbias fail CEF", "matbias fail pmat", "MEG" ]
			for vv in ["no_homolog", "no_data", "low_counts", "censored", "PEG", "patbias_fail_pmat_cutoff", "patbias_fail_CEF_cutoff", "patbias_fail_IF_cutoff", "fail_pval_cutoff", "matbias_fail_IF_cutoff", "matbias_fail_CEF_cutoff", "matbias_fail_pmat_cutoff", "MEG"]:
				y1.append(counts_s1_PEG[vv])
			c1 = ["#DF8344", "#D9D9D9", "#A6A6A6", "#7EAB55", "#4E71BE", "#6A99D0", "#A4C2E3", "#E0EBF6", "#FDF2D0", "#F9DBDF", "#F4B8B9", "#D28589", "#B02418" ]
			
		x2 = x1[4:]		# without no_homolog, no_data and low_counts (e.g. could be evaluated in both)
		y2 = y1[4:]
		c2 = c1[4:]
			
		def my_autopct(pct):
			return ('%.1f%%' % pct) if pct > 2 else ''

		fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))
		wedges, texts, autotexts = ax.pie(y1, autopct=my_autopct, textprops=dict(color="black"), colors=c1)
		ax.legend(wedges, x1, title="Status", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
		plt.savefig(outprefix+"_pie_"+types[i]+"_all.png", fmt='png', dpi=200)

		fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))
		wedges, texts, autotexts = ax.pie(y1, autopct='', textprops=dict(color="black"), colors=c1)
		ax.legend(wedges, x1, title="Status", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
		plt.savefig(outprefix+"_pie_"+types[i]+"_all_nolbl.png", fmt='png', dpi=200)

		fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))
		wedges, texts, autotexts = ax.pie(y2, autopct=my_autopct, textprops=dict(color="black"), colors=c2)
		ax.legend(wedges, x2, title="Status", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
		plt.savefig(outprefix+"_pie_"+types[i]+".png", fmt='png', dpi=200)

		fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))
		wedges, texts, autotexts = ax.pie(y2, autopct='', textprops=dict(color="black"), colors=c2)
		ax.legend(wedges, x2, title="Status", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
		plt.savefig(outprefix+"_pie_"+types[i]+"_nolbl.png", fmt='png', dpi=200)

	print "Done"
	
# do pathway analysis if --pathway file provided
if use_pathway == True:
	print "Analysing imprinted expression across genes in provided pathways. Summary below:"
	print "----------------------"
	print ""
	out2.write('pathway_ID\tspecies\tgeneID\tspecies2_homologs\tstatus\n')
	
	print "pathway_ID	species_1_bias	species_1_most_biased	species_2_bias	species_2_most_biased"

	# for each pathway, get favored direction in s1 and s2 (strongest imprinted gene)
	for pp in pathways:
	
		s1_best_mat_str = "no_data"
		s1_best_mat = ""
		s1_best_mat_val = -7

		s1_best_pat_str = "no_data"
		s1_best_pat = ""
		s1_best_pat_val = -7
		
		for s1_gene in pathways[pp]["species1"]:
			if s1_gene in s1_status_MEG:
#				print s1_gene,s1_status_MEG[s1_gene],s1_status_PEG[s1_gene]
				# grab s2 homolog if it exists
				hlog = ""
				if s1_gene in s1_to_s2:
					for s2_gene in s1_to_s2[s1_gene]:
						hlog = hlog + ',' + s2_gene
					hlog = hlog[1:]
				if s1_status_MEG[s1_gene] >= s1_status_PEG[s1_gene]:
					out2.write(pp+'\t'+species1+'\t'+s1_gene+'\t'+hlog+'\t'+values_rev_MEG[s1_status_MEG[s1_gene]]+'\n')		
				else:
					out2.write(pp+'\t'+species1+'\t'+s1_gene+'\t'+hlog+'\t'+values_rev_PEG[s1_status_PEG[s1_gene]]+'\n')		
						
				if s1_status_MEG[s1_gene] != -6 and s1_status_MEG[s1_gene] != 0:	# omits censored s1 genes, genes with no data
#					print s1_status_MEG[s1_gene],s1_status_PEG[s1_gene],s1_best_mat_val
					if s1_status_MEG[s1_gene] >= s1_status_PEG[s1_gene] and s1_status_MEG[s1_gene] > s1_best_mat_val:
						# mat bias or no bias
						s1_best_mat_val = s1_status_MEG[s1_gene]
						s1_best_mat = s1_gene
#						print "update mat, now:",s1_best_mat_val
					elif s1_status_PEG[s1_gene] > s1_best_pat_val:
						# pat bias
						s1_best_pat_val = s1_status_PEG[s1_gene]
						s1_best_pat = s1_gene
#						print "update pat, now:",s1_best_pat_val
											
		if s1_best_mat_val == -7 and s1_best_pat_val == -7:
			# no informative genes in s1
			s1_best_val = -7
			s1_best = ""
			s1_best_str = "no_informative_genes"		
		elif s1_best_mat_val > s1_best_pat_val:
			# overall maternal bias
			s1_best_val = s1_best_mat_val
			s1_best = s1_best_mat
			s1_best_str = values_rev_MEG[s1_best_val]
		elif s1_best_mat_val == s1_best_pat_val:
			# overall equal bias in both directions (this should be rare)
			s1_best_val = s1_best_mat_val
			s1_best = s1_best_mat
			s1_best_str = "both_"+values_rev_MEG[s1_best_val][8:]
		else:
			# overall paternal bias
			s1_best_val = s1_best_pat_val
			s1_best = s1_best_pat
			s1_best_str = values_rev_PEG[s1_best_val]
			
		if s1_best_val == 1:
			s1_best_str = "no_parental_bias"
			
#		print pp,s1_best_val,s1_best_str,s1_best
											
		# repeat for all the s2 homologs of the s2 genes
		s2_best_mat_str = "no_data"
		s2_best_mat = ""
		s2_best_mat_val = -7

		s2_best_pat_str = "no_data"
		s2_best_pat = ""
		s2_best_pat_val = -7
		
		for s2_gene in pathways[pp]["species2"]:
			if s2_gene in s2_status_MEG:
				hlog = ""
				if s2_gene in s2_to_s1:
					for s1_gene in s2_to_s1[s2_gene]:
						hlog = hlog + ',' + s1_gene
					hlog = hlog[1:]
#				print s2_gene,s2_status_MEG[s2_gene],s2_status_PEG[s2_gene]
				if s2_status_MEG[s2_gene] >= s2_status_PEG[s2_gene]:
					out2.write(pp+'\t'+species2+'\t'+s2_gene+'\t'+hlog+'\t'+values_rev_MEG[s2_status_MEG[s2_gene]]+'\n')		
				else:
					out2.write(pp+'\t'+species2+'\t'+s2_gene+'\t'+hlog+'\t'+values_rev_PEG[s2_status_PEG[s2_gene]]+'\n')		

				if s2_status_MEG[s2_gene] != -6 and s2_status_MEG[s2_gene] != 0:	# omits censored s2 genes, genes with no data

					if s2_status_MEG[s2_gene] >= s2_status_PEG[s2_gene] and s2_status_MEG[s2_gene] > s2_best_mat_val:
						# mat bias or no bias
						s2_best_mat_val = s2_status_MEG[s2_gene]
						s2_best_mat = s2_gene
					elif s2_status_PEG[s2_gene] > s2_best_pat_val:
						# pat bias
						s2_best_pat_val = s2_status_PEG[s2_gene]
						s2_best_pat = s2_gene

		if s2_best_mat_val == -7 and s2_best_pat_val == -7:
			# no informative genes in s2
			s2_best_val = -7
			s2_best = ""
			s2_best_str = "no_informative_genes"		
		elif s2_best_mat_val > s2_best_pat_val:
			# overall maternal bias
			s2_best_val = s2_best_mat_val
			s2_best = s2_best_mat
			s2_best_str = values_rev_MEG[s2_best_val]
		elif s2_best_mat_val == s2_best_pat_val:
			# overall equal bias in both directions (this should be rare)
			s2_best_val = s2_best_mat_val
			s2_best = s2_best_mat
			if s2_best_mat_val < 5:
				s2_best_str = "both_"+values_rev_MEG[s2_best_val][8:]
			else:
				s2_best_str = "both_MEGs_and_PEGs"
		else:
			# overall paternal bias
			s2_best_val = s2_best_pat_val
			s2_best = s2_best_pat
			s2_best_str = values_rev_PEG[s2_best_val]
		
		if s2_best_val == 1:
			s2_best_str = "no_parental_bias"

#		print pp,s2_best_val,s2_best_str,s2_best

		print pp+'\t'+s1_best_str+'\t'+s1_best+'\t'+s2_best_str+'\t'+s2_best
		
	print ""
	print "Done. Full results saved to",outprefix+"_pathways.txt"





