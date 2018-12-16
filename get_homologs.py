#!/usr/bin/env python

''' 
-------------------------
Usage: get_homologs.py query_species target_species genelist.txt outfile.txt [options]

Example: get_homologs.py "A. thaliana" "A. lyrata" mygenelist.txt myoutfile.txt

This script queries phytozome to map one set of gene IDs from query_species to
their homologs in target_species. 

Requires python module intermine (Smith et al. 2012; also see https://doi.org/10.14806/ej.17.1.200)

v.1.0	11/17/2018
by Colette Picard

Version history:
v.1.0 - initial build 11/17/2018

-------------------------
'''
 
import sys, os, re, argparse
try:
	from intermine.webservice import Service
except ImportError:
	print "Error: could not load module intermine. To install intermine, run:"
	print "pip install intermine"
	sys.exit(1)

if len(sys.argv) == 1:
	print "-------------------------"
	print "get_homologs v1.0		by Colette L. Picard, 11/17/2018"
	print "-------------------------"
	print """This is a simple script to query phytozome for homology data. User specifies
the query and target species, and provides a list of gene IDs (no header) from the query
species. This script obtains the target_species homolog of all the query_species gene IDs
provided in genelist.txt and outputs them. For example, given genelist.txt:

AT5G49160
AT1G69770
AT5G14620
AT2G32370
AT5G10140

running get_homologs.py "A. thaliana" "A. lyrata" genelist.txt testout.txt

produces the output:
A. thaliana	A. lyrata
AT5G10140	AL6G20600
AT5G10140	AL6G20630
AT1G69770	AL2G28960
AT5G49160	AL8G22690
AT5G14620	AL6G25400

NOTE: query_species genes without a homolog will not appear in the output file.
A gene in query_species may have more than one reported target_species homolog.

NOTE2: species name should be in double quotes, with first initial of genus and
full name of species (e.g. "A. thaliana")
"""
	print "Usage: get_homologs.py speciesA speciesB genelist.txt outfile.txt [options]"
	print "-------------------------"
	sys.exit(1)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('query_species', help = 'name of query species - gene IDs in genelist must be from this species')
parser.add_argument('target_species', help = 'name of target species - will obtain homologs of gene IDs in genelist in this species')
parser.add_argument('genelist', help = 'List of gene IDs, without header')
parser.add_argument('outfile', help = 'Name for output file')

args = parser.parse_args()

query_species = args.query_species
target_species = args.target_species
genelist = args.genelist
outfile = args.outfile

print "Running get_homologs v1.0		by Colette L. Picard, 11/17/2018"
print "-------------------------"
print "Query species:",query_species
print "Target species:",target_species
print "List of gene IDs to convert:",genelist
print "Saving result to output file:",outfile
print "-------------------------"

#-------------------------------------------------------------

# try to connect to phytozome
try:
	service = Service("https://phytozome.jgi.doe.gov/phytomine/service")
except Exception as e:
	print e
	print "Error connecting to phytomine"
	sys.exit(1)
		
# try to open list of gene IDs
num_homologs = {}
try:
	f = open(genelist, 'r') 
except IOError, e:
	print e
	print 'Could not open gene list',genelist,'...'
	sys.exit(2)
	
# try to open outfile
try:
	o = open(outfile, 'w') 
except IOError, e:
	print e
	print 'Could not create output file',outfile,'...'
	sys.exit(2)

o.write(query_species+'\t'+target_species+'\n')

glist = []
line = f.readline()
while line:	
	ll = line.strip().split()
	if len(ll) > 1:
		print "Error: input file has at least two columns; only one column (gene IDs) expected"
		sys.exit(1)
	glist.append(ll[0])
	line = f.readline()

# query database for all genes in the list
query = service.new_query("Homolog")
query.add_view("ortholog_gene.primaryIdentifier")
query.add_constraint("gene.primaryIdentifier", "ONE OF", glist, code = "A")
query.add_constraint("gene.organism.shortName", "=", query_species, code = "B")
query.add_constraint("ortholog_gene.organism.shortName", "=", target_species, code = "C")
query.set_logic("A and B and C")
query.add_view(
    "ortholog_gene.primaryIdentifier", "gene.primaryIdentifier"
)

for row in query.results(row="rr"):
	o.write(row["gene.primaryIdentifier"]+'\t'+row["ortholog_gene.primaryIdentifier"]+'\n')
	



