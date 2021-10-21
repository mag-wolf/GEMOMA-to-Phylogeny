#!/usr/bin/python3

#this is extract_scos_Annotation.py of the Annotation-to-PHYLOGENY pipeline. A pipeline
#for comprehensice, protein based phylogenomic studies based on the GeMoMa pipeline (Keilwagen 2016/2018).
#This script will go through a provided directory full of orthogroups and find single copy ortholog sequences (SCOS).
#First it will go through the list "orthogroups" and remove entries with paralogs in the list "orthogroupsinglecopies"
#Then it will go through the list "orthogroupsinglecopies" and remove files with too many missing species in the list "orthogroupsresult"
#Eventually, it will copy the files of this final list "orthogroupsresult" into new files with the *_SCOS.fasta ending.
#(you can then remove every file with the *.fa ending) 

#Good luck!

#import modules
import os, sys, os.path

#get path to Orthogroup directory
path = sys.argv[1]

# create lists of orthogroups were we remove files that doesn't fulfill the criteria of a SCOS
orthogroups = os.listdir( path )
orthogroupsinglecopies = os.listdir( path )
orthogroupsresult = os.listdir( path )

#setup directories, lists and counter
keepO = {}
keepSC = {}
listO = []
countO = 0
countF = 0
counttest = 0

####get total number of species###
for itemO in orthogroups:
	pathOG = os.path.join(path, itemO)
	filehandleO = open(pathOG)
	for line in filehandleO:
		if line.startswith('>'):
			species,origin = line.split("-",1)
			if species not in keepO:
				keepO[species] = [0]
				listO.append(species)
	filehandleO.close()

for itemO in listO:
	countO = countO + 1

###get minimal allowed number of species###
exludecount = countO // 4
finalcount = countO - exludecount

###gather a list of single copies###
for itemS in orthogroups:
	pathSC = os.path.join(path, itemS)
	filehandleSC = open(pathSC)
	for line in filehandleSC:
		if line.startswith('>'):
			species,origin = line.split("-",1)
			if species not in keepSC:
				 keepSC[species] = [0]
			else:
				orthogroupsinglecopies.remove(itemS)
				orthogroupsresult.remove(itemS)
				break
	filehandleSC.close()
	keepSC = {}

###exclude SCOS with too many missing species (25% missing)###
for itemF in orthogroupsinglecopies:
	pathSCOS = os.path.join(path, itemF)
	filehandleSCOS = open(pathSCOS)
	for line in filehandleSCOS:
		if line.startswith('>'):
			countF = countF + 1
	if countF <= finalcount:
		orthogroupsresult.remove(itemF)
		pass
	countF = 0
	filehandleSCOS.close()

###copy SCOS to new files###
for itemR in orthogroupsresult:
	pathR = os.path.join(path, itemR)
	filehandleR = open(pathR)
	filenameN,rest = itemR.split(".",1)
	filenameN2 = filenameN + "_SCOS.fasta"
	pathN = os.path.join(path, filenameN2)
	filehandleN = open(pathN,"a")
	for line in filehandleR:
		filehandleN.write(line)
	filehandleR.close()
	filehandleN.close()

