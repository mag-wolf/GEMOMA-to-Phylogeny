# GEMOMA-to-Phylogeny
-------------------------------

A pipeline created to automatize all my steps for phylogenetic analyses using the GEMOMA (Keilwagen et al. 
2019) annotation pipeline.

Does require perl5, python3 and anaconda environments!

by Magnus Wolf 2021 (magnus.wolf@senckenberg.de)
-------------------------------
M.Sc. Magnus Wolf

PhD-Student

Senckenberg Biodiversity and Climate Research Center (SBiK-F)

Group: Evolutionary Vertebrate Genomics

Georg-Voigt-Straße 14-16, floor 4, room 4.07

60325 Frankfurt am Main, Germany

Installation:
-------------------------------
1.) Download the tarball

2.) copy the tarball and rename it to what you desire

3.) extract the renamed tarball

    tar -xzvf renamed.tar.gz renamed
    
4.) Go into the extracted directory and extract all sub-directories:

    tar -xzvf ANNOTATIONS.tar.gz ANNOTATIONS
    tar -xzvf ASSEMBLIES.tar.gz ASSEMBLIES
    tar -xzvf Astral.tar.gz Astral
    tar -xzvf OrthoFinder_source.tar.gz OrthoFinder_source

4.) install dependencies via conda:

    conda create --name MAFFTenv
    conda install -n MAFFTenv -c bioconda mafft 
    conda create --name CLIPKITenv
    conda install -n CLIPKITenv -c jlsteenwyk clipkit
    conda create --name IQTREEenv
    conda install -n IQTREEenv -c bioconda iqtree
    conda create --name GEMOMAenv
    conda install -n GEMOMAenv -c bioconda gemoma
    conda create --name ORTHOFINDERenv                         #This environment doesn't actually contain the orthofinder binary but only sets up the environment to run it!
    conda install -n ORTHOFINDERenv -c bioconda mcl
    conda install -n ORTHOFINDERenv -c bioconda mafft
    conda install -n ORTHOFINDERenv -c bioconda iqtree
    conda install -n ORTHOFINDERenv -c bioconda diamond
    conda install -n ORTHOFINDERenv -c anaconda scipy

Now you are ready to start.

Usage:
-------------------------------
1.) Gather whole genome assemblies in "fasta" format for all species you want to have in your tree. Make sure to include an
outgroup! Copy them in the dedicated and empty "ASSEMBLIES" folder that you extracted previously. For convenience, I would 
recommend naming the fasta files after the scientific species name. E.g.: Balaenoptera_physalus.fasta, Balaenoptera_musculus.fasta etc.... 
Then try to gather as many whole proteome gff-files for the chosen species and put them into the 
dedicated and emtpy "ANNOTATIONS" folder that you extracted previously. Make sure that the files are named EXACTLY like
the assemblies just with a "gff" extension instead of a "fasta" extension. 

2.) Open the script GEMOMA-to-Phylogeny.sh with a text editor of your choice (e.g. nano).

3.) Edit general dependencies. Especially the working directory and the number of threads of parallel processes you wanne use.

4.) Edit GENOMA-to-Phylogeny dependencies and parameters. Especially the filter settings for too variable too conserved genes.

5.) Edit options for GEMOMA-to-PHYLOGENY. Call everything "TRUE" that you want to use. The pipeline 
contains 11 subparts that can be run independently if all other subparts are called "FALSE". By
leaving it as it is, everything will run one by one. 

6.) Now simply run:
    
    bash GEMOMA-to-Phylogeny.sh

I suggest piping screen outputs to an error log by adding: 

    2>&1 | tee error.log
 

Here a list of all subparts:
###

    rungemoma                   #run the genoma annotation pipeline for every assembly using all homology-information provided in the ANNOTATIONS directory
    orthocalling                #find orthologs between the annotations made by gemoma
    findsharedscos2             #find shared single copy orthologs (SSCO) with a 25% missing species cutoff
    makealignments2             #make alignments of all SSCO sequences using mafft
    trimmgenealignments2        #trimm the alignments using clipkit
    constructgenetrees2         #constructing gene trees of every single SSCO alignment using iqtree
    filtergenetrees2            #filter genetrees and alignments based on the maximum likelihood distance (throw out too conserved and too variable genes)
    concatgenealignments2       #concatenate gene SSCO into one big matrix using FASconCAT
    constructsupermtree2        #constructing a tree from the concatenated alignment using iqtree
    constructsuperttree2        #constructing a consensus tree based on all constructed SSCO trees using Astral 

###


Good luck!
-------------------------------
