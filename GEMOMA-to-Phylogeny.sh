#!/bin/bash

#!/bin/bash
#
# A pipeline created to automatize all my steps for phylogenetic analyses using the GEMOMA annotation pipeline.
# In the end you'll have a concatenated supermatrix tree, a consensus supertrees and individual genetrees for
# all kinds of tree visualization tools!
#
#
# Does require specific conda environments!
#
# by Magnus Wolf 2021 (magnus.wolf@senckenberg.de)
#
# In the following part you need to provide paths and desired settings for one or all of the different scripts as well as some general
# tree constructing steps:
#
#################################general dependencies and parameters#####################################################################
source ~/anaconda3/etc/profile.d/conda.sh       	#make sure to be able to switch between conda envs
WORKINGDIR=/path/to/your/workingdir/			#path to the Working directory
ulimit -n 2048						#get more doable jobs-temp files, I suggest not to change that.
THREADS=20                                      	#number of threads
TASKS=20                                                #number of parallel task to run alignment or trimming processes
MAFFTenv=MAFFTenv                                       #name of the MAFFT conda env
CLIPKITenv=CLIPKITenv					#name of the CLIPKIT conda env
FASCONCAT=$WORKINGDIR"/FASconCAT-G_v1.04.pl"            #path to the FASconCAT-G perl script by Patrick Kueck from Koenig Museum Bonn
IQTREEenv=IQTREEenv                                     #name of the IQTREE conda env
ASTRAL=$WORKINGDIR"/Astral/astral.5.7.3.jar"            #path to the Astral java script for supermatrix trees
#
#
##################################GEMOMA-to-Phylogeny dependencies and parameters######################################################
GEMOMAenv=GEMOMAenv					#name of a conda env that includes the GEMOMA annotation pipeline
ORTHOFINDERenv=ORTHOFINDERenv				#name of a conda env that includes MAFFT, IQTREE, MCL, DIAMOND and scipy
FINDSHAREDOG=$WORLINGDIR"/../"
EXTRACTSCOS2=$WORKINGDIR"/extract_scos_Annotation.py"    #path to my python script extract_scos_Annotation.py, no change necessary when the script is in the working directory!
minPIperbrachen=1                                       #cutoff to kick out too conserved genes. The cutoff value should be a rough estimate of how many mutations per brach you want.
#                                                       #I used 1 but more is always better. The less species you use the higher you can go. I got ~500 genes with 30 species and a
#                                                       #PI-cutoff of 3. I would try a couple of different ones.
maxVariance=0.40                                        #cutoff to kick out too variavle genes. This should be a percentile cutoff and a rough estimation when a gene is to variable.
#                                                       #I chose 25% in variable sides based on my datasets. But again, I would try a couple of different cutoffs. The lower the better.
#                                                       #I would suggest to different combinations of both cutoffs used here and aim for ~500 resulting genes.
#
# Now that you provided the dependencies, you can choose which part of which script you wanne run, everything you'll call #TRUE# will run!
# You can find a detailed discription for each step in a makeshift square down below. Additionaly, this discription will pop up if you run the respective
# step.
#
##################################Options for GEMOMA-to-Phylogeny###################################################################
rungemoma=TRUE
orthocalling=TRUE
findsharedscos2=TRUE
makealignments2=TRUE
trimmgenealignments2=TRUE
filteralignments2=TRUE
concatgenealignments2=TRUE
trimmsupermatrix2=TRUE
constructgenetrees2=TRUE
constructsupermtree2=TRUE                                #this will take a loooooong time!
constructsuperttree2=TRUE
#
#########################################################################################################################################
# Now you finished. Go into the Working Directory and hit #GEMOMA-to-Phylogeny.sh#. Then everything you chosed will run one by one.
# I recommend using some sort of log pipeing since long pipelines like this one tend to produce a lot of hidden errors. Try e.g.
# #2>&1 | tee error.log# behind the bash command.
#
# After this line I wouldn't touch anything...

if [[ "$rungemoma" = TRUE ]]
        then
        echo "########################start rungemoma############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is rungemoma of GEMOMA-to-Phylogeny, a script for           #"
        echo "#comprehensive, protein based phylogenomic studies based on the   #"
        echo "#GeMoMa pipeline (Keilwagen 2016/2018). This part will run the    #"
	echo "#actual GeMoMa pipeline to get annotations for all of your        #"
	echo "#provided assemblies. Make sure that the genome assemblies of all #"
	echo "#species you wanne use are in the #ASSEMBLIES# directory. Now     #"
	echo "#you need to find as many gff or gff3 files for those species and #"
	echo "#place them in the #ANNOTATIONS# diretory in the working dir.     #"
	echo "#These gff files should contain the whole proteom and should      #"
	echo "#be named in the same fashion like the assemblies just with the   #"
	echo "# #gff# or #gff3# as a filename extension. Make sure to have      #"
	echo "#GeMoMa corretly installed in the GEMOMAenv! For now it will      #"
	echo "#always use 100Gb of RAM and the number of threads you provided   #"
	echo "#above. Will output one directory per assembly species in the     #"
	echo "# #GEMOMA# directory.                                             #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir GEMOMA
	cd GEMOMA
        conda activate $GEMOMAenv
	for file in ./../ASSEMBLIES/*.fasta;
	do
		ln -s $file
	done
        for file in ./../ANNOTATIONS/*.gf*;
        do
                ln -s $file
        done
        for file in *.gf*;
        do
		filename=$(echo $file | cut -d "." -f1)
                echo $filename".xyz "$file" "$filename".fasta" >> filelist.tmp.txt
        done
	for file in *.fasta;
	do
		filename=$(echo $file | cut -d "." -f1)
		seqcommandpart3=$(cat filelist.tmp.txt | grep -v "$filename" | tr "fasta\n" "fasta " | sed "s/xyz /xyz a=/g" | sed "s/gff /gff g=/g" | sed "s/gff3 /gff3 g=/g" | sed "s/fasta /fasta s=own i=/g" | sed "s/.xyz / /g" | sed "s/..$//g")
		GeMoMa GeMoMaPipeline -Xmx100g threads=$THREADS outdir=$WORKINGDIR"/"GEMOMA"/"$filename"_anno" GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$file i=$seqcommandpart3 2>&1 | tee "GeMoMa_"$filename"_anno.log"
	done
	rm filelist.tmp.txt
	conda deactivate
	cd ./../
        wait
        fi

if [[ "$orthocalling" = TRUE ]]
        then
        echo "########################start orthocalling#########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is orthocalling of GEMOMA-to-Phylogeny, a script for        #"
        echo "#comprehensive, protein based phylogenomic studies based on the   #"
        echo "#GeMoMa pipeline (Keilwagen 2016/2018). This part will run        #"
	echo "#OrthoFinder (Emms and Kelly 2019) to get orthologs. The          #"
	echo "#installation comes with this tool, make sure that the path is    #"
	echo "#correct. Apart of that you need a conda env with dependencies.   #"
	echo "#This conda env needs mafft, iqtree, mcl, diamond and scipy       #"
	echo "#correctly installed. Will output a txt file containing a list    #"
	echo "#of all found single copy orthologs.                              #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir ORTHOFINDER
        cd ORTHOFINDER
	mkdir PROTEOMS
	cd PROTEOMS
        for dir in ./../../GEMOMA/*_anno/;
        do
                dirname=$(echo $dir | cut -d "/" -f5)
                echo $dirname
                specname=$(echo $dirname | tr "/" "_" | cut -d "_" -f1,2)
                echo $specname
                cp $dir"predicted_proteins.fasta" ./
                mv predicted_proteins.fasta $dirname"_predicted_proteins.fasta"
                sed -i "s/>/>$specname-/g" $dirname"_predicted_proteins.fasta"
        done
	cd ./../
        conda activate $ORTHOFINDERenv
	python $WORKINGDIR"/OrthoFinder_source/orthofinder.py" -f PROTEOMS/ -t $THREADS -a $THREADS -M msa -A mafft -T iqtree -os 2>&1 | tee orthocalling.log
	conda deactivate
        cd ./../
        wait
        fi

if [[ "$findsharedscos2" = TRUE ]]
        then
        echo "########################start findsharedscos2######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is findsharedscos2 of GEMOMA-to-Phylogeny, a script for     #"
        echo "#comprehensive, protein based phylogenomic studies based on the   #"
        echo "#GeMoMa pipeline (Keilwagen 2016/2018). This part will search     #"
	echo "#single copy ortholog sequences (SCOS) in the Orthogroups found   #"
	echo "#by OrthoFinder. Make sure that the path to the                   #"
	echo "#python script #extract_scos_Annotation.py# is set correctly to   #"
	echo "#Working directory! Will allow 25% missing species representation #"
	echo "#per SCOS! If you want to change it, change the line              #"
	echo "#...exludecount = counter // 4... to whatever you desire. 4 means #"
	echo "#1/4 are allowed to be absent.                                    #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd ORTHOFINDER
	cp -r ./PROTEOMS/OrthoFinder/Results_*/Orthogroup_Sequences ./
	mv Orthogroup_Sequences SCOS
	rm -r ./PROTEOMS/OrthoFinder
	python $EXTRACTSCOS2 ./SCOS
	rm ./SCOS/*.fa
	cd ./../
        wait
        fi

if [[ "$makealignments2" = TRUE ]]
        then
        echo "########################start makealignments2######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is makealignments2 of GEMOMA-to-Phylogeny, a script for     #"
        echo "#comprehensive, protein based phylogenomic studies based on the   #"
        echo "#GeMoMa pipeline (Keilwagen 2016/2018). This part will run the    #"
	echo "#alignment building process using mafft and filters them a bit.   #"
        echo "#Requires a working conda env with mafft > v7 installed. Since    #"
        echo "#gene sequences are short, mafft is using the most accurate       #"
        echo "#option with 1000 iterative refinements incorporating local       #"
        echo "#pairwise alignment information. Afterwards it will also remove   #"
        echo "#every character of the IUPAC code that represents more then one  #"
        echo "#amino acid and will remove every MSA with unequal length.        #"
        echo "#FASTCONCAT can't handle those...This was the easiest method to   #"
        echo "#implement and there might be better ones!                        #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $MAFFTenv
        mkdir ALIGNMENTS2
        cd ALIGNMENTS2
        for file in ./../ORTHOFINDER/SCOS/*.fasta;
        do
                ln -s $file
        done
        task7(){
                 genename=$(echo $1 | cut -d "." -f1)
                 mafft --maxiterate 1000 --localpair $1 > $genename"_alig.fasta";
        }
        N=$TASKS
        (
        for file in *.fasta;
        do
                ((i=i%N)); ((i++==0)) && wait
                task7 "$file" &
        done
        )
        wait
        for file in *_alig.fasta;
        do
                filename=$(echo $file | cut -d "_" -f1,2)
                while IFS= read -r line;
                do
                        if [[ $line == ">"* ]];
                        then
                                echo $line | cut -d "-" -f1 >> $filename"_alig_noIUPAC.fasta"
                        else
                                echo $line | sed "s/J/-/g" | sed "s/Z/-/g" | sed "s/B/-/g" | sed "s/X/-/g" >> $filename"_alig_noIUPAC.fasta"
                        fi
                done < "$file"
        done
        for file in *_alig_noIUPAC.fasta;
        do
                filename=$(echo $file | cut -d "_" -f1,2)
                seqlength=0
                adder=0
                while IFS= read -r line;
                do
                        if [[ $line == ">"* ]];
                        then
                                echo $seqlength >> $filename"_noIUPAC_temp.txt"
                                seqlength=0
                        else
                                adder=$(echo $line | wc -c)
                                seqlength=$(($seqlength + $adder))
                        fi
                done < "$file"
                echo $seqlength >> $filename"_noIUPAC_temp.txt"
                uniqlengths=$(cat $filename"_noIUPAC_temp.txt" | uniq | wc -l)
                if [ $uniqlengths -gt 2 ];
                then
                        echo $file
                        echo $uniqlengths
                        rm $file
                fi
        done
        rm *_noIUPAC_temp.txt
        rm *_alig.fasta
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$trimmgenealignments2" = TRUE ]]
        then
        echo "########################start trimmgenealignments2#################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is trimmgenealignments2 of GEMOMA-to-Phylogeny, a           #"
	echo "#script for comprehensive, protein based phylogenomic studies     #"
	echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
	echo "#will  trimm all individual SCOS alignments in the ALIGNMENTS2    #"
        echo "#directory. Make sure to provide a working conda env with ClipKIT #"
        echo "#installed. The program will only keep parsimony-informative      #"
        echo "#sites and completely conserved sites. It also incorporates a     #"
        echo "#dynamic gap filtering. Will output timmed.fasta files.           #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd ALIGNMENTS2
        conda activate $CLIPKITenv
        task8(){
                 genename=$(echo $1 | cut -d "_" -f1,2)
                 clipkit $1 -m kpic-smart-gap -o $genename"_alig_trimmed.fasta"
                 clipkit $1 -m kpi-smart-gap -o $genename"_parinfosite.fasta"
        }
        N=$TASKS
        (
        for file in *alig_noIUPAC.fasta;
        do
                ((i=i%N)); ((i++==0)) && wait
                task8 "$file" &
        done
        )
        wait
        cd ./../
        wait
        fi

if [[ "$filteralignments2" = TRUE ]]
        then
        echo "########################start filteralignments2####################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is filteralignments2 of GEMOMA-to-Phylogeny, a              #"
        echo "#script for comprehensive, protein based phylogenomic studies     #"
        echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
        echo "#will infere how many individuals you used for your analysis      #"
        echo "#and then determine the minimal number of parsimony informative   #"
        echo "#sites you need to infere a meaninful tree. Afterwards it will    #"
        echo "#remove all aligments with a lower number of informative site.    #"
        echo "#This is especially important when using conserved sequences!     #"
        echo "#Eventually, it will remove temporary files called                #"
        echo "# #._parinfosite.fasta and test if you have enough genes to       #"
        echo "#proceed. Otherwise the pipeline will crash here!                 #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd ALIGNMENTS2
        counter=0
        for file in *_parinfosite.fasta;
        do
                genename=$(echo $file | cut -d "_" -f1,2)
                ntaxa=$(cat $file | grep ">" | wc -l)
                minsites1=$((2 * $ntaxa))
                minsites2=$(($minsites1 - 3))
                minsites3=$(($minsites2 * $minPIperbrachen))
                echo "the minimal number of parsimony informative sites neccessary for "$file" is "$minsites3"."
                onlyseq=$(cat $file | grep -v ">" | tr "\n" "_" | sed "s/_//g")
                nsites=$(echo $onlyseq | wc -c)
                nsites2=$(($nsites / $ntaxa))
                echo "we found "$nsites2" parsimony informative sites in "$file"."
                if [ "$nsites2" -lt "$minsites3" ];
                then
                        counter=$((counter+1))
                        echo "removing "$genename" from further analysis because it was too conserves."
                        rm $genename"_alig_trimmed.fasta" $genename"_parinfosite.fasta" $genename"_alig_noIUPAC.fasta"
                fi
                onlytotalsites=$(cat $genename"_alig_trimmed.fasta" | grep -v ">" | tr "\n" "_" | sed "s/_//g")
                ntotalsites=$(echo $onlytotalsites | wc -c)
                ntotalsites2=$(($ntotalsites / $ntaxa))
                maxsites=$(echo "scale=4; $ntotalsites2 * $maxVariance" | bc)
                maxsitesint=${maxsites%.*}
                echo "the maximal number of parsimony informative sites for "$file" is "$maxsitesint"."
                if [ "$nsites2" -gt "$maxsitesint" ] ;
                then
                        counter=$((counter+1))
                        echo "removing "$genename" from further analysis because it was too variable."
                        rm $genename"_alig_trimmed.fasta" $genename"_parinfosite.fasta" $genename"_alig_noIUPAC.fasta"
                fi
        done
        filecount=$(ls -lq *_parinfosite.fasta | wc -l)
        echo "filecount :"$filecount.
        echo "removed files :"$counter
        if [[ "$filecount" -lt 50 ]];
        then
                echo "you got too less sequences to proceed. This means that either"
                echo "your missing data filtering in the findsharedbusco step was to"
                echo "strict or that either the minPIperbrachen or the maxVariance"
                echo "cutoffs were to strict. Retry with less strict cutoffs or"
                echo "exclude some species."
                while true; do
                        read -p "Do you wish to proceed anyways?" yn
                        case $yn in
                                [Yy]* ) echo "nice try"; break;;
                                [Nn]* ) exit;;
                                * ) echo "Please answer yes or no.";;
                        esac
                done
        fi
        rm *_parinfosite.fasta
        cd ./../
        wait
        fi

if [[ "$concatgenealignments2" = TRUE ]]
        then
        echo "########################start concatgenealignments2################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is concatgenealignments2 of GEMOMA-to-Phylogeny, a          #"
        echo "#script for comprehensive, protein based phylogenomic studies     #"
        echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
        echo "#will concatenate all gene alignments in the ALIGNMENTS           #"
        echo "#directory using the FASconCAT_v1.11.pl perl script written by    #"
        echo "#Patrick Kueck from the Museum Koenig (https://github             #"
        echo "#.com/PatrickKueck/FASconCAT) with default options. Make sure to  #"
        echo "#provide the path to the script above. Will output a so called    #"
        echo "#FcC_supermatrix.fas in a dedicated SUPERMATRIX directory.        #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir SUPERMATRIX2
        cd SUPERMATRIX2
        for file in ./../ALIGNMENTS2/*alig_noIUPAC.fasta;
        do
                ln -s $file
        done
        perl $FASCONCAT -s
        rm *.fasta
        cd ./../
        wait
        fi

if [[ "$trimmsupermatrix2" = TRUE ]]
        then
        echo "########################start trimmsupermatrix2####################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is trimmsupermatrix2 of GEMOMA-to-Phylogeny, a              #"
        echo "#script for comprehensive, protein based phylogenomic studies     #"
        echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
        echo "#will trimm the supermatrix made in the step above using the      #"
        echo "#same clipkit command as used to trimm the genealignments. This   #"
        echo "#part is only necessary when no trimming was applied to the       #"
        echo "#genealignments. Otherwise it might be redundant. Use if desired. #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd SUPERMATRIX2
        conda activate $CLIPKITenv
        clipkit FcC_supermatrix.fas -m kpic-smart-gap -o FcC_supermatrix_trimmed.fas
        conda deactivate
        cd ./../
        wait
        fi

if [[ "$constructgenetrees2" = TRUE ]]
        then
        echo "########################start constructgenetrees2##################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is constructgenetrees2 of GEMOMA-to-Phylogeny, a            #"
        echo "#script for comprehensive, protein based phylogenomic studies     #"
        echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
        echo "#will construct individual gene trees for each of the             #"
        echo "#alignments in the ALIGNMENTS2 directory using IQtree. This       #"
        echo "#requires a working conda env with IQtree installed. Make sure    #"
        echo "#to provide the name of the env above. Will make 1000 bootstrap   #"
        echo "#replications and will use 5 threads since more will not increase #"
        echo "#the speed with so short sequences. You'll find constructed trees #"
        echo "#in the GENETREES directory.                                      #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $IQTREEenv
        mkdir GENETREES2
        cd GENETREES2
        for file in ./../ALIGNMENTS2/*trimmed.fasta;
        do
                ln -s $file
        done
        for file in *.fasta;
        do
                iqtree -s $file -bb 1000 -nt 5
        done
        rm *.fasta
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$constructsupermtree2" = TRUE ]]
        then
        echo "########################start constructsupermtree2#################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is constructsupermtree2 of GEMOMA-to-Phylogeny, a           #"
        echo "#script for comprehensive, protein based phylogenomic studies     #"
        echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
        echo "#will construct a supermatrix tree using IQtree and the           #"
        echo "#supermatrix alignment constructed above. Make sure that the      #"
        echo "#alignment is trimmed and that you provided the name of the       #"
        echo "#working conda env above. Will make 1000 bootstrap replications.  #"
        echo "#Constructed tree is lated placed in the SUPERMATRIX2 directory.  #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        conda activate $IQTREEenv
        cd SUPERMATRIX2
        iqtree -s FcC_supermatrix_trimmed.fas -bb 1000 -nt $THREADS
        cd ./../
        conda deactivate
        wait
        fi

if [[ "$constructsuperttree2" = TRUE ]]
        then
        echo "########################start constructsuperttree2#################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is constructsuperttree2 of GEMOMA-to-Phylogeny, a           #"
        echo "#script for comprehensive, protein based phylogenomic studies     #"
        echo "#based on the GeMoMa pipeline (Keilwagen 2016/2018). This part    #"
        echo "#will construct a consensus supertree tree from the gene          #"
        echo "#trees constructed in the steps above. This will be done by using #"
        echo "#Astral-III, a java script that takes a concatenated newick file  #"
        echo "#and can be downloaded at https://github.com/smirarab/ASTRAL.     #"
        echo "#Make sure to provide the path to the java script above. If your  #"
        echo "#java version isn't up to date create a conda env with a new      #"
        echo "#java installation and then run the overall script.               #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir SUPERTREE2
        cd SUPERTREE2
        for file in ./../GENETREES2/*.treefile;
        do
                ln -s $file
        done
        for file in *.treefile;
        do
                filename=$(echo $file | cut -d "_" -f1)
                mv $file $filename"_1000bt.nex"
        done
        cat *_1000bt.nex > GENETREES_concat.nex
        java -jar $ASTRAL -i GENETREES_concat.nex -o SUPERTREE.nex
        rm *1000bt.nex
        cd ./../
        wait
        fi
