#!/bin/bash

#Requires bowtie2 in the path
#Make a bowtie2 database, given a set of folders of genome sequences. Each folder is named for the organism. All of the FASTA files in that folder are chromosomes or plasmids for that organism
#Using parallel to speed things up
hash bowtie2 2>/dev/null || { echo >&2 "I require bowtie2 but it's not in the path.  Aborting."; exit 1;}
hash parallel 2>/dev/null || { echo >&2 "I require GNU parallel but it's not in the path.  Aborting."; exit 1;}
hash nucmer 2>/dev/null || { echo >&2 "I require nucmer but it's not in the path.  Aborting."; exit 1;}
hash show-coords 2>/dev/null || { echo >&2 "I require show-coords but it's not in the path.  Aborting."; exit 1;}


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #The scripts directory, which is in the same folder as the folder containing the raw genomes, database, etc.
#cd $DIR/.. 
nproc=24 #Number of concurrent processes to run

if (( ${#1} > 0 )); then #Optionally specify the database folder
	folder="$1"
else
	folder="./data/db"
fi
#folder="data/db"
if [ ! -d $folder ]; then 
echo "$folder not found"
echo "Please place the raw FASTA files in the desired output folder in a subfolder called 'raw_genomes'"; exit; fi 
raw_genomes=$folder/raw_genomes
if [ ! -d $raw_genomes ]; then 
echo "$raw_genomes not found"
echo "Please place the raw FASTA files in the desired output folder in a subfolder called 'raw_genomes'"; exit; fi 

all=$folder/all.fasta #Put all of the sequences in a single file (renamed for the organism)
names=$folder/all.txt #List of names
if [ ! -d $folder/temp ]; then mkdir $folder/temp; fi
#sim=$folder/sim
#if [ ! -d $sim ]; then mkdir $sim; fi

#Check to make sure that no names are duplicated
echo "Checking for duplicated names"
for f in $raw_genomes/*fasta; do f=${f#*raw_genomes/}; echo ${f%.fasta}; done > $names
dups=`sort $names | uniq -c  | awk '($1>1){print $2}'`
#for strain in `cat $names`; do
#	if (( `grep -c $strain $names` > 1 )); then
#		dups="$dups
#`grep $strain $names`"
#	fi
if (( ${#dups} > 0 )); then
	echo "Error in set of reference genome names: redundancy"
	echo $dups
	exit
fi

#Check to make sure that there aren't any empty sequences in the set of reference genomes
echo "Checking for empty records"
for f in $raw_genomes/*fasta; do 
	empties=`grep -B 1 "^$" $f|grep ">"`
	if (( ${#empties} > 0 )); then
		echo "Empty sequence found in ${f##*/}"
		echo "Specific record: $empties"
		exit
	fi
done

if [ -s $folder/all.null ]; then
	echo -n "Database already exists. Overwrite or exit (o/E)? "
	read bool
	if [ $bool == "O" ] || [ $bool == "o" ]; then echo "Overwriting"; rm -f $all $folder/db.1.bt2 $folder/all.null; else exit; fi 
fi

echo '#!/bin/bash
#Find the proportion of the two genomes that is unique
g1=$1
g2=$2
out=$3
nucmer -p $out $g1 $g2 2>/dev/null
show-coords -T -r -c -l -I 90 $out.delta > $out.coords 2>/dev/null #Minimum identity set to 90%

python << END
Rcoverage={} #Bases covered
Rlengths={} #Length of reference
Riden={} #Identity of region aligned
Qcoverage={}
Qlengths={} 
Qiden={}
f=open("$out.coords", "r")
def mean(x):
	if len(x) == 0:
		return(0)
	return(float(sum(x))/len(x))
for line in f.readlines():
	line=line.strip().split("\t")
	if len(line) == 13:
		Rchr=line[11]
		Qchr=line[12]
		if Rchr not in Rlengths.keys():
			Rlengths[Rchr]=int(line[7])
			Rcoverage[Rchr]=set([])
			Riden[Rchr]=[]
		if Qchr not in Qlengths.keys():
			Qlengths[Qchr]=int(line[8])
			Qcoverage[Qchr]=set([])
			Qiden[Qchr]=[]
		for q in range(int(line[0]), int(line[1])): 
			Rcoverage[Rchr].add(q)
			Riden[Rchr].extend([float(line[6])/100])
		for q in range(int(line[2]), int(line[3])): 
			Qcoverage[Qchr].add(q) 
			Qiden[Qchr].extend([float(line[6])/100])
if sum(Rlengths.values()) > 0:
	Rprop=sum([mean(Riden[x])*len(Rcoverage[x]) for x in Rcoverage.keys()])/float(sum(Rlengths.values()))
else:
	Rprop=0
if sum(Qlengths.values()) > 0:
	Qprop=sum([mean(Qiden[x])*len(Qcoverage[x]) for x in Qcoverage.keys()])/float(sum(Qlengths.values()))
else:
	Qprop=0
f.close()
f=open("$out", "w")
f.write(str(max(Rprop, Qprop)))
f.close()
END

rm $out.delta
'>unique_coverage.sh

#For each genome, compare it to the rest and calculate a single similarity value (proportion overlap)
align_all_genomes(){
	echo "Aligning all genomes"
	for d in $raw_genomes/*fasta; do 
		for r in $raw_genomes/*fasta; do
			if [[ $d < $r ]]; then 
				if [ ! -s $folder/temp/${d##*/}.${r##*/} ]; then 
					echo bash unique_coverage.sh $d $r $folder/temp/${d##*/}.${r##*/}
				fi
			fi
		done
	done | parallel -P $nproc #Process in parallel
}
align_all_genomes

#For each genome, check to see if it is unique enough to keep in the database
#For those that are not unique enough, merge together
echo "Finding redundant genomes"
declare -A links #Links of similar genomes
for d in $raw_genomes/*fasta; do 
	for r in $raw_genomes/*fasta; do
		if [[ $d < $r ]]; then 
			if [ -s $folder/temp/${d##*/}.${r##*/} ]; then 
				p=$(cat $folder/temp/${d##*/}.${r##*/})
			else
				p=0
			fi
			if [ ${p%e*} == $p ]; then
				if (( $(echo "$p > 0.99"|bc -l) )) ; then #Unique coverage below 1%
					echo "${r##*/} is too similar to ${d##*/}, merging"
					links["$r"]="${links["$r"]} $r $d"
					links["$d"]="${links["$d"]} $r $d"
				fi
			fi
		fi
	done
done

#Extend links to a second level of connection
for r in ${!links[@]}; do links["$r"]=$(for d in ${links["$r"]}; do echo ${links["$d"]}; done | tr ' ' '\n' | sort -u | tr '\n' ' ' ); done
#Extend links to a third level of connection
for r in ${!links[@]}; do links["$r"]=$(for d in ${links["$r"]}; do echo ${links["$d"]}; done | tr ' ' '\n' | sort -u | tr '\n' ' ' ); done
#Make names for each group
for r in ${!links[@]}; do
	if [ "$r" != "${links["$r"]}" ]; then
		fp=${r%/*}/
		genus=${r##*/}; genus=${genus%%_*}
		species=${r##*/}; species=${species#*_}; species=${species%%_*}
		new_name=$fp$genus\_$species\_$(for d in ${links["$r"]}; do echo ${d#$fp$genus\_$species\_}|sed 's/\.fasta//'; done | sed '/^$/d' | tr '\n' '-')
		new_name=${new_name%-}.fasta #Remove trailing hyphen
		cp $r $new_name
		mv $r $r\_OMIT
	fi
done

for r in ${!gps[@]}; do
	echo $r
	echo ${gps["$r"]}
done

for r in ${!gps[@]}; do
	if [ -s $r ]; then
		mv $r $folder/temp/${gps["$r"]} #Rename the file
	fi
done
align_all_genomes #Because the renamed genomes should be aligned as well

#Calculate the species and strain-unique proportion of each retained reference
echo "Calculating strain- and species-unique regions"
python << END
import os
strains=[f for f in os.listdir('$raw_genomes') if f.endswith('fasta')]
species=['_'.join(f.split('_')[0:2]) for f in strains]
def calc_unique(ref, queries): #Given one reference and a list of queries, find the amount of the reference that is unique
	coverage={} #Bases covered
	lengths={} #Length of reference
	iden={} #Percent identities of the bases aligned
	for q in queries:
		fp=''
		if os.path.exists('$folder/temp/'+ref+'.'+q+'.coords'): 
			fp='$folder/temp/'+ref+'.'+q+'.coords'
			pos=0 #The fields to grab from the coords file
			refname=11
			rlen=7
		if os.path.exists('$folder/temp/'+q+'.'+ref+'.coords'): 
			fp='$folder/temp/'+q+'.'+ref+'.coords'
			pos=2
			refname=12
			rlen=8
		if fp != '':
			f=open(fp, 'r')
			for line in f.readlines():
				line=line.strip().split('\t')
				if len(line) == 13:
					r=line[refname]
					if r not in lengths.keys():
						lengths[r]=int(line[rlen])
						coverage[r]=set([])
						iden[r]=1
					align_len=int(line[pos+1])-int(line[pos])
					align_iden=float(line[6])/100
					if align_len > 0 and align_iden > 0:
						iden[r]=iden[r]*len(coverage[r])+(align_iden*align_len)
						iden[r]=iden[r]/float(len(coverage[r])+align_len)
						coverage[r].update(set(range(int(line[pos]), int(line[pos+1]))))
	if len(lengths)==0:
		return(1)
	return(1-sum([iden[x]*len(coverage[x]) for x in coverage.keys()])/float(sum(lengths.values())))
#def mean(x):
#	if len(x) > 0:
#		if sum(x) == 0: return(0)
#		return(float(sum(x))/len(x))
#	return(0)
null=open('$folder/all.null', 'w')
for strain in sorted(strains):
	sp='_'.join(strain.split('_')[0:2]) 
	neighbor_strains=[ strains[q] for q in range(0, len(strains)) if species[q] == sp and strains[q] != strain ] #Same species
	neighbor_species=[ strains[q] for q in range(0, len(strains)) if species[q] != sp and strains[q] != strain ] #Different species
	if len(neighbor_strains) > 0:
		strain_p = calc_unique(strain, neighbor_strains)
	else:
		strain_p = float(1)
	if len(neighbor_species) > 0:
		species_p = calc_unique(strain, neighbor_species)
	else:
		species_p = float(1)
	strain=strain.replace('.fasta', '')
	sp=sp.replace('.fasta', '')
	null.write('\t'.join([strain, str(strain_p), sp, str(species_p)])+'\n')
END

echo '#!/bin/sh
d=$1
n=${d%.fasta} #Extension MUST be ".fasta"
n=${n#*raw_genomes/}
c=1 #Chromosome counter
cat $d|while read line; do 
	if (( ${#line} >0 )); then
		if [ ${line:0:1} == ">" ]; then
			echo ">$n.$c" #Rename the chromosomes for the organism they are from, plus a random counter
			c=`expr $c + 1`
		else
			echo $line
		fi
	fi
done
'>format_fasta.sh

if [ ! -s $all ]; then
	for d in $raw_genomes/*fasta; do echo $d #Each reference genome is in a single FASTA file (can have multiple sequences in the same file), named for the organism
		done | parallel -P $nproc -k bash format_fasta.sh > $all
fi

echo "Building bowtie index"
bowtie2-build $all $folder/db > /dev/null

if [ ! -s $folder/db.1.bt2 ]; then
	echo "Error building bowtie2 index"; exit; fi

#Clean up
#rm -r $folder/temp
rm -f $folder/all.fasta
rm -f $folder/*rpkm
rm -f $folder/*tab
rm -f format_fasta.sh
rm -f unique_coverage.sh
echo "Done"

exit
