# STRAIN IDENTIFICATION BY ALIGNMENT TO NEAR NEIGHBORS (SIANN)

SIANN was created to quickly identify a targeted set of organisms from a shotgun sequencing dataset. SIANN provides a fast answer and errs on the side of specificity rather than sensitivity. SIANN will only detect the organisms at >1% abundance that are also present within the reference database created by the user.

## Implementation
SIANN is written in Python but also calls out to an outside aligner, bowtie2. The reference database is formatted for use by bowtie2, and a crucial aspect of the reference genomes is that they must be named according to a strict convention of “Genus_species_strain”. That taxonomic information informs the algorithm when a read is aligned to a species- or strain-specific level, rather than spreading across multiple species or genera. 

## Usage
Running SIANN takes an input of a FASTQ file (or a pair of FASTQ files) and outputs a table with the top organisms that are present in that sample. 

Executing `siann.py` gives the following help message:

```
usage: siann.py [-h] [-d DB] [-t THREADS] [--paired PAIRED] [--report]
                [--reads_out] [--keep_sam] --reads READS --out OUT

required arguments:
--reads	Set of reads (FASTQ/FASTA) to be processed
--out	Prefix for output files

optional arguments:
  -h, --help				show this help message and exit
  -d DB, --db DB			database of reference genomes to use
  -t THREADS, --threads THREADS
                			number of threads to use for alignment (all by default)
  --paired PAIRED			second set of reads in pair (if any)
  --report              	turn off the generation of a report
  --reads_out           	turn on the output of species- and strain-specific reads
  --keep_sam            	retain the aligned reads in SAM format
```
 
Creating a reference database is a bit more involved. The following executables must be in the PATH:
  *	bowtie2
  *	parallel
  *	nucmer (part of mummer)
  *	show-coords (also part of mummer) 


## Installing

First, a set of reference genomes must be assembled, all in their own FASTA file named for the organism (as described above) in a folder named raw_genomes inside a folder named for the database. So, to create a database called ‘database,’ first create siann/database then create siann/database/raw_genomes and then add the FASTA files into that folder. Finally, execute `siann/scripts/make_database.sh database`. This will do many things, such as checking to make sure that there are no duplicate names, checking to make sure there are no empty sequences, whether such a database already exists, etc. Check the documentation inside SIANN/scripts/make_database.sh for more details. 
The process of making a database does two things. The trivial step is to create a bowtie database index. The harder step is to establish how much of each genome is species-specific and how much is strain-specific by aligning every genome against each other. The crucial pieces of the resulting database is a bowtie index and a file called all.null, which contains the proportion of each reference that is species- and strain-specific. Without those files, the database cannot be used. 

## Requirements
Can be installed with conda, using python2 (for now)

`conda create -n siann -c bioconda bowtie2 parallel mummer scipy python=2`

or:

`conda create -n siann -c bioconda siann`


## Database Construction
### Structure
Database directory can be specified at time of construction.  the database MUST contain a subdirectory named:
```
raw_genomes
```
Each genome to be used in the database must be named with 3 named parts (e.g. for bacteria, Genus, species, strain), separated by \_.  For example:

`Genus_species_strain.fasta`

Improperly named fastas will result in errors during analysis.

