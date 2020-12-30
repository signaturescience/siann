# Example workflow to use SIANN
## Background
SIANN is intended to be used to determine lineage of a member of a metagenomic community at the species or strain level.  It is expected that alternative methods (e.g. METSCALE) have been used to classify the community at the genus or species level.  Once the species has been determined, if further detail is desired, SIANN can be used in a 2step process:
## Database construction

1) Place all downloaded genomes (named Genus_Species_Strain.fasta) into a subfolder named data\db\raw_genomes
  ```
  mkdir -p data\db\raw_genomes 
  cd data\db\raw_genomes
  ```
2) Download all finished isolates of bacterial species to be classified. 
  * Downloaded references from NCBI will have to be renamed. 
  * Species with < 1000 sequenced isolates can be downloaded directly with ftp, e.g.
    *
    ```
    lftp -c 'mget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Mycobacterium_aquaticum/*'
    ```
  * Species with > 1000 sequenced isolates require parsing assembly_summary.txt file.  It is also advisable to reduce total number of genomes, to reduce database generation time. 
     * For example, to get the first 50 entries from assembly_summary.txt (in this case Staphylococcus aureues): 
     ```
     lftp -c "mget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/assembly_summary.txt" 
     head -n 52 assembly_summary.txt | tail -50 | awk -F\\t '{ gsub("[: =]","",$9); gsub("strain","",$9); printf("%s\t%s\t%s\t%s\n",$1,$9,$16,$20)}' > selectedStrains.txt 
     while IFS=$'\t' read -r accession strain asm location;
     do 
     lftp -c "mget ${location}/${accession}_${asm}_genomic.fna.gz"; 
     strain=${strain//[^[:alnum:]]/}
     zcat ${accession}_${asm}*.gz > Staphylococcus_aureus_${strain}.fasta
     done < selectedStrains.txt

     ```
3) Format DB:
```
cd ../../
./make_database.sh
```
   
# Test with data
1) Test with HMP_mock_even data:
```
mkdir data/reads
cd data/reads
wget https://osf.io/7cbvh/download 
cd ../../
./siann.py -d data/db --reads data/reads/HMP_even.trimQ2.fastq.gz --out siann_test_1 -t 20
```
2) Data output will be in directory.

  
