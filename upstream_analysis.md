
#pipeline for upstream QC analysis 

## Index the reference genome 

* download the updated assmbly (Mmul_10) AND annotation in GTF format here: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/`
* Run STAR as follows:

### Start interactive session with enough memory

```
srun -p highmem --pty -c 8 --mem 300G -t 0-05:00 /bin/bash 
```

### activate your conda environment and print working directory

```
conda activate [environment-name]
/n/data1/hms/dbmi/farhat/ba157/reference/Mmul_10
```

### unzip the reference assembly and annotation

```
gunzip *gz
```

### run STAR to generate genome index

``` 
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles GCF_003339765.1_Mmul_10_genomic.fna --sjdbGTFfile GCF_003339765.1_Mmul_10_genomic.gtf
```

## Align reads to reference

### Make a folder for the alignment results

```
mkdir /n/data1/hms/dbmi/farhat/ba157/Hansen-2018-bulkRNAseq/STAR-results
```


### Run STAR with wrap on all fastq files

A trial run took 1:22 hours on 1 core with 28G memory for oen sample.


```
for R1 in `ls /n/data1/hms/dbmi/farhat/ba157/Hansen-2018-bulkRNAseq/fastqc_trimmed/*_1_*`; do FILE=${R1//_1_trimmed.fastq.gz/}; PREFIX=${FILE//\/n\/data1\/hms\/dbmi\/farhat\/ba157\/Hansen-2018-bulkRNAseq\/fastqc_trimmed\/SRR/}; sbatch --mem=50000 -t 0-03:00:00 -c 1 -p short --wrap="
STAR --runThreadN 1 --genomeDir /n/data1/hms/dbmi/farhat/ba157/reference/Mmul_10 --readFilesIn ${FILE}_1_trimmed.fastq.gz ${FILE}_2_trimmed.fastq.gz --readFilesCommand gunzip -c --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix SRR${PREFIX}"; done
```


### Concatenate gene counts to one matrix 

* we chose the unstranded raw gene counts (i.e. second column of the ReadsPerGene output)
* we use awk to add the isolate name in the first row

```
for file in *ReadsPerGene.out.tab; do awk 'NR==1 { print "column", FILENAME }; 1' "$file" | cut -f1,2 > temp && mv temp "$file"; done
```


* we use sed to remove the file extension

```
sed -e 's/_ReadsPerGene.out.tab//g;s/column[[:space:]]//g' -i *ReadsPerGene.out.tab
```

* make new files with only the counts to combined them in the next step

```
for i in `ls *ReadsPerGene.out.tab`; do cut -f2 $i > $i.counts; echo $i; done
```

* combine all files into one large matrix 

```
paste *counts > Hansen-gene-counts.tsv
paste Hansen-gene-counts_column_names.txt Hansen-gene-counts.tsv > Hansen-gene-counts.final.tsv
```

* final gene count table is here

```
/n/data1/hms/dbmi/farhat/ba157/Hansen-2018-bulkRNAseq/gene_counts
```
