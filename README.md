

# Sequence Aligner 


### Description :
***
Bash script that uses STAR to align pair end reads to indexed reference genome and generates counts for the use of downstream analysis. Fastq files with respect to directions(R1 and R2), lane and flowcell are merged together and are convert to BAM files using STAR. Reads are compiled together by individual samples and are constructed into a matrix file that can be read by R. 

### References :
***
#### STAR
* [Installation](https://github.com/alexdobin/STAR)
* [Docuementation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

### Definitions: 
***
#### Inputs 
- @Param : $1 - `Meta Data` - Text file with List of barcode Ids/Samples 
- @Param : $2 - `Results Dir` - Output folder for results 
- @Param : $3 - `Data Dir` - Text file with paths to folder. New line for each flow cell. 
- @Param : $4 - `genome dir` - Path to indexed genome
- @Param : $5 - `Threads` - number of threads used in alignmnet 

#### Outputs 
- @return results_dir/results
- @return results_dir/counts 

### USAGE: 
***
First make sure you have an indexed genome. In order to do this go to [Ensemble](https://uswest.ensembl.org/index.html). Search for the species of interest and download the FASTA file and GTF file of interest under Gene annotation file. For the FASTA file, access the FTP site and navigate to the /dna directory. You wnat to to use the genome with all the chromosomes so use the primary assembly file.  For example if we are interested in the Xenopus genome be use [this](ftp://ftp.ensembl.org/pub/release-101/fasta/xenopus_tropicalis/dna/). 

You can downloaded it by copying the link of the file and using the `wget` command in terminal. 

```sh
$ wget ftp://ftp.ensembl.org/pub/release-101/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.d na.primary_assembly.1.fa.gz
```

Decompress the file 
```sh
$ gunzip  Xenopus_tropicalis.Xenopus_tropicalis_v9.1.d na.primary_assembly.1.fa.gz
```

Next download the GTF file in the Gene Annotation tabs in ensemble. There are several to choose from use the README to select which is best. Download the GTF using the FTP link 

```sh
$ wget ftp://ftp.ensembl.org/pub/release-101/gtf/xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.101.gtf.gz
```

Decompress the file 
```sh
$ gunzip  Xenopus_tropicalis.Xenopus_tropicalis_v9.1.101.gtf.gz
```

Next make sure STAR is properly installed see instructions [here](https://github.com/alexdobin/STAR) and you exported the command as to PATH

You can do this by 
```sh
$ export PATH=path/to/STAR-2.7.6a/source:$PATH
```

Create index genome by first making a folder to store our index 

```sh
mkdir index
```


Use STAR to create indexed genome. 
```sh
STAR --runThreadN 6 \
     --runMode genomeGenerate \
     --genomeDir index \
    --genomeFastaFiles path/to/fasta \
    --sjdbGTFfile path/to/gtf \
    --sjdbOverhang 99
```

You may increase the number of threads you use based of the machine. Check the number of threads you can use by;
```sh
$ ps -eo nlwp | wc -l
```

The number of threads are use will also require STAR to open more files and the `ulimit` will need to increase. To chekc the max limit use; 

```sh
$ ulimit -Hn
```

To increase the `ulimit` use;
```sh
ulimit -n Number you want to increase to
```

Before you run the aligner make sure Meta file is a list of the sample IDs that match the data and each new line is a new sample with no empty lines, like the following example;

```
sample_1_ACTATCTACT
sample_2_ACTATGTACT
sample_3_ACTATCTATT
```

Ensure the same with the data file is the same, that each new line is a path to the next flowcell; 
```
path/to/flowcell1
path/to/flowcell2
path/to/flowcell3
```

Create an output folder; 

```sh
mkdir example_results
```

Run the following 
```sh
bash  aligner_formal.sh meta.txt \
                        demo_res  \
                        data_file.txt \
                        path/to/index \
                        number_of_threads \ 
```

Once the script is running you will see prompts on the progress of the alignment. This may take sometime depending on size of the reads, the amount of samples and properties of your machine. 

```
Alignment Start  
------------------------------------

sample_1_ACTATCTACT
Organizing Data by Samples from each flow cell given
Unzipping and merging lanes
Using STAR to Align
Nov 13 16:48:08 ..... started STAR run
Nov 13 16:48:08 ..... loading genome
Nov 13 16:48:10 ..... started mapping
Nov 13 17:02:20 ..... finished mapping
Nov 13 17:02:21 ..... started sorting BAM
Nov 13 17:02:22 ..... finished successfully
```




