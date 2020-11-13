#!/bin/bash

# Align Broad Smart Seq2 Data 

echo "Alignment Start"
echo "------------------------------------"


#$1 - Meta data, list of barcodes/samples 

#list files given in meta

Meta_data=$1
results_dir=$2
data_dir=$3
genome_dir=$4
threads=$5

echo " "
mkdir "${results_dir}/results"
cp $1 "${results_dir}/results"
cd "${results_dir}/results"

#Read Meta file  and  go through each sample
while read p; do
  echo "$p"
  #Make folder for sample name
  mkdir "${results_dir}/results/${p}"
  #copy files over to folder
  regrex="${p}.unmapped.[0-9]"

  echo  "Organizing Data by Samples from each flow cell given"

  while read d; do 
  #Flowcell #1 
  	cd $d 
# cd /home/vasanthchandrasekhar_g_harvard_e/scData/data/H2KH7BGXG/get.broadinstitute.org/pkgs/SN0203980
  	ls | grep "$regrex" | xargs -I '{}' cp '{}' "${results_dir}/results/${p}"
  done < $data_dir



 
#Flowcell #2
# cd /home/vasanthchandrasekhar_g_harvard_e/scData/data/H2JKVBGXG/get.broadinstitute.org/pkgs/SN0203975
# ls | grep "$regrex" | xargs -I '{}' cp '{}' "${results_dir}/results/${p}"

  

  #Merge files from different lanes
  cd "${results_dir}/results/${p}"


  echo "Unzipping and merging lanes"
  gunzip *
  merge_R1="merge_${p}.1.unmapped.fastq"
  merge_R2="merge_${p}.2.unmapped.fastq"
  cat *unmapped.1.fastq > "$merge_R1"
  cat *unmapped.2.fastq > "$merge_R2"

 #Align sequences using STAR

  echo "Using STAR to Align"
  #mkdir "${results_dir}/results/${p}_star"
  STAR --runThreadN $threads \
  --genomeDir $genome_dir \
  --readFilesIn "$merge_R1"  "$merge_R2" \
  --outFileNamePrefix "${results_dir}/results/${p}/${p}_star/" \
  --quantMode GeneCounts \
  --outSAMtype BAM SortedByCoordinate 

  mv "${results_dir}/results/${p}/${p}_star/ReadsPerGene.out.tab" ${p}_ReadsPerGene.out.tab  
  cd "${results_dir}/results" 
done <$Meta_data

echo "out of loop"

#Creating counts
#cd "${results_dir}/results" 
#mkdir "${results_dir}/results/counts" 
#for x in */*ReadsPerGene.out.tab; do  s=`basename $x | cut -f1 -d.`;  echo $s;  cat $x | tail -n +5 | cut$


