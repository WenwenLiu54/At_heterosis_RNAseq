#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# RNAseq_part2_salmon.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-11-3
# Usage: nohup sh RNAseq_part2_salmon.sh /data2/usr/LiuWW/project_1/RNAseq_part1_quality_control_result > RNAseq_part2_salmon.log 2>&1 &

#------------------------------------------------------------------------------------

## parameters required to change (based to different machines' limitation)

THREADS="20"

#------------------------------------------------------------------------------------

## Locations require to change (based to different computers) && Assign the I/O variables

RAW_LOCATION=${1}   # /data2/usr/LiuWW/project_1/RNAseq_part1_quality_control_result

INDEX="/data2/usr/LiuWW/project_1/try/AtRTD2_index_2"

TOTAL_OUTPUT="/data2/usr/LiuWW/project_1/try/RNAseq_part2_salmon_output_2"

#------------------------------------------------------------------------------------

## mkdir new directories
if [ ! -d ${INDEX} ]
    then mkdir -p ${INDEX}
fi

if [ ! -d ${TOTAL_OUTPUT} ]
    then mkdir -p ${TOTAL_OUTPUT}
fi

#------------------------------------------------------------------------------------

## Build index (--type The type of index to build: the only option is "puff" in this version 1.0.0 of salmon)

salmon index \
	-t /data2/usr/LiuWW/project_1/try/AtRTD2/AtRTDv2_QUASI_19April2016.fa \
    -i ${INDEX}/AtRTDv2_QUASI_19April2016_index \
	--type puff \
	-k 31 \
	> ${INDEX}/salmon_index.log 2>&1

wait && echo "Build index done!" 

## Quantification

for i in `cat /data2/usr/LiuWW/project_1/sample_name_list.txt`
do

SAMPLE_NAME=${i}

FASTQ1="${RAW_LOCATION}/${i}/clean.${i}_R1.fq.gz"
FASTQ2="${RAW_LOCATION}/${i}/clean.${i}_R2.fq.gz"

#------------------------------------------------------------------------------------

## Locations require to change (based to different computers) && Assign the I/O variables

OUTPUT="${TOTAL_OUTPUT}/${SAMPLE_NAME}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

## mkdir new directories

if [ ! -d ${OUTPUT} ]
    then mkdir -p ${OUTPUT}
fi
    
#------------------------------------------------------------------------------------

salmon quant \
	-i ${INDEX}/AtRTDv2_QUASI_19April2016_index \
	-l IU \
	-1 ${FASTQ1} \
	-2 ${FASTQ2} \
	-p ${THREADS} \
	--numBootstraps 30 \
    --seqBias \
    --validateMappings \
    --writeUnmappedNames \
    -o ${OUTPUT}/${SAMPLE_NAME}_transcripts_quant \
	> ${OUTPUT}/${SAMPLE_NAME}_salmon_quantification.log 2>&1
done	

wait && echo " Quantification is done!"	
	
