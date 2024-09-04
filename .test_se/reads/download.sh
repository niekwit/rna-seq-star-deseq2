#!/usr/bin/env bash
#https://www.science.org/doi/full/10.1126/science.1227919
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/SRR592966/SRR592966.fastq.gz -o 5hr_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/SRR592965/SRR592965.fastq.gz -o 72hr_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/SRR592963/SRR592963.fastq.gz -o 5hr_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/SRR592968/SRR592968.fastq.gz -o 72hr_2.fastq.gz

#Only first 1E6 lines are kept for all samples (head -1000000)
