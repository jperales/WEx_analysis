#!/bin/bash

bed=$2;
bam=$1;

bedtools coverage -d -abam $bam -b $bed | \
	awk -F"\t" 'BEGIN{c=0;l=0;b=0;print "Avg_cov\tPerc_basesCov";}{c=c+$6;l=l+1;if($6!=0){b=b+1;}}END{print c/l,"\t",b/l*100;}';
