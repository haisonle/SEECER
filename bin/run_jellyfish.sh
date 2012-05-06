#!/bin/bash
JF=$1
LCOUNT=$4

if [ "$#" -eq "4" ];
then
$JF count -m $3 -o jf_tmp -c 3 -s 10000000 -t 32 --both-strands $5
else
$JF count -m $3 -o jf_tmp -c 3 -s 10000000 -t 32 --both-strands $5 $6
fi;

# merge
N_TMP=`ls -1 jf_tmp_* | wc -l`
if [ $N_TMP -eq 1 ]
then
    mv jf_tmp_0 jf_merged_$3
else
    $JF merge jf_tmp_* -o jf_merged_$3
    rm jf_tmp_*
fi

$JF dump --lower-count=$LCOUNT -o $2 -c jf_merged_$3
rm jf_merged_$3
