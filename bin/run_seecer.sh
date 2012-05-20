#/usr/bin/bash
#
# Copyright (C) 2012  Hai-Son Le (haisonle@gmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# This script runs the SEECER pipeline of 4 steps:
#
# 1. Replace Ns and strip off read IDs (to save memory).
# 2. Run JELLYFISH to count kmers.
# 3. Correct errors with SEECER.
# 4. Clean up and put back original read IDs.
#


BINDIR='./bin/'
JF="../jellyfish-1.1.4/bin/jellyfish"
K=17
SEECER_PARAMS=""
SeecerStep=1
LCOUNT=3
TMPDIR=''

usage=$(cat << EOF
   # This script runs the SEECER pipeline of 4 steps:
   #
   # 1. Replace Ns and strip off read IDs (to save memory).
   # 2. Run JELLYFISH to count kmers.
   # 3. Correct errors with SEECER.
   # 4. Clean up and put back original read IDs.
   
   run_seecer.sh [options] read1 read2

   read1, read2: are Fasta/Fastaq files.
	  If only read1 is provided, the reads are considered singles.
          Otherwise, read1 and read2 are paired-end reads.

   Options:
      -t <v> : *required* specify a temporary working directory.
      -k <v> : specify a different K value (default = 17).
      -j <v> : specify the location of JELLYFISH binary (default = $JF).
      -p <v> : specify extra SEECER parameters (default = '').
      -s <v> : specify the starting step ( default = 1). Values = 1,2,3,4.
      -h : help message
EOF
);

while getopts ":j:p:k:s:t:h" opt; do
  case $opt in
    t)
      TMPDIR=$OPTARG
      ;;  
    k)
      K=$OPTARG
      ;;
    j)
      JF=$OPTARG
      ;;
    p)
      SEECER_PARAMS=$OPTARG
      ;;
    s)
      SeecerStep=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "$usage"
      exit 1;
      ;;
    h)
       echo "$usage"
       exit 1;
       ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1;
      ;;
  esac
done

shift $(($OPTIND - 1))

if [ -z "$TMPDIR" ];
then
    echo "Missing -t parameter: please specify a temporaty working directory to store JELLYFISH output.";
    exit 1;
fi


Read1=$1
Read2=$2

if [ ! -z "$1" ];
then
    Read1_N="${Read1}_N"
    RS_ARGS="$Read1,${Read1}_N"
    Reads="$Read1"
    Reads_N="${Read1}_N"
    Reads_O="${Read1}_corrected.fa"
else
    echo "$usage"
    exit 1;
fi;

if [ ! -z "$Read2" ];
then
    Read2_N="${Read2}_N"
    RS_ARGS="$RS_ARGS $Read2,${Read2}_N"
    Reads="$Reads,$Read2"
    Reads_N="${Reads_N},${Read2}_N"
    Reads_O="${Reads_O},${Read2}_corrected.fa"
fi;

echo "K is set to $K" >&2
echo "JELLYFISH bin is set to $JF" >&2
echo "SEECER parameters: $SEECER_PARAMS" >&2
echo "Starting SEECER from step: $SeecerStep" >&2

echo "------------------------" >&2
echo "Read 1 = $1" >&2
echo "Read 2 = $2" >&2
echo "------------------------" >&2



##
# Seecer runs in several steps, output of each step could be saved
# so that Seecer can starts at the next step
#
##

# 1. Remove read IDs to save memory
# Output: for each read files a temp file
if [ $SeecerStep -le 1 ];
then
    echo "++ Step 1: Replacing Ns ... and stripping off read IDs"
    echo
    ${BINDIR}/random_sub_N $RS_ARGS
fi;

# 2. Running JELLYFISH to count kmers
# Output: counts_$SeecerK file

if [ $SeecerStep -le 2 ];
then
    echo "++ Step 2: Running JELLYFISH to count kmers ..."
    echo
    bash ${BINDIR}/run_jellyfish.sh $JF $TMPDIR/counts_${K}_${LCOUNT} $K $LCOUNT $TMPDIR $Read1_N $Read2_N
fi;

# 3. Seecer main part: correcting errrors
# Output: corrected.fa
if [ $SeecerStep -le 3 ];
then
    echo "++ Step 3: Correcting errors with SEECER ... Your reads are in good hands!"
    echo "-----------------------------------------------------------------------"
    echo " *** Start time: " `date`;

    ${BINDIR}/seecer $Read1_N $Read2_N $SEECER_PARAMS --kmer $K -k $TMPDIR/counts_${K}_${LCOUNT} -o $TMPDIR/corrected.fasta
    echo " *** End time: " `date`;
    echo "-----------------------------------------------------------------------"
    echo
fi;

# 4. Put back the original read IDs
# Output:
if [ $SeecerStep -le 4 ];
then
    echo "++ Step 4: Cleaning and putting back original read IDs ... We finish soon!"
    ${BINDIR}/replace_ids $TMPDIR/corrected.fasta $Reads $Reads_N $Reads_O
    rm $TMPDIR/corrected.fasta
fi;

