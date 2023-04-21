#!/bin/bash
INPDIR=inputs

make

EXENAME=seqMDS

OUTDIR=out$INPDIR-$EXENAME-`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR
touch $OUTDIR/time.txt

#
# RUN SEQ MDS
#

#for f in X-n1001-k43 X-n979-k58 Golden_12 Golden_16 Brussels2 Flanders2 P-n101-k4 CMT5 # Flanders1 Flanders2  # X-n1001-k43 Golden_16 Flanders2 CMT5 # `ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
for file in `ls -Sr $INPDIR/*.vrp` 
do
  #file=$INPDIR/$f.vrp
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  
  #~ ROUNDING is enabled for all instances except Golden+CMT 
  isROUND=1  
  
  # ROUNDING is disabled if Golden or CMT
  if [[ $fileName = Golden* ]] || [[ $fileName = CMT* ]]; then
    isROUND=0
  fi
  
  ## EXECUTION if require run more times
  ./$EXENAME.out $file -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  #./$EXENAME.out $file -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  #./$EXENAME.out $file -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  #./$EXENAME.out $file -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  #./$EXENAME.out $file -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  
  echo $file - Done $isROUND
done
sort $OUTDIR/time.txt

#exit

# 
# RUN PAR MDS
# 

EXENAME=parMDS

OUTDIR=out$INPDIR-$EXENAME-`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

touch $OUTDIR/time.txt

# Num OMP Threads is default to #CPU CORES ! 
nTHREADS=`nproc --all`

# Else pick from args!
if [ $# -gt 0 ]; then
  nTHREADS=$1
fi

for file in `ls -Sr $INPDIR/*.vrp` 
do
  # GET the filename with extension removing foldername prefix.
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  
  #~ ROUNDING is enabled for all instances except Golden+CMT 
  isROUND=1  
  
  # ROUNDING is disabled if Golden or CMT
  if [[ $fileName = Golden* ]] || [[ $fileName = CMT* ]]; then
    isROUND=0
  fi
  
  # EXECUTION
  ./$EXENAME.out $file -nthreads $nTHREADS -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  
  echo $file - Done $isROUND
done

sort $OUTDIR/time.txt

# For my references
# DEL Later. Run smaller ones first. # C* X* P* *O*
#`ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`

