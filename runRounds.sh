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

for file in `ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
do
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  
  #~ ROUNDING is enabled for all instances except Golden+CMT 
  isROUND=1  
  
  # ROUNDING is disabled if Golden or CMT
  if [[ $fileName = Golden* ]] || [[ $fileName = CMT* ]]; then
    isROUND=0
  fi
  
  ## EXECUTION
  ./$EXENAME.out $file -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  
  echo $file - Done
done
sort $OUTDIR/time.txt

EXENAME=parMDS

# 
# RUN PAR MDS
# 

OUTDIR=out$INPDIR-$EXENAME-`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

touch $OUTDIR/time.txt

# Num OMP Threads is default to 20! 
nTHREADS=20

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
