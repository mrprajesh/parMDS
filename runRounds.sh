#!/bin/bash
INPDIR=inputs

make

EXENAME=seqMDS

OUTDIR=out$INPDIR$EXENAME`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

#
# RUN SEQ MDS
#
#~ for file in `ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
#~ do
  #~ fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  #~ ./$EXENAME.out $file > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  #~ echo $file - Done
#~ done

#~ sort $OUTDIR/time.txt

EXENAME=parMDS

# 
# RUN PAR MDS
# 

OUTDIR=out$INPDIR-$EXENAME-`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

#`ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
for file in $INPDIR/*.vrp
do
  # GET the filename with extension removing foldername prefix.
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  
  # ROUNDING is enable for all instances except Golden+CMT 
  if [[ $fileName = Golden* ]] || [[ $fileName = CMT* ]]; then
    echo $fileName round 0
    ./$EXENAME.out $file -nthreads 20 -round 0 > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  else
    ./$EXENAME.out $file -nthreads 20 -round 1 > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  fi
  echo $file - Done
done


sort $OUTDIR/time.txt
