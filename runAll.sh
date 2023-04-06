#!/bin/bash
INPDIR=inputs

make

EXENAME=seqMDS

OUTDIR=out$INPDIR$EXENAME`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

for file in `ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
do
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  ./$EXENAME.out $file > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  echo $file - Done
done

sort $OUTDIR/time.txt

EXENAME=parMDS

OUTDIR=out$INPDIR$EXENAME`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

for file in `ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
do
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  ./$EXENAME.out $file 20 > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
  echo $file - Done
done

sort $OUTDIR/time.txt
