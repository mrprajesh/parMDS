#!/bin/bash
INPDIR=inputs

make

EXENAME=parMDS

# 
# RUN PAR MDS
# 

OUTDIR=out$INPDIR-$EXENAME-`date +%d-%b-%Y-%H%M%S`
mkdir -p $OUTDIR

touch $OUTDIR/time.txt

#~ # Num OMP Threads is default to 20! 
#~ nTHREADS=20

#~ # Else pick from args!
#~ if [ $# -gt 0 ]; then
  #~ nTHREADS=$1
#~ fi

#~ X-n1001-k43 X-n979-k58 Golden12 Golden16 Brussels2 Flanders2 P-n101-k4 CMT5
#~ X-n1001-k43 Golden16 Flanders2 CMT5

for f in X-n1001-k43 Golden16 Flanders2 CMT5 #`ls -Sr $INPDIR/*.vrp` 
do
  file =$INPDIR/$f.vrp 
  # GET the filename with extension removing foldername prefix.
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  
  #~ ROUNDING is enabled for all instances except Golden+CMT 
  isROUND=1  
  
  # ROUNDING is disabled if Golden or CMT
  if [[ $fileName = Golden* ]] || [[ $fileName = CMT* ]]; then
    isROUND=0
  fi
  
  for nTHREADS in 1 2 4 8 16 24 32 40
  do
    # EXECUTION
    ./$EXENAME.out $file -nthreads $nTHREADS -round $isROUND > $OUTDIR/$fileName.sol 2>> $OUTDIR/time.txt
    cat $OUTDIR/time.txt | tail -1 | cut -f1,5,9 -d' ' >> $OUTDIR/comm-time.txt
  done
  echo $file - Done $isROUND
done

sort $OUTDIR/time.txt

# For my references
# DEL Later. Run smaller ones first. # C* X* P* *O*
#`ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
