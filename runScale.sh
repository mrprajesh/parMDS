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

#~ X-n1001-k43 X-n979-k58 Golden_12 Golden_16 Brussels2 Flanders2 P-n101-k4 CMT5
#~ X-n1001-k43 Golden_12 Flanders2 CMT5

for f in X-n1001-k43 Golden_12 Flanders2 CMT5 # Antwerp1  Antwerp2 Brussels1  Brussels2 Flanders1  Flanders2 Leuven1  Leuven2 Ghent1 Ghent2 # X-n1001-k43 X-n979-k58 Golden_12 Golden_16 Brussels2 Flanders2 P-n101-k4 CMT5
do
  file=$INPDIR/$f.vrp 
  # GET the filename with extension removing foldername prefix.
  fileName=$(echo $file | awk -F[./] '{print $(NF-1)}')
  
  #~ ROUNDING is enabled for all instances except Golden+CMT 
  isROUND=1  
  
  # ROUNDING is disabled if Golden or CMT
  if [[ $fileName = Golden* ]] || [[ $fileName = CMT* ]]; then
    isROUND=0
  fi
  
  for nTHREADS in 01 02 04 08 16 24 32 40
  do
    # EXECUTION
    ./$EXENAME.out $file -nthreads $nTHREADS -round $isROUND > $OUTDIR/$fileName-$nTHREADS.sol 2>> $OUTDIR/time.txt
    cat $OUTDIR/time.txt | tail -1 | cut -f1,5,9,11 -d' ' >> $OUTDIR/comm-time.txt
  done
  echo $file - Done $isROUND
done

cat $OUTDIR/comm-time.txt 
#sort $OUTDIR/time.txt #this sort will be required if inputs were not run in alphabetical order

# For my references
# DEL Later. Run smaller ones first. # C* X* P* *O*
#`ls -Sr $INPDIR/*.vrp`  # ls -Sr [CXP0] # `ls -Sr $INPDIR/*.vrp`
