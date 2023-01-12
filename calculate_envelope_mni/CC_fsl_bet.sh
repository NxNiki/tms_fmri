#!/bin/bash

#inputdir=/media/sf_Linux_Shared/CausCon/CC_Porcupine_T1
inputdir=/media/sf_Linux_Shared/CausCon/T1/Selected
outputdir=$inputdir/FSL_BET

# example
#/usr/local/fsl/bin/bet $inputdir/CAUSCON_1013_RC1849_T1 $outputdir/CAUSCON_1013_RC1849_T1_brain -A -o

# loop through all available files
# for file in $inputdir/*_T1.nii
for file in $inputdir/*_t1.nii.gz
do

   filename="$(basename $file)"
   echo $filename

   #outputname=${filename/'.nii'/"_brain"}
   outputname=${filename/'.nii.gz'/"_brain"}
   #echo $outputname

   /usr/local/fsl/bin/bet $inputdir/$filename $outputdir/$outputname -A -o

done
