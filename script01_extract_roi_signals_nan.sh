#!/bin/sh
## use bash to run the script.

data_dir="Data_betavalues"
mask_dir="IndEnvMask"

out_dir="out01_roi_signals"

mkdir -p $out_dir

while read line; do
    echo "$line"
    subj=$(echo $line | cut -d '_' -f 2)
    roi=$(echo $line | cut -d '_' -f 4)_$(echo $line | cut -d '_' -f 5)
    echo $subj
    echo $roi

    img=$(ls $data_dir/*/*/$line)
    #echo $img

    subjid=${subj: -4}
    echo $subjid

    mask=$(ls $mask_dir/$roi/${subjid}_10mm.nii)

    fslmeants -i $img -o $out_dir/tmp_${subj}_${roi}_10mm.txt -m $mask
    
done < TEHC_filename.txt
#done < NTHC_filename.txt


