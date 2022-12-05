#!/bin/sh
## use bash to run the script.


in_dir="out01_roi_signals"
out_dir="out02_roi_signals_merge"

#out_file="NTHC_roi_signals_6mm.txt"
out_file="NTHC_roi_signals_10mm.txt"

mkdir -p $out_dir

echo "subj,group,roi,beta_value" > $out_dir/$out_file

#for file in $in_dir/tmp_*_6mm.txt; do
for file in $in_dir/tmp_*_10mm.txt; do
    echo $file
    file_name=$(echo $file | cut -d '/' -f 2)
    subj=$(echo $file_name | cut -d '_' -f 2)
    group=${subj:0:4}
    subjid=${subj: -4}
    roi=$(echo $file_name | cut -d '_' -f 3)_$(echo $file_name | cut -d '_' -f 4)
    echo $subjid,$group,$roi,$(cat $file)>>$out_dir/$out_file
done
