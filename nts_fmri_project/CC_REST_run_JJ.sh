#!/bin/sh
##############SETUP ENVIRONMENT###################
export ETKINLAB_DIR=/share/PI/aetkin/
export ANALYSIS_PIPE_DIR=${ETKINLAB_DIR}/toolboxes/analysis_pipeline_cd/
export SPM8DIR=${ETKINLAB_DIR}/toolboxes/spm8_sge/
export PATH=${PATH}:${ETKINLAB_DIR}
export PATH=${PATH}:${ANALYSIS_PIPE_DIR}

######################
TARGETDIR="rest"
DATADIR=$SCRATCH/CausalConnectome/
TASKDIR=$DATADIR/rest/
STRUCTDIR=$DATADIR/structural/
MASKDIR=$HOME/scripts/CausalConnectome/ppi_masks_standard
masks="Hippo_AAL_atlas"



cd $TASKDIR


for i in `ls -r ${PWD}/CausCon_*.nii.gz`  ; do

    s=`basename $i | awk -F _ '{ print $2 }' | sed 's/CausCon//g'`
    echo ${s}

i=$(basename $i .nii.gz)
# if analyzed dir not exists, then do analysis below
if [ ! -d ${i}.${TARGETDIR} ] ; then

#find structural image and dir
T1=`ls ${STRUCTDIR}*${s}*_t1.nii.gz | sed -n '1p'`
T1=`remove_ext $T1`
if [ `imtest $T1` = 0  ]  || [ ! -d ${T1}.struct_only ]; then
echo "invalid image : $T1"
# exit 1
fi

if [ `imtest $T1` = 1 ] && [ -d ${T1}.struct_only ]; then
  for m in ${masks}; do
  sbatch -p aetkin,normal --time=3:00:00 --wrap=" analysis_pipeline_Jing.sh -resting -no_resting_gm_mask -func_data $i -t1 ${T1} -tr 2 \
  -deleteVolumes 6 -reg_info ${T1}.struct_only \
  -fc_rois_mni ${MASKDIR}/${m}.nii.gz -atlas_conn_opts --useAllLabels,--doCalcSingleLabel_FC \
  -motion -output_extension ${TARGETDIR}"
    echo RUNNING $TARGETDIR on $i $m
  done
fi

 fi #-rest_conn_only # -doGBC (global connectivity)
done
