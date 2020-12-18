#!/bin/bash

### Andre Hechler, 29.05.2019
### registration of structural images to DWI

Usage() {
    cat <<EOF

Usage: pipe_6_tractography <session_list.csv>

  <input> is the session as csv file
  script needs exactly one input file

EOF
    exit 1
}

# set error message if code fails within "set -e"
abort() {
    echo "***Tractography was stopped due to an error. Please revise code.***"
    exit 1
}

trap 'abort' 0 

[ "$#" -ne 1 ] && Usage
#Checks if number of inputs is unequal to 1 (one csv file). If this evaluates to #true, the usage message above is returned

#####################################################
#### Inputs #########################################
#####################################################

session_list=$1

while IFS=',' read SUBJECT SESSION; do

PROJECT_DIRECTORY=/home/andre/Desktop/Thesis/01_data

SUBJECT_PATH=${PROJECT_DIRECTORY}/derivatives/sub-${SUBJECT}
SESSION_PATH=${PROJECT_DIRECTORY}/derivatives/sub-${SUBJECT}/ses-${SESSION}

# anat folder and anat files
ANAT_PATH=${SESSION_PATH}/anat
DWI_PATH=${SESSION_PATH}/dwi

DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif")
FOD_FILE=$(find ${DWI_PATH} -maxdepth 1 -name "*dwi_fod.mif")
PARC_FILE=$(find ${ANAT_PATH} -name "aparc+aseg_regDWI_labels.mif")
SEG_FILE=$(find ${ANAT_PATH} -name "*T1w_seg_DWI.mif")
GMWMI_FILE=$(find ${ANAT_PATH} -name "*gmwmi_regDWI.mif")

TCK_OUT="${DWI_FILE%dwi.mif}tracts.tck"
TCKSFT_OUT="${DWI_FILE%dwi.mif}tracts-sift.tck"

CONN_OUT="${DWI_FILE%dwi.mif}connectome_zeros.csv"

#####################################################

# anatomically constrained tractography, seeding at GM-WM interface with 20m streamlines
tckgen ${FOD_FILE} ${TCK_OUT} -act ${SEG_FILE} -seed_gmwmi ${GMWMI_FILE} -select 20M

echo "#########################################################"
echo "tckgen for ${SUBJECT} completed"
echo "#########################################################"

# filtering down to 1/4 - 5m streamlines
tcksift ${TCK_OUT} ${FOD_FILE} ${TCKSFT_OUT} -act ${SEG_FILE} -term_number 5m

echo "#########################################################"
echo "tcksift for ${SUBJECT} completed"
echo "#########################################################"

# create connectome matrix
# optional: -zero_diagonal and -symmetric
tck2connectome -symmetric -zero_diagonal ${TCKSFT_OUT} ${PARC_FILE} ${CONN_OUT}

echo "#########################################################"
echo "connectome for ${SUBJECT} completed"
echo "#########################################################"

trap : 0

done < $session_list
