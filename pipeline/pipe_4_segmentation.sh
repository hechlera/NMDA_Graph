#!/bin/bash

### Andre Hechler, 29.05.2019
### Segmentation + registration of structural images to DWI

Usage() {
    cat <<EOF

Usage: pipe_3_segmentation-registration <session_list.csv>

  <input> is the session as csv file
  script needs exactly one input file

EOF
    exit 1
}

# set error message if code fails within "set -e"
abort() {
    echo "***Longitudinal registration was stopped due to an error. Please revise code.***"
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
#README_FILE=${ANAT_PATH}/README.txt

MPRAGE_FILE=$(find ${ANAT_PATH} -name "*T1w.nii.gz") 

echo MPRAGE_FILE ${MPRAGE_FILE}

# define output files for BET
MPRAGE_SEG="${MPRAGE_FILE%.nii.gz}_seg.nii.gz"
MPRAGE_GMWMI="${MPRAGE_FILE%.nii.gz}_seg-gmwmi.nii.gz"

set -e
#######################################

# Segmentation
5ttgen fsl ${MPRAGE_FILE} ${MPRAGE_SEG}

# gray matter-white matter interface mask
MPRAGE_SEG=$(find ${ANAT_PATH} -name "*T1w_seg.nii.gz") 
5tt2gmwmi ${MPRAGE_SEG} ${MPRAGE_GMWMI}


trap : 0

done < $session_list
