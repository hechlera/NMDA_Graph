#!/bin/bash

### Andre Hechler, 29.05.2019
### Preprocessing of DWI images to tractogram construction

Usage() {
    cat <<EOF

Usage: pipe_2_dwipreproc.sh <session_list.csv>

  <input> is the session as csv file
  script needs exactly one input file

EOF
    exit 1
}

# set error message if code fails within "set -e"
abort() {
    echo "***DWI preprocessing was stopped due to an error. Please revise code.***"
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

# folders and files
ANAT_PATH=${SESSION_PATH}/anat
DWI_PATH=${SESSION_PATH}/dwi

MPRAGE_FILE=$(find ${ANAT_PATH} -name "*T1w.nii.gz") 
DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif") 
JSON_FILE=$(find ${DWI_PATH} -name "*.json") 

echo MPRAGE_FILE ${MPRAGE_FILE}
echo DWI_FILE ${DWI_FILE}

# define initial output file names
MPRAGE_SEGMENT="${MPRAGE_FILE%.nii.gz}-seg.nii.gz"
DWI_DENOISE="${DWI_FILE%.mif}_dn.mif"
DWI_NOISEMAP="${DWI_FILE%.mif}_noise.mif"


############# DWI Preprocessing (Tournier) #####################
## Run this script for (single shell, ) single tissue CSD ######
################################################################

# denoise
dwidenoise ${DWI_FILE} ${DWI_DENOISE} -noise ${DWI_NOISEMAP}

# preprocess
DWI_DN=$(find ${DWI_PATH} -name "*dwi_dn.mif") 
DWI_PREPROC="${DWI_DN%.mif}-preproc.mif"
dwipreproc ${DWI_DN} ${DWI_PREPROC} -rpe_header -json_import ${JSON_FILE}

# biascorrect with ANTs
DWI_PP=$(find ${DWI_PATH} -name "*dwi_dn-preproc.mif") 
DWI_BCOR="${DWI_PP%.mif}_-bcor.mif"
dwibiascorrect ${DWI_PP} ${DWI_BCOR} -ants

# Estimate response function
RFE="${DWI_FILE%.mif}_wm-rfe.txt"
dwi2response tournier ${DWI_BCOR} ${RFE}
###optional: check RFE - shview wm_response.txt


trap : 0

done < $session_list
