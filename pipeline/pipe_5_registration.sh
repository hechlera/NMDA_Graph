#!/bin/bash

### Andre Hechler, 29.05.2019
### registration of structural images to DWI

Usage() {
    cat <<EOF

Usage: pipe_5_registration <session_list.csv>

  <input> is the session as csv file
  script needs exactly one input file

EOF
    exit 1
}

# set error message if code fails within "set -e"
abort() {
    echo "***Registration was stopped due to an error. Please revise code.***"
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
FREES_PATH=/home/andre/Desktop/Thesis/03_data_freesurfer/sub-${SUBJECT}_ses-${SESSION}re/mri

DWI_FILE=$(find ${DWI_PATH} -name "*dwi.nii.gz")
GMWMI_FILE=$(find ${ANAT_PATH} -name "*gmwmi.nii.gz") 

PARC_FILE=$(find ${FREES_PATH} -name "aparc+aseg.mgz")
T1_MASK=$(find ${FREES_PATH} -name "brainmask.mgz")

# define output files names
DWI_MASK="${DWI_FILE%.nii.gz}_mask.nii.gz"
DWI_3D="${DWI_MASK%.nii.gz}-3D.nii.gz"

set -e
#######################################

# Create folder for transformation matrices
echo "Creating output folder for registration matrices"
mkdir ${ANAT_PATH}/registration_files
REG_FOLDER=${ANAT_PATH}/registration_files

# Copy freesurfer files
cp ${PARC_FILE} ${ANAT_PATH}/
cp ${T1_MASK} ${ANAT_PATH}/

# Convert everything to nii.gz
PARC_FILE=$(find ${ANAT_PATH} -name "aparc+aseg.mgz")
T1_MASK=$(find ${ANAT_PATH} -name "brainmask.mgz")
mrconvert ${PARC_FILE} "${PARC_FILE%.mgz}.nii.gz"
mrconvert ${T1_MASK} "${T1_MASK%.mgz}.nii.gz"

# BET on DWI file
## BET isn't really successful in skullstripping DWI files. Probably not optimised for 
## the low intensities of DWI skull signal. 
bet ${DWI_FILE} ${DWI_MASK} -R

# Get 3D version of 4D DWI file
mrmath ${DWI_MASK} mean -axis 3 ${DWI_3D}

# Get transform matrix from flirt based on freesurfer brainmask and 3D DWI mask
T1_MASK=$(find ${ANAT_PATH} -name "brainmask.nii.gz")
flirt -in ${T1_MASK} -ref ${DWI_3D} -omat ${REG_FOLDER}/T12DWI.mat -dof 6

# Convert flirt matrix to mrtrix matrix
transformconvert ${REG_FOLDER}/T12DWI.mat ${T1_MASK} ${DWI_3D}  flirt_import ${REG_FOLDER}/T12DWI.mrtrix

# Apply transform matrix to parcellated image and GMWMI file
PARC_FILE=$(find ${ANAT_PATH} -name "aparc+aseg.nii.gz")
mrtransform ${PARC_FILE} -linear ${REG_FOLDER}/T12DWI.mrtrix -interp nearest "${PARC_FILE%.nii.gz}_regDWI.nii.gz"
mrtransform ${GMWMI_FILE} -linear ${REG_FOLDER}/T12DWI.mrtrix -interp nearest "${GMWMI_FILE%.nii.gz}_regDWI.nii.gz"


trap : 0

done < $session_list
 
