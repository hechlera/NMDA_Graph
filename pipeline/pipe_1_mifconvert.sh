#!/bin/bash

[ "$#" -ne 1 ] && Usage
#Checks if number of inputs is unequal to 1 (one csv file). If this evaluates to #true the usage message above is returned

#####################################################
#### Inputs #########################################
#####################################################

###ATTENTION: neurodial_get_project_parameters.sh has to be set to the correct path!!!###

session_list=$1


while IFS=',' read SUBJECT SESSION; do
    PROJECT_DIRECTORY=/home/andre/Desktop/Thesis/01_data

if [ ! -d ${PROJECT_DIRECTORY} ]; then
	echo -e "\n\n###ERROR: Project directory ${PROJECT_DIRECTORY} does not exist. Please check call within neurodial_get_project_parameters.sh###\n\n"
	break
fi

# Get folderpath
INPUT_DWI=${PROJECT_DIRECTORY}/derivatives/sub-${SUBJECT}/ses-${SESSION}/dwi
#

if [ ! -d ${INPUT_DICOM} ]; then
	echo -e "\n\n###ERROR: DICOM files at ${INPUT_DICOM} not found. Please check path###\n\n"
	continue
fi

cd ${INPUT_DWI}

DWI_FILE=$(find ${INPUT_DWI} -name "*dwi.nii.gz") 
BVEC_FILE=$(find ${INPUT_DWI} -name "*dwi.bvec") 
BVAL_FILE=$(find ${INPUT_DWI} -name "*dwi.bval")

mrconvert ${DWI_FILE} sub-${SUBJECT}_ses-${SESSION}_dwi.mif -fslgrad ${BVEC_FILE} ${BVAL_FILE}

done < ${session_list}
