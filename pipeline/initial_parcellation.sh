#!/bin/bash

[ "$#" -ne 1 ] && Usage

session_list=$1

while IFS=',' read SUBJECT SESSION; do
	
	#Define paths to anatomical files
	PROJECT_DIRECTORY=/home/andre/Desktop/Thesis/01_data/

	DERIVATIVES_FOLDER=${PROJECT_DIRECTORY}/derivatives
	ANAT_FOLDER=${DERIVATIVES_FOLDER}/sub-${SUBJECT}/ses-${SESSION}/anat
	MPRAGE_FILE=$(find ${ANAT_FOLDER} -name "*T1w.nii.gz") 

	if [ -z ${MPRAGE_FILE} ]; then
		echo "No MPRAGE file found for subject ${SUBJECT} session ${SESSION}. Continuing with next subject..."
		continue
	fi
	
	SUBID=sub-${SUBJECT}_ses-${SESSION}

	### for testing
	#echo "sub ID is ${SUBID}"
	#echo "MPRAGE at ${MPRAGE_FILE}"

	recon-all -all -s ${SUBID} -i ${MPRAGE_FILE}

done < $session_list
