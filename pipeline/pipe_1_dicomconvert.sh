#!/bin/bash

### Andre Hechler, 14.11.2018
### based on dcm2nii

Usage() {
    cat <<EOF

Usage: pipe_1_dicomconvert <session_list.csv>

  <input> is the session as csv file
  script needs exactly one input file

EOF
    exit 1
}

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
INPUT_DICOM=${PROJECT_DIRECTORY}/dicoms/NMDA/${SUBJECT}/${SUBJECT}_${SESSION}
#

if [ ! -d ${INPUT_DICOM} ]; then
	echo -e "\n\n###ERROR: DICOM files at ${INPUT_DICOM} not found. Please check path###\n\n"
	continue
fi

BIDS_CONFIG="/home/andre/Desktop/Thesis/04_scripts/bids_config.json"

if [ ! -e ${BIDS_CONFIG} ]; then
	echo -e "\n\n###ERROR: BIDS config file at ${BIDS_CONFIG} does not exist. Please check path###\n\n"
	break
fi

SOURCE_FOLDER=${PROJECT_DIRECTORY}/sourcedata
DERIVATIVES_FOLDER=${PROJECT_DIRECTORY}/derivatives

    #echo dicoms at ${INPUT_DICOM}

# Run conversion
dcm2bids -d ${INPUT_DICOM} -p ${SUBJECT} -s ${SESSION} -c ${BIDS_CONFIG} -o ${SOURCE_FOLDER}

# Copy files to derivates folder for processing and analysis
cp -r ${SOURCE_FOLDER}/sub-${SUBJECT} ${DERIVATIVES_FOLDER}
# Create readme, write conversion time
touch ${SOURCE_FOLDER}/sub-${SUBJECT}/ses-${SESSION}/README.txt
echo "$(date): Files converted" >> ${SOURCE_FOLDER}/sub-${SUBJECT}/ses-${SESSION}/README.txt
touch ${DERIVATIVES_FOLDER}/sub-${SUBJECT}/ses-${SESSION}/README.txt
echo "$(date): Files converted" >> ${DERIVATIVES_FOLDER}/sub-${SUBJECT}/ses-${SESSION}/README.txt

done < ${session_list}
