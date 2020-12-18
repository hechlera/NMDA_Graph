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
DWI_PATH=${SESSION_PATH}/dwi
RFE_PATH=${PROJECT_DIRECTORY}/derivatives/response_functions

DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif")

echo DWI_FILE ${DWI_FILE}

###################### CSD (single tissue) ##################
## Performs CSD for (single shell,) single tissue DWI data ##
#############################################################

# Pool and average response function
RFE_FOLDER=${PROJECT_DIRECTORY}/derivatives/response_functions_tournier
#OUT_ARFE="avgwm-rfe.txt"
#cd ${RFE_FOLDER}
#average_response ./*.txt ${OUT_ARFE}

############################
# Optional: get and dilate brain mask (to ensure that the mask 
# doesn't cut into relevant tissue for streamline termination
### added maxdepth in case backup files lurk somewhere
DWI_BCOR=$(find ${DWI_PATH} -maxdepth 1 -name "*-bcor.mif")
DWI_MASK="${DWI_FILE%.mif}-preproc_mask.mif"
dwi2mask ${DWI_BCOR} ${DWI_MASK}

DWI_MASK=$(find ${DWI_PATH} -maxdepth 1 -name "*_mask.mif")
DWI_MASK_D="${DWI_MASK%.mif}-dil2.mif"
maskfilter ${DWI_MASK} dilate -npass 2 ${DWI_MASK_D}
#############################

# Perform constrained spherical deconvolution
RFE_FILE=$(find ${RFE_FOLDER} -name "*avgwm-rfe.txt") 
CSD="${DWI_FILE%.mif}_fod.mif"
dwi2fod csd ${DWI_BCOR} ${RFE_FILE} ${CSD} -mask ${DWI_MASK_D}
###optional: check fiber orientation density image - mrview fod.mif -odf.load_sh fod.mif



######################################## CSD (msmt) [not complete] ############
## In addition to HC's, This pipeline is also appropriate for PATIENTS,      ##
## especially with white matter damage that complicates tissue classification##
##                                                                           ##
## WARNING: To use this pipeline in MRTRIX, multi-shell data is necessary    ##
## (at least 3 b-values). Single shell, multi-tissue CSD is not yet implemen-##
## ted. Experimentally, the msmt algorithm can be run with "single" shell    ##
## data (b=0, b=1000), only supplying WM and CFS, but it is not citable yet  ##
###############################################################################

# Average response function
#average_response

# get and dilate brain mask (to ensure that the mask doesn't cut into 
# relevant tissue for streamline termination
#DWI_BC=$(find ${DWI_PATH} -name "*-bcor.mif")
#DWI_MASK="${DWI_FILE%.mif}-preproc_mask.mif"
#dwi2mask ${DWI_BCOR}

#DWI_MASK=$(find ${DWI_PATH} -name "*_mask.mif")
#DWI_MASK_D="${DWI_MASK%.mif}-dil2.mif"
#maskfilter ${DWI_MASK} dilate -npass 2 ${DWI_MASK_D}

# Perform constrained spherical deconvolution
#RFE_FILE=$(find ${DWI_PATH} -name "*wm-rfe.txt") 
#CSD="${DWI_FILE%.mif}_fod.mif"
#dwi2fod csd ${DWI_BCOR} ${RFE_FILE} ${CSD}
###optional: check fiber orientation density image - mrview fod.mif -odf.load_sh fod.mif


trap : 0

done < $session_list
