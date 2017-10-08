#!/bin/bash

# This will echo what you give it, and write it to a variable named
# LOG in your current environment. This is done, not an argument, to
# make it super easy to write it in and more readable.
log(){
  echo "\n $1 \n" | tee -a ${LOG}
}

# This will submit your string as an sbatch, with dependencies,
# resource specifications, log it to the right place, and will save
# the job id back as `JOBID_${STEP_NAME}`
sbatchify(){
  export LOCAL_DIR_LOGS=$1
  export STEP_NAME=$2
  export DEPENDS_ON=$3
  export STEP_RESOURCES=$4
  export STEP_STRING=$5
#
  if [ -z "${DEPENDS_ON}" ] ; then 
    RUNSTRING="
      sbatch
        --mail-type=BEGIN,END,FAIL --mail-user=${USER}@nyu.edu
        --job-name=${STEP_NAME}
        --output=${LOCAL_DIR_LOGS}/%A_${STEP_NAME}.out 
        --error=${LOCAL_DIR_LOGS}/%A_${STEP_NAME}.err 
        ${STEP_RESOURCES}
        --wrap=\"${STEP_STRING}\"
    ;"
  else 
    RUNSTRING="
      sbatch
        --mail-type=BEGIN,END,FAIL --mail-user=${USER}@nyu.edu
        --job-name=${STEP_NAME}
        --dependency=afterok:${DEPENDS_ON}
        --output=${LOCAL_DIR_LOGS}/%A_${STEP_NAME}.out 
        --error=${LOCAL_DIR_LOGS}/%A_${STEP_NAME}.err 
        ${STEP_RESOURCES}
        --wrap=\"${STEP_STRING}\"
    ;"
  fi
#
  log "Trying to run:
${RUNSTRING}"
#
  RESPONSE=( $(eval ${RUNSTRING} ) )
#
  eval JOBID_${STEP_NAME}="'${RESPONSE[-1]}'"
#
}

