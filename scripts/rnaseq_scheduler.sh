#!/bin/bash

source ./scripts/pipeline_utilities.sh

### Log directory
LOCAL_DIR_LOGS="./tmp/logs"
if [ ! -d ${LOCAL_DIR_LOGS} ]; then mkdir -p ${LOCAL_DIR_LOGS}; fi
LOG="${LOCAL_DIR_LOGS}/captainslog_$( date +%Y%m%d%H%M%S ).txt"
# You need to define LOG as a global var. 

#
#
# Preliminaries and niceties 
#
#

# What is the id of your experiment? Give me something short and
# without weird characters, non-descriptive is great (ie dme211).
EXPERIMENT_NAME="testing"

# Are these Trumiseq adapters? Not supported
isTrumiSeq=0

# Module definitions - which ones are you using? 
# You probably don't want to touch these.
GFFREAD="gffread/intel/0.9.8c"
SAMTOOLS="samtools/intel/1.3.1"
HISAT2="hisat2/intel/2.0.5"
TRIM_GALORE="trim_galore/0.4.4"
CUTADAPT="cutadapt/intel/1.12"
#fastqc="fastqc/0.11.5"
R="r/intel/3.3.2"
QUALIMAP="qualimap/2.2.1"

log "You're running the RNAseq scheduler script on experiment 
\"${EXPERIMENT_NAME}\". I've made a log file at ${LOG} for you.
...
great"

#
#
# Moving files around #
#
#

# Then, we make local copies of everything. Why? 
# Reproducibility, via ASF1 ( Atomic Semantic Fusion 1 folder ).
# Local data directory, so in this folder you find the local copy
# of the data, configuration files, indicies, other crap.
LOCAL_DIR_DATA="./data"
if [ ! -d ${LOCAL_DIR_DATA} ] ; then mkdir -p ${LOCAL_DIR_DATA} ; fi
# For scratch files
LOCAL_DIR_TMP="./tmp"
if [ ! -d ${LOCAL_DIR_TMP} ] ; then mkdir -p ${LOCAL_DIR_TMP} ; fi
# For output
LOCAL_DIR_OUT="./out"
if [ ! -d ${LOCAL_DIR_OUT} ] ; then mkdir -p ${LOCAL_DIR_OUT} ; fi
# Input data, so this should aim at the place on the gencore scratch
INPUT_DIR_FASTQS='/scratch/dhm267/someTestDatar'
LOCAL_DIR_FASTQS="${LOCAL_DIR_DATA}/$( basename ${INPUT_DIR_FASTQS})"
cp -r ${INPUT_DIR_FASTQS} ${LOCAL_DIR_FASTQS}
# This should aim at the folder of HISAT2 indicies of your reference
REF_DIR_HISAT2INDEX="/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic_hisat_index"
LOCAL_DIR_HISAT2INDEX="${LOCAL_DIR_DATA}/$( basename $REF_DIR_HISAT2INDEX)"
LOCAL_HISAT2INDEX="${LOCAL_DIR_HISAT2INDEX}/genome"
cp -r ${REF_DIR_HISAT2INDEX} ${LOCAL_DIR_HISAT2INDEX}
# This should aim at an appropriate GFF for counting
REF_GFF="/genomics/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic.gff"
LOCAL_GFF="${LOCAL_DIR_DATA}/$( basename ${REF_GFF})"
cp -r ${REF_GFF} ${LOCAL_GFF}
# This is the name of the local GTF file. We make it in step 1
LOCAL_GTF=$( echo -n ${LOCAL_GFF} | sed 's/\.gff/.gtf/' )

log "Okay, I copied
  ${INPUT_DIR_FASTQS} to ${LOCAL_DIR_FASTQS} ,
  ${REF_DIR_HISAT2INDEX} to ${LOCAL_DIR_HISAT2INDEX} ,
  ${REF_GFF} to ${LOCAL_GFF} , 
  and I expect us to make a ${LOCAL_GTF}."

#
#
# We begin work.
#
#

##############################
# Generate gtf file from gff #
##############################

sbatchify "${LOCAL_DIR_LOGS}" \
  "GFFREAD" "" \
  "--nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00" \
  "module purge;\
   module load ${GFFREAD};\
   gffread -E ${LOCAL_GFF} -T -o- > ${LOCAL_GTF};\
  "
log "Started GFF to GTF conversion as jobid ${JOBID_GFFREAD}"

#####################
# Define the inputs #
#####################

# This is interpreting the ls command as an array
LOCAL_FASTQS=( $( ls $LOCAL_DIR_FASTQS/*fastq.gz ) )
# When we spit the array out and use David's perl one-liner to
# cut out just the sample ids, then look for unique ones.
LOCAL_SAMPLEIDS=( $( echo ${LOCAL_FASTQS[@]} | \
  perl -pe 's/[^\s]+n0[0-9]_(.+?)\.fastq\.gz/$1/g' | sort | uniq 
  ) )

log "Looking in ${LOCAL_DIR_FASTQS}, I see IDs:
  ${LOCAL_SAMPLEIDS[@]}
I'm going to start launching jobs for each of these. "

# Next, we iterate over each sample ID
for THIS_SAMPLE_ID in ${LOCAL_SAMPLEIDS[@]} ; do

  log "Scheduling for ID ${THIS_SAMPLE_ID}"

# Here, I'm only supporting reads that are 01 and 02.
# We find all the FASTQ files of that ID.
  THIS_SAMPLE_FILES=( $( \
    ls ${LOCAL_DIR_FASTQS}/*n0[12]_${THIS_SAMPLE_ID}.fastq.gz ) )

# And that let's us calculate if it's paired end or not.
  isPE=$(( ${#THIS_SAMPLE_FILES[@]} - 1 ))

# Then we report back what we found
  log "I see files ${THIS_SAMPLE_FILES[@]} with that sample ID."
  log "Since there are ${#THIS_SAMPLE_FILES[@]}, I conclude that
    the files are..."
  if [ "${isPE}" = 0 ] ; then
    log "   single ended reads."
  elif [ "${isPE}" = 1 ] ; then
    log "   paired ended reads."
  else
# incase we find something really weird
    log "... err I dunno. Giving up on this sample."
    (>&2 log "ERROR: Sample ${THIS_SAMPLE_ID} saw ${#THIS_SAMPLE_FILES[@]} 
    files, so gave up. Suggested action?")
    continue
  fi

# Now that we know what they are, we can slot them accordingly into
# variables.
  if [ "${isPE}" = 1 ] ; then
    THIS_SAMPLE_FASTQ_PAIR1=$( \
      ls ${LOCAL_DIR_FASTQS}/*01_${THIS_SAMPLE_ID}.fastq.gz )
    log "Paired end file 1 is : ${THIS_SAMPLE_FASTQ_PAIR1}"
    THIS_SAMPLE_FASTQ_PAIR2=$( \
      ls ${LOCAL_DIR_FASTQS}/*02_${THIS_SAMPLE_ID}.fastq.gz )
    log "Paired end file 2 is : ${THIS_SAMPLE_FASTQ_PAIR2}"
  elif [ "${isPE}" = 0 ] ; then
    THIS_SAMPLE_FASTQ=$( \
      ls ${LOCAL_DIR_FASTQS}/*01_${THIS_SAMPLE_ID}.fastq.gz )
    log "Single end file is : ${THIS_SAMPLE_FASTQ}"
  fi

####################
# Adapter trimming #
####################
# seqz from https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

  if [ "${isTrumiSeq}" = 1 ] ; then
    log "well good luck then"
    continue
# TODO
# This needs to be rewritten to take in an adapter sheet and a 
# sample sheet. Why? Because every Trumiseq sample has a different
# adapter. Darach has code for this. But it would require some 
# more work to generalize it...

  elif [ "${isTrumiSeq}" = 0 ] ; then

    if [ "${isPE}" = 1 ] ; then

      sbatchify "${LOCAL_DIR_LOGS}" \
        "TRIMMING" "" \
        "--nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:15:00"\
        "module purge; \
         module load ${CUTADAPT}; \
         cutadapt \
           -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
           -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
           -o ${LOCAL_DIR_TMP}/hisat2trimmed_${THIS_SAMPLE_ID}_1.fastq \
           -p ${LOCAL_DIR_TMP}/hisat2trimmed_${THIS_SAMPLE_ID}_2.fastq \
           ${THIS_SAMPLE_FASTQ_PAIR1} \
           ${THIS_SAMPLE_FASTQ_PAIR2} \
         ;"
      log "Started to trim paired end, jobid ${JOBID_TRIMMING}"

    elif [ "${isPE}" = 0 ] ; then

      sbatchify "${LOCAL_DIR_LOGS}" \
        "TRIMMING" "" \
        "--nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:15:00"\
        "module purge; \
         module load ${CUTADAPT}; \
         cutadapt \
           -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
           -o ${LOCAL_DIR_TMP}/hisat2trimmed_${THIS_SAMPLE_ID}.fastq \
           ${THIS_SAMPLE_FASTQ} \
         ;"
      log "Started to trim single end, jobid ${JOBID_TRIMMING}"

    fi

  fi

#########################
# Alignment with HISAT2 #
#########################

  if [ "${isPE}" = 1 ] ; then

    sbatchify "${LOCAL_DIR_LOGS}" \
      "HISAT2_ALIGN" "${JOBID_TRIMMING}" \
      "--nodes=1 --ntasks-per-node=8 --mem=30GB --time=00:30:00"\
      "module purge ; \
       module load ${HISAT2} ; \
       hisat2 -q -x ${LOCAL_HISAT2INDEX} \
         -1 ${LOCAL_DIR_TMP}/hisat2trimmed_${THIS_SAMPLE_ID}_1.fastq \
         -2 ${LOCAL_DIR_TMP}/hisat2trimmed_${THIS_SAMPLE_ID}_2.fastq \
         -S ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.sam \
       ;"
    log "Started aligning with HISAT2, jobid ${JOBID_HISAT2_ALIGN}"

  elif [ "${isPE}" = 0 ] ; then

    sbatchify "${LOCAL_DIR_LOGS}" \
      "HISAT2_ALIGN" "${JOBID_TRIMMING}" \
      "--nodes=1 --ntasks-per-node=8 --mem=30GB --time=00:30:00"\
      "module purge ; \
       module load ${HISAT2} ; \
       hisat2 -q -x ${LOCAL_HISAT2INDEX} \
         -U ${LOCAL_DIR_TMP}/hisat2trimmed_${THIS_SAMPLE_ID}.fastq \
         -S ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.sam \
       ;"
    log "Started aligning with HISAT2, jobid ${JOBID_HISAT2_ALIGN}"

  fi
    
  sbatchify "${LOCAL_DIR_LOGS}" \
    "HISAT2_CLEAN" "${JOBID_HISAT2_ALIGN}" \
    "--nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00"\
    "module purge; \
     module load ${SAMTOOLS} ; \
     samtools sort \
       ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.sam > \
       ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
     ;
     samtools index \
       ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
     ;"
  log "Started to clean up after HISAT2, using samtools,
    jobid ${JOBID_HISAT2_CLEAN}"
    
#
#	umi_tools dedup --umi-separator=: --method=directional --output-stats=${ID}.stats -I ${ID}.bam -S ${ID}.dedup.bam -L ${ID}.dedup.log
#

############################
# Generate table of counts #
############################

  if [ "${isPE}" = 1 ] ; then
 
    sbatchify "${LOCAL_DIR_LOGS}" \
      "RSUBREAD_QUANT" "${JOBID_HISAT2_CLEAN}" \
      "--nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00"\
      "module purge; \
       module load ${R} ; \
       Rscript ./scripts/count_exons_by_gene_id.R --args \
         ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
         ${LOCAL_DIR_OUT}/quantified_${THIS_SAMPLE_ID}.tsv \
         ${LOCAL_GTF} 0 1 \
       ;"
    log "Started to Rsubread count, jobid ${JOBID_RSUBREAD_QUANT}"

  elif [ "${isPE}" = 0 ] ; then

    sbatchify "${LOCAL_DIR_LOGS}" \
      "RSUBREAD_QUANT" "${JOBID_HISAT2_CLEAN}" \
      "--nodes=1 --ntasks-per-node=1 --mem=30GB --time=00:05:00"\
      "module purge; \
       module load ${R} ; \
       Rscript ./scripts/count_exons_by_gene_id.R --args \
         ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
         ${LOCAL_DIR_OUT}/quantified_${THIS_SAMPLE_ID}.tsv \
         ${LOCAL_GTF} 0 0 \
       ;"
    log "Started to Rsubread count, jobid ${JOBID_RSUBREAD_QUANT}"

  fi
    

##############
# QC metrics #
##############
#Run Qualimap on individual BAM file
#http://qualimap.bioinfo.cipf.es/doc_html/command_line.html#command-line
#NB. for qualimap rnaseq set -pe if paired end sequencing used 

  if [ "${isPE}" = 1 ] ; then

    sbatchify "${LOCAL_DIR_LOGS}" \
      "HISAT2_CLEAN" "${JOBID_HISAT2_ALIGN}" \
      "--nodes=1 --ntasks-per-node=8 --mem=30GB --time=00:05:00"\
      "module purge; \
       module load ${QUALIMAP} ; \
       qualimap bamqc \
         -bam ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
         -gff ${LOCAL_GFF} \
         -nt 8 \
         -outdir ${LOCAL_DIR_OUT} \
         -outfile ${LOCAL_DIR_OUT}/qualimap_bamReport_${THIS_SAMPLE_ID}.pdf \
       ;
       qualimap rnaseq \
         -bam ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
         -gtf ${LOCAL_GTF} \
         --paired -nt 8 \ 
         -outdir ${LOCAL_DIR_OUT} \
         -outfile ${LOCAL_DIR_OUT}/qualimap_rnaseqReport_${THIS_SAMPLE_ID}.pdf \
       ;"
    log "Started to clean up after HISAT2, using samtools,
      jobid ${JOBID_HISAT2_CLEAN}"
    
  elif [ "${isPE}" = 0 ] ; then
    
    sbatchify "${LOCAL_DIR_LOGS}" \
      "HISAT2_CLEAN" "${JOBID_HISAT2_ALIGN}" \
      "--nodes=1 --ntasks-per-node=8 --mem=30GB --time=00:05:00"\
      "module purge; \
       module load ${QUALIMAP} ; \
       qualimap bamqc \
         -bam ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
         -gff ${LOCAL_GFF} \
         -nt 8 \
         -outdir ${LOCAL_DIR_OUT} \
         -outfile ${LOCAL_DIR_OUT}/qualimap_bamReport_${THIS_SAMPLE_ID}.pdf \
       ;
       qualimap rnaseq \
         -bam ${LOCAL_DIR_TMP}/hisat2aligned_${THIS_SAMPLE_ID}.bam \
         -gtf ${LOCAL_GTF} \
         -nt 8 \ 
         -outdir ${LOCAL_DIR_OUT} \
         -outfile ${LOCAL_DIR_OUT}/qualimap_rnaseqReport_${THIS_SAMPLE_ID}.pdf \
       ;"
    log "Started to clean up after HISAT2, using samtools,
      jobid ${JOBID_HISAT2_CLEAN}"
    
  fi

  exit

done


#
#
### This is going to have to take custom adapter for using TrumiSeq
### Darach can write that part
##WORKDIR=$(pwd) 
##mkdir ${WORKDIR}"/tmp/dme211"
### data directory, obviously, pointing at directory with the fastqs demultiplexed
###   by the core
#
##module purge
##module load cutadapt/intel/1.12
##
### the below reads in an index, parses it into bash arrays, humor me okay?
##
##unset indicies
##declare -A indicies
##unset adapterName
##declare -A adapterName
##IFS=$'\n';
##for i in $(tail -n +2 data/dme211/dme211barcodeIndex.csv );do 
##  thisSampleName=$(echo -n $i | perl -ne '/^(.+?),(.+?),(.+?),(.+?)$/;print $1;' ); 
##  thisAdapterName=$(echo -n $i | perl -ne '/^(.+?),(.+?),(.+?),(.+?)$/;print $3;' ); 
##  adapterIndex=$(echo -n $i | perl -ne '/^(.+?),(.+?),(.+?),(.+?)$/;print $4;' ); 
##  indicies["${thisAdapterName}"]="${adapterIndex}";
##  adapterName["${thisSampleName}"]="${thisAdapterName}";
##done;
##echo "Read in:"
##echo ${!adapterName[@]}
##echo "as mapping to"
##echo ${adapterName[@]}
##echo "and"
##echo ${!indicies[@]}
##echo "as mapping to"
##echo "${indicies[@]}"
##echo 
##
##unset adaptSeq
##declare -A adaptSeq
##IFS=$'\n';
##for i in $(tail -n +2 data/dme211/trumiseqAdapters.csv );do 
##  thisSampleName=$(echo -n $i | perl -ne '/^(DG.+)_P7,(.+?)$/;print $1;' ); 
##  if [ "$thisSampleName" != "" ]; then
##    thisAdapterSeq=$(echo -n $i | perl -ne '/^(DGseq_.+?)_P7,(.+?)$/;print $2;' ); 
##    adaptSeq["$thisSampleName"]="${thisAdapterSeq}";
##  fi
##done;
##echo "Read in:"
##echo ${!adaptSeq[@]}
##echo "as mapping to"
##echo ${adaptSeq[@]}
##echo 
##
##for i in $(/bin/ls $DATADIR | grep "_[wq]1\?[0-9]_"); do
##  echo `date`
##  thisSampleName=$(echo -n $i | gawk -F _ '{print $3}');
##  thisAdapterName=${adapterName["$thisSampleName"]}
##  echo "doing file $i, which is $thisSampleName, which is $thisAdapterName"
##  thisAdaptSeq=${adaptSeq["$thisAdapterName"]}
##  thisIndex=${indicies["${adapterName["$thisSampleName"]}"]}
### below we grab the non-empty lines from the fastq, as there was a
### bioinformatics hiccup previously
##  runstring="
##cat ${DATADIR}$i | grep --regexp ^$ --invert-match > tmp/dme211/fixed$i;
##cutadapt -a A${thisAdaptSeq} --cut=0 -o tmp/dme211/dme211.${thisSampleName}.${adapterName["$thisSampleName"]}.$thisIndex.adapterTrimmed.fastq tmp/dme211/fixed$i"
##
#
#

#
#
#
#
#
## This archived for 
##SBATCH --array=0-23
##SBATCH --output=logs/%A_%a.o
##SBATCH --error=logs/%A_%a.e
##SBATCH --job-name=%a
##SBATCH --nodes=1
##SBATCH --tasks-per-node=1
##SBATCH --mem=20G
##SBATCH --time=24:00:00
##SBATCH --mail-user=${USER}@nyu.edu
##SBATCH --mail-type=FAIL
#



