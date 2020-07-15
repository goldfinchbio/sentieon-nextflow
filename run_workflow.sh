#!/bin/bash

SAMPLE=$1
FASTQ1=$2
FASTQ2=$3
S3OUTPUT=$4
S3WORKDIR=$5
S3LOGDIR=$6
WES=$7
AWSREGION=$8
BATCH_JOBDEF=$9

CODE_VER=$(</usr/src/app/code_vers)
BASEOUTPUT=s3://${S3OUTPUT}/
WORKDIR=s3://${S3WORKDIR}/
WORKFLOWLOGPATH=s3://${S3LOGDIR}/

nextflow run sentieon_dnaseq.nf \
  -with-report ${SAMPLE}_report.html \
  -with-trace ${SAMPLE}_trace.txt \
  -with-timeline ${SAMPLE}_timeline.html \
  -with-dag ${SAMPLE}_flowchart.png \
  --fastqr1 ${FASTQ1} \
  --fastqr2 ${FASTQ2} \
  --sample_id ${SAMPLE} \
  --outpath ${BASEOUTPUT} \
  --code_ver ${CODE_VER} \
  -bucket-dir ${WORKDIR}/${SAMPLE} \
  --output ${BASEOUTPUT}/${SAMPLE} \
  --jobdef ${BATCH_JOBDEF} \
  --wes ${WES}

retVal=$?

aws s3 cp ${SAMPLE}_report.html ${WORKFLOWLOGPATH}/${SAMPLE}/ --region ${AWSREGION}
aws s3 cp ${SAMPLE}_trace.txt ${WORKFLOWLOGPATH}/${SAMPLE}/ --region ${AWSREGION}
aws s3 cp ${SAMPLE}_timeline.html ${WORKFLOWLOGPATH}/${SAMPLE}/ --region ${AWSREGION}
aws s3 cp ${SAMPLE}_flowchart.png ${WORKFLOWLOGPATH}/${SAMPLE}/ --region ${AWSREGION}
aws s3 cp .nextflow.log ${WORKFLOWLOGPATH}/${SAMPLE}/ --region ${AWSREGION}

# Return Code for Batch Job Status
if [ $retVal -ne 0 ]; then
    echo "The nextflow process has failed."
    exit $retVal
else
    echo "The nextflow process is successful."
fi
