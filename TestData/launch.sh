#!/bin/bash

SETUPDIR=../SETUP/
DATADIR=../RUNDATA/
INPUT=testinput
OUTPUT=testoutput

export OMP_NUM_THREADS=8
export GPU_LIST="0 1 2"
cd $DATADIR
${SETUPDIR}nbody6.gpu < ${SETUPDIR}${INPUT} > ${SETUPDIR}${OUTPUT}
