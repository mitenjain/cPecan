#!/bin/bash

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

#Get directory containing shell script
BASEDIR=$(dirname $0)/../

#Set the python path to just the relevant directories
export PYTHONPATH=${BASEDIR}

#Set the path to just the relevant directories
export PATH=${BASEDIR}/sonLib/bin:${BASEDIR}/jobTree/bin:${PATH}

#Run described script
$@

