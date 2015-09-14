#!/bin/bash

USAGE='USAGE: ./run.sh <param_file> [<extra_args>]'

if [ -z $1 ]
  then
    echo $USAGE
    exit
fi

export CADO_HOME=/usr/local/src/cado-nfs
export PYTHONPATH=$CADO_HOME/scripts/cadofactor/

PARAMS_FILE=$1
shift
python3 factor.py $PARAMS_FILE $@
