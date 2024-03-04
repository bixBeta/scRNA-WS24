#!/bin/bash

export TZ=America/New_York
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

echo "start jupyter \$1=$1"
notebook_dir=$1
if [[ -z $notebook_dir ]]; then
  notebook_dir=/workdir
fi
preferred_dir=$2
if [[ -z $preferred_dir ]]; then
   preferred_dir=$notebook_dir
fi

echo "calling jupyter-lab --notebook-dir=$notebook_dir --preferred-dir=$preferred_dir --allow-root --no-browser --ip=0.0.0.0"
jupyter-lab --notebook-dir=$notebook_dir --preferred-dir=$preferred_dir --allow-root --no-browser --ip=0.0.0.0 


