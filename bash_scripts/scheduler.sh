#!/bin/bash

source /broad/software/scripts/useuse
reuse UGER
reuse Python-3.6

cd ~/method_comparison
python3 job_scripts/initializer.py $@
