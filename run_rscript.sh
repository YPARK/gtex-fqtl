#!/bin/bash -l

jobs=$1

source /broad/software/scripts/useuse > /dev/null
reuse -q R-3.3
reuse -q GCC-5.2
reuse -q .icc-2015
reuse -q Java-1.8

export MKROOT=/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/intel/mkl

export LD_LIBRARY_PATH=${LIBRARY_PATH}:${LD_LIBRARY_PATH}
export LD_PRELOAD=${MKLROOT}/lib/intel64/libmkl_core.so:${MKLROOT}/lib/intel64/libmkl_sequential.so

jobname=$(echo $jobs | awk '{ gsub("/","_"); print }')
command=$(zless $jobs | head -n $SGE_TASK_ID | tail -n 1)

mkdir -p log/rscript/${jobname}/
log=log/rscript/${jobname}/$SGE_TASK_ID.log
[ -f $log ] && rm -f $log

printf "[%s] Running ... \n$command" "$(date)" >> $log
printf "\n\n" >> $log

exec $command >> $log 2>&1

printf "[%s] \nDone\n\n" "$(date)" >> $log

