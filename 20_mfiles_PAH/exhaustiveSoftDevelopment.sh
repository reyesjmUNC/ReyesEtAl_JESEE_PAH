##!/bin/bash
for ((i=1;i<=11;i++))
do
  bsub -x -q day -n 12 -R "span[hosts=1]" /nas02/apps/matlab-2012b/matlab -nodesktop -nosplash -singleCompThread -r "exhaustiveSoftDevelopment(${i})" -logfile "exhaustiveSoftDevelopment${i}.out"
  sleep 30
done