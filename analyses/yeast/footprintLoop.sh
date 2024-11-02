#!/bin/bash
for i in 10 20 30 50 80
do
   sbatch runFootprinting.sh $i
done