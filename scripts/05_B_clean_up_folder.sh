#!/bin/sh

# This script is useful to clean up the simulated reads folder after running 06_align_simulation.py
# Run it carefully as it was written very quickly

cd /g/typas/Personal_Folders/Neeka/Model/data/10_whole_pipeline/ERP105624/simulated_reads
for basefolder in $(ls)
do 
    cd $basefolder
    for file in $(find . | grep 01.fq$)
        do mv $file ${basefolder}_1.fq
        cd ..
    done
done

for basefolder in $(ls)
do 
    cd $basefolder
    for file in $(find . | grep 02.fq$)
        do mv $file ${basefolder}_2.fq
        cd ..
    done
done