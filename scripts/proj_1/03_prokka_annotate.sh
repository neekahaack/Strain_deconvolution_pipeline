#micromamba activate prokka_env

#should maybe navigate to directory where unzipped fastas are? uneccessary
#cd /g/typas/Personal_Folders/Neeka/Model/data/01_sequences/proj_name/proj_name_unzipped

#need to make a "for do" loop to run prokka on all files
for file in /g/typas/Personal_Folders/Neeka/Model/data/01_sequences/PRJNA903559/PRJNA903559_unzipped/*.fna
    do 
    base=$(basename "$file" .fna)
    prokka \
    --cpus 32 \
    --outdir /g/typas/Personal_Folders/Neeka/Model/data/03_annotation/PRJNA903559/$base \
    --prefix $base \
    "$file"
    done

#accidentally didnt direct prokka to use all 32 cpu, so session might end before all genomes are annotated
#might need to figure out how to direct it to start where it left off



#notes
#in bash, for info on function, type help "function"
#‘basename’: Strip directory and suffix from a file name
#"$" is command substitution, runs function and captures output as string
#echo is print for bash !!! 
#no " " around "=" in bash
#$file → expands the variable and then the shell splits it on spaces and special characters.
#"$file" → expands the variable but keeps it as one single argument, even if it contains spaces, tabs, or special characters.

file="/g/typas/Personal_folders/Neeka/Model/data/test_one_genome/annotation/second try.txt"
base=$(basename $file .txt)
echo $base
echo "$file"
basename "$file" .txt
echo $file
basename $file .txt #doesnt work
echo file
echo "file"