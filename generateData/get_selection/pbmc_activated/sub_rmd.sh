#!/bin/sh
#INITIALISE FOLDERS
my_folder=/nfs/research/marioni/alsu
out_folder=${my_folder}/clust_out/geneBasis
err_folder=${my_folder}/clust_err/geneBasis

#SELECT SCRIPT
#If you change this, you MUST update the wrapper's grep
script_name=genes_pbmc

#CHOOSE PARAMETERS
#RAM in megabytes
memory=500000
r_command="rusage[mem=${memory}]"
#num_processors
nproc=4

smg=/hps/software/users/marioni/alsu/singularity/alsu_image.simg
script=/nfs/research/marioni/alsu/geneBasis/am_geneBasis/generateData/get_selection/pbmc_activated/run_rmd.R

bsub -q bigmem -e ${err_folder}/${script_name} \
-o ${out_folder}/${script_name} \
-M $memory -R $r_command -n $nproc -J ${script_name} \
"singularity exec $smg Rscript $script"
