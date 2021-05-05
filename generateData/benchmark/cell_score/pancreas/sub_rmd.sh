#INITIALISE FOLDERS
my_folder=/nfs/research1/marioni/alsu
out_folder=${my_folder}/clust_out/geneBasis
err_folder=${my_folder}/clust_err/geneBasis

#SELECT SCRIPT
#If you change this, you MUST update the wrapper's grep
script_name=cs_pan
#CHOOSE PARAMETERS
#RAM in megabytes
memory=300000
r_command="rusage[mem=${memory}]"
#num_processors
nproc=5

smg=/nfs/research1/marioni/alsu/singularity/R1.simg
script=/nfs/research1/marioni/alsu/geneBasis/am_geneBasis/generateData/benchmark/cell_score/pancreas/run_rmd.R

bsub -q research-rh74 -e ${err_folder}/${script_name} \
-o ${out_folder}/${script_name} \
-M $memory -R $r_command -n $nproc -P bigmem -J ${script_name} \
"singularity exec $smg Rscript $script"
