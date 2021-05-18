#INITIALISE FOLDERS
my_folder=/nfs/research1/marioni/alsu
out_folder=${my_folder}/clust_out/geneBasis
err_folder=${my_folder}/clust_err/geneBasis

#CHOOSE PARAMETERS
#RAM in megabytes
memory=50000

#num_processors
#nproc=4

#out_folder=/nfs/research1/marioni/alsu/geneBasis/data/processed_time/scmer/spleen

#cd $out_folder

n_samples_init=(1 2 3 4 5 6)
min_disp=(0 0.25 0.5)
n_threads=(1 3 6 9)
lambdas=(0.0002 0.0004 0.0006 0.0008 0.001 0.0015 0.002 0.0025 0.003)

for n_samples_discard_blood in 0 1 2 3 4
do
    for lambda in 0.0004 0.0006 0.0008 0.001 0.0015 0.002 0.0025
    do
        script_name=scmer_bat_eff
        bsub -e ${err_folder}/${script_name} \
        -o ${out_folder}/${script_name} \
        -M $memory -J ${script_name} \
        "python /nfs/research1/marioni/alsu/geneBasis/am_geneBasis/generateData/benchmark_method/batch_correction/scmer/batch_effect.py $n_samples_discard_blood $lambda"
    done
done
