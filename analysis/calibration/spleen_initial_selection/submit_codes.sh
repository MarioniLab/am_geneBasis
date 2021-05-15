#INITIALISE FOLDERS
my_folder=/nfs/research1/marioni/alsu
out_folder=${my_folder}/clust_out/geneBasis
err_folder=${my_folder}/clust_err/geneBasis

#CHOOSE PARAMETERS
#RAM in megabytes
memory=45000

#num_processors
#nproc=4

for num_markers in 25 50 100
do
    for eps in 0.01 0.1 0.5 1 2 5 10 100
    do
        script_name=scgene_cal_spl
        bsub -e ${err_folder}/${script_name} \
        -o ${out_folder}/${script_name} \
        -M $memory -J ${script_name} \
        "python /nfs/research1/marioni/alsu/geneBasis/am_geneBasis/analysis/calibration/spleen_initial_selection/spleen_scGeneFit_options.py $num_markers $eps"
    done
done
