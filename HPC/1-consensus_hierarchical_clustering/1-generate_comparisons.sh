for simul_study_id in $(seq 2 1 2)
do

nrows=$(expr $(cat ../Simulation_parameters/Simulation_parameters_list_$simul_study_id.txt | wc -l) - 1)

echo ID of simulation study: $simul_study_id

for j in $(seq 1 1 $nrows)
do
echo $j
sed "s/{simul_study_id_input}/${simul_study_id}/g" template_comparisons.sh > run1.sh
sed "s/{params_id_input}/${j}/g" run1.sh > run.sh
qsub run.sh
done
done

rm run1.sh
