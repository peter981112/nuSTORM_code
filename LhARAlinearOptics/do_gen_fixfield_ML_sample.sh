#run fixfield for each element
FFA_type_num=50 #$1
part_num=50
echo FFA_type_num  : ${FFA_type_num}

dir0=$(pwd)
data_dir=${dir0}/data

mkdir -p ${data_dir}/ML_input_sample ${data_dir}/ML_output_sample ${data_dir}/ML_FFA_sample
input_file_dir=${data_dir}/ML_input_sample
output_file_dir=${data_dir}/ML_output_sample

FFA_def_file_name="FFA_ML_test"
input_particle_file_name="ML_test_particle"
ML_input_particle_sample_file_name="ML_test_particle_input"

output_particle_file_name="ML_test_particle_out"

python3 fixfield_ML_sample_gen.py ${FFA_type_num} ${part_num} "data/ML_input_sample/"
echo particle generated

date
echo Starting fixfield ML sample gen simulation

for kk in $(seq 1 ${FFA_type_num})
do
    ((FFA_num=kk-1))
    echo ============================================================================================================================
    echo FFA num : ${FFA_num}
    latt_file="${FFA_def_file_name}-${FFA_num}.latt"

    particle_input="${input_particle_file_name}-${FFA_num}.dat"

    particle_output="${output_particle_file_name}-${FFA_num}.dat"

    echo start simulating latt file :  ${input_file_dir}/${latt_file}, particle input : ${input_file_dir}/${particle_input}, particle output : ${output_file_dir}/${particle_output}
    
    ./fixfield.exe ${input_file_dir}/${latt_file} ${input_file_dir}/${particle_input} ${output_file_dir}/${particle_output}

    echo ${kk}th FFA simulated
    date
done

echo simulation finished
date

ML_input_output_lst=()

for ii in ${input_file_dir}/${ML_input_particle_sample_file_name}*.txt
do
    ML_input_output_lst+=(${ii})
done

for ii in ${output_file_dir}/${output_particle_file_name}*.dat
do
    ML_input_output_lst+=(${ii})
done

#echo ML_input_output_lst : ${ML_input_output_lst[@]}

echo start merging for ML format

python3 fixfield_ML_merge.py ${ML_input_output_lst[@]}

cp ${data_dir}/*/*merged* ${data_dir}/ML_FFA_sample

echo formatting finished
date