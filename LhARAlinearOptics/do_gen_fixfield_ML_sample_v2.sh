#run fixfield for each element
FFA_type_num_1=5 #$1
FFA_type_num_2=5 #$1
part_num=5 #0
echo FFA_type_num_1  : ${FFA_type_num_1}   FFA_type_num_2  : ${FFA_type_num_2}

dir0=$(pwd)
data_dir=${dir0}/data

mkdir -p ${data_dir}/ML_input_sample ${data_dir}/ML_output_sample ${data_dir}/ML_FFA_sample
input_file_dir=${data_dir}/ML_input_sample
output_file_dir=${data_dir}/ML_output_sample

FFA_def_file_name="FFA_ML_test"
FFA_input_data_file_name="ML_test_FFA_input"
input_particle_file_name="ML_test_particle_input"
ML_input_particle_sample_file_name="ML_test_particle_input"
output_particle_file_name="ML_test_particle_out"

python3 fixfield_ML_sample_gen_v2.py ${FFA_type_num_1} ${FFA_type_num_2} ${part_num} "data/ML_input_sample/"
echo particle generated

date
echo Starting fixfield ML sample gen simulation

FFA_index=0
for kk in $(seq 1 ${FFA_type_num_1})
do
    for jj in $(seq 1 ${FFA_type_num_2})
    do
        ((FFA_num_1=kk-1))
        ((FFA_num_2=jj-1))
        echo ============================================================================================================================
        echo FFA num 1 : ${FFA_num_1} FFA num 2 : ${FFA_num_2} FFA_index : ${FFA_index}
        latt_file_1="${FFA_def_file_name}-1-${FFA_index}.latt"
        latt_file_2="${FFA_def_file_name}-2-${FFA_index}.latt"

        particle_input="${input_particle_file_name}-${FFA_index}.dat"

        particle_output_1="${output_particle_file_name}-1-${FFA_index}.dat"

        particle_output_2="${output_particle_file_name}-2-${FFA_index}.dat"

        echo start simulating latt file :  ${input_file_dir}/${latt_file}, particle input : ${input_file_dir}/${particle_input}, particle output : ${output_file_dir}/${particle_output}
        
        ./fixfield.exe ${input_file_dir}/${latt_file_1} ${input_file_dir}/${particle_input} ${output_file_dir}/${particle_output_1}
        ./fixfield.exe ${input_file_dir}/${latt_file_2} ${output_file_dir}/${particle_output_1} ${output_file_dir}/${particle_output_2}

        ((FFA_index=FFA_index+1))

        echo ${kk}th FFA simulated
        date
    done
done

echo simulation finished
date

ML_input_output_lst=()

for ii in ${input_file_dir}/${input_particle_file_name}*.dat
do
    ML_input_output_lst+=(${ii})
done

for ii in ${output_file_dir}/${output_particle_file_name}*.dat
do
    ML_input_output_lst+=(${ii})
done

for ii in ${input_file_dir}/${FFA_input_data_file_name}*.txt
do
    ML_input_output_lst+=(${ii})
done



#echo ML_input_output_lst : ${ML_input_output_lst[@]}

echo start merging for ML format

python3 fixfield_ML_merge_v2.py ${ML_input_output_lst[@]}

cp ${data_dir}/*/*merged* ${data_dir}/ML_FFA_sample

rm ${data_dir}/*input*/*merged* ${data_dir}/*output*/*merged*

echo fixedfield data formatting for ML finished
date