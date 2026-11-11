#run fixfield for each element
Beamline_type_num=6 #$1
part_num=5 #0
echo Beamline_type_num  : ${Beamline_type_num}

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

python3 fixfield_ML_sample_gen_v3.py ${Beamline_type_num} ${part_num} "data/ML_input_sample/"
echo particle generated

date
echo Starting fixfield ML sample gen simulation


for kk in $(seq 1 ${Beamline_type_num})
do
    ((Beamline_num=kk-1))
    echo ============================================================================================================================
    echo Beamline num : ${Beamline_num} #FFA_index : ${FFA_index}
    particle_input="${input_particle_file_name}-${Beamline_num}.dat"

    FFA_num_in_beamline=8
    FFA_index=0
    latt_file_lst=()
    particle_output_lst=()
    particle_input_name=data/custum_test_out_bdsim_muons.dat #${input_file_dir}/${particle_input}
    for jj in $(seq 1 ${FFA_num_in_beamline})
    do
        latt_file_name="${input_file_dir}/${FFA_def_file_name}-${FFA_index}-${Beamline_num}.latt"
        latt_file_lst+=(${latt_file_name})
        particle_output_name=" ${output_file_dir}/${output_particle_file_name}-${FFA_index}-${Beamline_num}.dat"
        particle_output_lst+=(${particle_output_name})
        echo ./fixfield.exe ${latt_file_name} ${particle_input_name} ${particle_output_name}
        ./fixfield.exe ${latt_file_name} ${particle_input_name} ${particle_output_name}
        ((FFA_index=FFA_index+1))
        particle_input_name=${particle_output_name}
    done

    #echo start simulating latt file :  ${input_file_dir}/${latt_file}, particle input : ${input_file_dir}/${particle_input}, particle output : ${output_file_dir}/${particle_output}

    echo ${kk}th FFA simulated
    date
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

python3 fixfield_ML_merge_v3.py ${ML_input_output_lst[@]}

cp ${data_dir}/*/*merged* ${data_dir}/ML_FFA_sample

rm ${data_dir}/*input*/*merged* ${data_dir}/*output*/*merged*

echo fixedfield data formatting for ML finished
date