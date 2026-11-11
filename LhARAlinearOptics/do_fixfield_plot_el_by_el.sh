#run fixfield for each element
particle_type=$1
echo particle type : ${particle_type}

dir0=$(pwd)
data_dir=${dir0}/data

if [[ ${particle_type} == *"muon"* ]] || [[ ${particle_type} == *"Muon"* ]] || [[ ${particle_type} == *"mu"* ]]
then
    particle_input=${muon_init_input} #"custum_test_out_bdsim_pions0.dat"
    particle_type="muon"
    particle_list=(${data_dir}/custom_test_muon_bdsim.txt)
fi

if [[ ${particle_type} == *"pion"* ]] || [[ ${particle_type} == *"Pion"* ]] || [[ ${particle_type} == *"pi"* ]]
then
    particle_input=${pion_init_input} #"custum_test_out_bdsim_pions0.dat"
    particle_type="pion"
    particle_list=(${data_dir}/custom_test_pion_bdsim.txt)
fi

date
echo Start listing particle data files

for ii in $data_dir/custum_test_out_bdsim_*.dat
do
#if data is the given particle type and #if data doesn't have any numbers in file name string
    if [[ $ii == *$particle_type* ]] && [[ $ii == *[0123456789]* ]]
    then
        particle_list+=(${ii})
    fi
    

done

echo particle_list : ${particle_list[@]}

echo start plotting

python3 fixfield_simple_anal_ele.py ${particle_list[@]}

echo plotting finished
date