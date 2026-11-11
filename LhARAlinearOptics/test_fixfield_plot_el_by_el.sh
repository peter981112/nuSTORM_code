#run fixfield for each element
particle_type=$1
echo particle type : ${particle_type}

dir0=$(pwd)
data_dir=${dir0}/data

muon_out0="custum_test_out_bdsim_muon1.dat"
muon_out1="custum_test_out_bdsim_muon2.dat"
muon_out2="custum_test_out_bdsim_muon3.dat"
muon_out3="custum_test_out_bdsim_muon4.dat"
muon_out4="custum_test_out_bdsim_muon5.dat"

if [[ ${particle_type} == *"muon"* ]] || [[ ${particle_type} == *"Muon"* ]] || [[ ${particle_type} == *"mu"* ]]
then
    particle_input=${muon_init_input} #"custum_test_out_bdsim_pions0.dat"
    particle_type="muon"
fi

if [[ ${particle_type} == *"pion"* ]] || [[ ${particle_type} == *"Pion"* ]] || [[ ${particle_type} == *"pi"* ]]
then
    particle_input=${pion_init_input} #"custum_test_out_bdsim_pions0.dat"
    particle_type="pion"
fi

lattice_num=5

date
echo Start listing particle data files
#./fixfield.exe data/nustorm_hyb_no_quad.latt data/custom_test_muon_bdsim.txt data/custum_test_out.dat

muon_list=(${data_dir}/custom_test_muon_bdsim.txt)

muon_list+=(${data_dir}/${muon_out0})
muon_list+=(${data_dir}/${muon_out1})
muon_list+=(${data_dir}/${muon_out2})
muon_list+=(${data_dir}/${muon_out3})
muon_list+=(${data_dir}/${muon_out4})

echo muon_list : ${muon_list[@]}

python3 fixfield_simple_anal_ele.py ${muon_list[@]}

echo plotting finished
date