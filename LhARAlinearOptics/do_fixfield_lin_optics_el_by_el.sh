#run fixfield for each element
particle_type=$1
echo particle type : ${particle_type}

dir0=$(pwd)
data_dir=${dir0}/data

pion_init_input="custom_test_pion_bdsim.txt"
muon_init_input="custom_test_muon_bdsim.txt"
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

lattice_num=11

date
echo Starting fixfield element by element simulation
#./fixfield.exe data/nustorm_hyb_no_quad.latt data/custom_test_muon_bdsim.txt data/custum_test_out.dat

for kk in $(seq 1 ${lattice_num})
do
    latt_file=nustorm_single_FFA${kk}.latt

    #beam_input_file=${data_dir}/custum_test_out_bdsim_pions$((kk-1)).dat

    particle_output="custum_test_out_bdsim_${particle_type}${kk}.dat"

    echo start simulating latt file : ${latt_file}, particle input : ${particle_input}, particle output : ${particle_output}
    
    ./fixfield.exe ${data_dir}/${latt_file} ${data_dir}/${particle_input} ${data_dir}/${particle_output}

    particle_input=${particle_output}

    echo ${kk}th element simulated
    date
done

echo simulation finished
date