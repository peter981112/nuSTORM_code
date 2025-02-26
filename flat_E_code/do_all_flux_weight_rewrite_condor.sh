
flat_E_flux_file_name="flat_E_test.txt"
curr_path=$(realpath $(pwd))
echo curr paht - $curr_path
txt_prefix="$curr_path/nuSTORM_flux_txt_files/"
txt_ext=".txt"
root_ext=".root"
slash="/"

#export NUGENTKI=$(pwd)
#export LOCALBIN=${NUGENTKI}/bin
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NUGENTKI}/style
#export PATH=$PATH:${LOCALBIN}

#######################################################################################
mkexe.sh weight_rewrite_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

mkdir -p root_target

rm flux_weight_rewrite_condor_arg.txt #write argumets to this txt file and submited to condor. delete txt before writing new one

for ii in $curr_path/nuSTORM_flux_txt_files/*.txt
do

    echo flux txt file - $ii
    
    temp_txt_file=$(echo $ii | cut -c$((${#txt_prefix}+1))- | rev | cut -c$((${#txt_ext}+1))- | rev)
    echo temp_txt_file - $temp_txt_file
    for jj in $curr_path/flat_E_source/*/*.root
    do
        temp_root_file=$(basename ${jj})
        temp_root_file=$(echo $temp_root_file | cut -c$((${#slash}))- | rev | cut -c$((${#root_ext}+1))- | rev)
        echo temp_root_file - $temp_root_file

        if [[ "$ii" == *"Numu"* ]] || [[ "$ii" == *"numu"* ]] || [[ "$ii" == *"munu"* ]]
        then
            if [[ "$jj" == *"Numu"* ]] || [[ "$jj" == *"numu"* ]] || [[ "$jj" == *"munu"* ]]
            then

                if [[ "$jj" == *"outAna9"* ]] || [[ "$jj" == *"outana9"* ]]
                then
                    printf "$flat_E_flux_file_name $ii $jj $curr_path/root_target/see_outana9_${temp_txt_file}_${temp_root_file}\n" >> flux_weight_rewrite_condor_arg.txt
                fi
                if [[ "$jj" == *"outAna7"* ]] || [[ "$jj" == *"outana7"* ]]
                then
                    printf "$flat_E_flux_file_name $ii $jj $curr_path/root_target/see_outana7_${temp_txt_file}_${temp_root_file}\n"  >> flux_weight_rewrite_condor_arg.txt
                fi
            fi
        fi

        if [[ "$ii" == *"Nue"* ]] ||[[ "$ii" == *"nue"* ]] ||[[ "$ii" == *"enu"* ]]
        then
            if [[ "$jj" == *"Nue"* ]] ||[[ "$jj" == *"nue"* ]] ||[[ "$jj" == *"enu"* ]]
            then
                if [[ "$jj" == *"outAna9"* ]] || [[ "$jj" == *"outana9"* ]]
                then
                    printf "$flat_E_flux_file_name $ii $jj $curr_path/root_target/see_outana9_${temp_txt_file}_${temp_root_file}\n" >> flux_weight_rewrite_condor_arg.txt
                fi
                if [[ "$jj" == *"outAna7"* ]] || [[ "$jj" == *"outana7"* ]]
                then
                    printf "$flat_E_flux_file_name $ii $jj $curr_path/root_target/see_outana7_${temp_txt_file}_${temp_root_file}\n" >> flux_weight_rewrite_condor_arg.txt
                fi
            fi
        fi
        
    done
done

condor_submit flux_weight_rewrite_sub.sub #submit to condor with the sub file

echo done
