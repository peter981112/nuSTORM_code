

export numu_ChargedPi="flat_E_source/numu/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_charged_pion_only.root"


export nue_ChargedPi="flat_E_source/nue/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_charged_pion_only.root"

export flat_E_flux_file_name="flat_E_test.txt"

export txt_prefix="nuSTORM_flux_txt_files/"
export txt_ext=".txt"
export root_ext=".root"
export slash="/"

#######################################################################################
mkexe.sh weight_rewrite_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

mkdir -p root_target
for ii in nuSTORM_flux_txt_files/*.txt
do

    echo flux txt file - $ii
    
    temp_txt_file=$(echo $ii | cut -c$((${#txt_prefix}+1))- | rev | cut -c$((${#txt_ext}+1))- | rev)
    echo temp_txt_file - $temp_txt_file
    for jj in flat_E_source/*/*.root
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
                    echo reweighting numu outana9
                    nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $jj > root_target/see_outana9_${temp_txt_file}_${temp_root_file}.log &
                fi
                if [[ "$jj" == *"outAna7"* ]] || [[ "$jj" == *"outana7"* ]]
                then
                    echo reweighting numu outana7
                    nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $jj > root_target/see_outana7_${temp_txt_file}_${temp_root_file}.log &
                fi
            fi
        fi

        if [[ "$ii" == *"Nue"* ]] ||[[ "$ii" == *"nue"* ]] ||[[ "$ii" == *"enu"* ]]
        then
            if [[ "$jj" == *"Nue"* ]] ||[[ "$jj" == *"nue"* ]] ||[[ "$jj" == *"enu"* ]]
            then
                if [[ "$jj" == *"outAna9"* ]] || [[ "$jj" == *"outana9"* ]]
                then
                    echo reweighting nue outana9
                    nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $jj > root_target/see_outana9_${temp_txt_file}_${temp_root_file}.log &
                fi
                if [[ "$jj" == *"outAna7"* ]] || [[ "$jj" == *"outana7"* ]]
                then
                    echo reweighting nue outana7
                    nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $jj > root_target/see_outana7_${temp_txt_file}_${temp_root_file}.log &
                fi
            fi
        fi
        
    done
done

echo done
