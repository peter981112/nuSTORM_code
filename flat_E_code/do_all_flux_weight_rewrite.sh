export numu_nopi="flat_E_source/numu/outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"

export numu_PIZERO="flat_E_source/numu/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"


export nue_nopi="flat_E_source/nue/outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"

export nue_PIZERO="flat_E_source/nue/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"

export flat_E_flux_file_name="flat_E_test.txt"

#export flux_file_name="E7SpectraMuSig556Nue.txt"
#######################################################################################
mkexe.sh weight_rewrite_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

for ii in nuSTORM_flux_txt_files/*.txt
do

    echo flux txt file - $ii

    if [[ "$ii" == *"Numu"* ]]
    then
        echo reweighting numu
        printf ""

        nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $numu_nopi  #> see${opt}.log &

 #opt=${tag}GFSPIZEROa7t4nuC;            
        nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $numu_PIZERO   #> see${opt}.log &
    fi

    if [[ "$ii" == *"Nue"* ]]
    then
        echo reweighting nue
        printf ""

        nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $nue_nopi  #> see${opt}.log &

 #opt=${tag}GFSPIZEROa7t4nuC;            
        nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $nue_PIZERO   #> see${opt}.log &
    fi
    #nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $nopi  #> see${opt}.log &

 #opt=${tag}GFSPIZEROa7t4nuC;            
    #nohup ./do_weight_rewrite.sh $flat_E_flux_file_name $ii $PIZERO   #> see${opt}.log &
    printf ""

done

echo done
