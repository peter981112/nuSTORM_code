export GFSa9t1nuH=
export nopi="flat_E_source/outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"

export PIZERO="flat_E_source/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"

export flat_E_flux_file_name="flat_E_test.txt"

export flux_file_name="E7SpectraMuSig556Nue.txt"

export tag=MINERvALE_GiBUU_test_
#######################################################################################
mkexe.sh weight_rewrite_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

#exit

#void weight_rewrite(string flat_E_flux_file, string flux_file, TString fn)
#./weight_rewrite_arg "flat_E_test.txt" "E7SpectraMuSig556Nue.txt" "flat_E_source/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"
#./weight_rewrite_arg "flat_E_test.txt" "E7SpectraMuSig556Nue.txt" "flat_E_source/outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"


 #opt=${tag}GFSa9t1nuC;               
 nohup ./weight_rewrite_arg $flat_E_flux_file_name $flux_file_name $nopi  #> see${opt}.log &

 #opt=${tag}GFSPIZEROa7t4nuC;            
 nohup ./weight_rewrite_arg $flat_E_flux_file_name $flux_file_name $PIZERO   #> see${opt}.log &

exit