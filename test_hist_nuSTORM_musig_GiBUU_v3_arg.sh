export E3_outAna9="E3SpectraMuSig560Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E4_outAna9="E4SpectraMuSig566Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E5_outAna9="E5SpectraMuSig558Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E6_outAna9="E6spectraMuSig557Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E7_outAna9="E7SpectraMuSig556Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"

export E3_outAna7="E3SpectraMuSig560Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E4_outAna7="E4SpectraMuSig566Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E5_outAna7="E5SpectraMuSig558Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E6_outAna7="E6spectraMuSig557Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E7_outAna7="E7SpectraMuSig556Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"



#export root_prefix=${nopi%/*}/
#echo root_prefix $root_prefix

#export root_ext=".root"

#######################################################################################
mkexe.sh hist_nuSTORM_musig_GiBUU_v3_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

#temp_root_file=$(echo $nopi | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
#echo temp_root_file - $temp_root_file

 #opt=${tag}GFSa9t1nuC;               
nohup ./hist_nuSTORM_musig_GiBUU_v3_arg $E3_outAna9 $E4_outAna9 $E5_outAna9 $E6_outAna9 $E7_outAna9 > outAna9_see_plot_numu.log & #${root_prefix}

#temp_root_file=$(echo $PIZERO | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
#echo temp_root_file - $temp_root_file

 #opt=${tag}GFSPIZEROa7t4nuC;            
nohup ./hist_nuSTORM_musig_GiBUU_v3_arg $E3_outAna7 $E4_outAna7 $E5_outAna7 $E6_outAna7 $E7_outAna7 > outAna7_see_plot_numu.log & #${root_prefix}

exit