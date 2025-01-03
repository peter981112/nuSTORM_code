export E3_outAna9="root_target/numu/E3SpectraMuSig560Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E4_outAna9="root_target/numu/E4SpectraMuSig566Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E5_outAna9="root_target/numu/E5SpectraMuSig558Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E6_outAna9="root_target/numu/E6spectraMuSig557Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E7_outAna9="root_target/numu/E7SpectraMuSig556Numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"

export E3_outAna7="root_target/numu/E3SpectraMuSig560Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E4_outAna7="root_target/numu/E4SpectraMuSig566Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E5_outAna7="root_target/numu/E5SpectraMuSig558Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E6_outAna7="root_target/numu/E6spectraMuSig557Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"
export E7_outAna7="root_target/numu/E7SpectraMuSig556Numu_outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root"

export E3_nue_flux="nuSTORM_flux_txt_files/E3SpectraMuSig560Nue.txt"
export E3_numu_flux="nuSTORM_flux_txt_files/E3SpectraMuSig560Numu.txt"
export E4_nue_flux="nuSTORM_flux_txt_files/E4SpectraMuSig566Nue.txt"
export E4_numu_flux="nuSTORM_flux_txt_files/E4SpectraMuSig566Numu.txt"
export E5_nue_flux="nuSTORM_flux_txt_files/E5SpectraMuSig558Nue.txt"
export E5_numu_flux="nuSTORM_flux_txt_files/E5SpectraMuSig558Numu.txt"
export E6_nue_flux="nuSTORM_flux_txt_files/E6spectraMuSig557Nue.txt"
export E6_numu_flux="nuSTORM_flux_txt_files/E6spectraMuSig557Numu.txt"
export E7_nue_flux="nuSTORM_flux_txt_files/E7SpectraMuSig556Nue.txt"
export E7_numu_flux="nuSTORM_flux_txt_files/E7SpectraMuSig556Numu.txt"


export root_path="root_target/numu/"


#export root_prefix=${nopi%/*}/
#echo root_prefix $root_prefix

#export root_ext=".root"

#######################################################################################
mkexe.sh hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

#temp_root_file=$(echo $nopi | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
#echo temp_root_file - $temp_root_file

 #opt=${tag}GFSa9t1nuC;               
nohup ./hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg $E3_outAna9 $E4_outAna9 $E5_outAna9 $E6_outAna9 $E7_outAna9 $root_path $E3_numu_flux $E4_numu_flux $E5_numu_flux $E6_numu_flux $E7_numu_flux > $root_path/outAna9_see_plot_numu.log & #${root_prefix}

#temp_root_file=$(echo $PIZERO | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
#echo temp_root_file - $temp_root_file

 #opt=${tag}GFSPIZEROa7t4nuC;            
nohup ./hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg $E3_outAna7 $E4_outAna7 $E5_outAna7 $E6_outAna7 $E7_outAna7 $root_path $E3_numu_flux $E4_numu_flux $E5_numu_flux $E6_numu_flux $E7_numu_flux > $root_path/outAna7_see_plot_numu.log & #${root_prefix}

exit