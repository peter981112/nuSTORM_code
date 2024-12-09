export nopi="root_target/numu/outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_E7SpectraMuSig556Numu.root"

export PIZERO="root_target/numu/outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_E7SpectraMuSig556Numu.root"

export root_prefix=${nopi%/*}/
echo root_prefix $root_prefix

export root_ext=".root"

#######################################################################################
mkexe.sh hist_2d_prod_no_data_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

temp_root_file=$(echo $nopi | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
echo temp_root_file - $temp_root_file

 #opt=${tag}GFSa9t1nuC;               
nohup ./hist_2d_prod_no_data_arg $nopi $root_prefix  > ${root_prefix}see_plot_${temp_root_file}.log &

temp_root_file=$(echo $PIZERO | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
echo temp_root_file - $temp_root_file

 #opt=${tag}GFSPIZEROa7t4nuC;            
nohup ./hist_2d_prod_no_data_arg $PIZERO $root_prefix  > ${root_prefix}see_plot_${temp_root_file}.log &

exit