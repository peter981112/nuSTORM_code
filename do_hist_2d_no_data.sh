date

rootfile=$1

root_path=$2

echo rootfile ${rootfile}
echo root_path ${root_path}

#######################################################################################

#void hist_2d_prod_no_data(TString fn, string root_file_path)            
./hist_2d_prod_no_data_arg ${rootfile} ${root_path}  #> see${opt}.log &


date
exit