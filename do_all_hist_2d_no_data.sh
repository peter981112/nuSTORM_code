export root_path="root_target"

export root_ext=".root"

#######################################################################################
mkexe.sh hist_2d_prod_no_data_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

for ii in $root_path/*/*.root
do

    echo root file - $ii
    export root_prefix=${ii%/*}/
    temp_root_file=$(echo $ii | cut -c$((${#root_prefix}+1))- | rev | cut -c$((${#root_ext}+1))- | rev)
    echo temp_root_file - $temp_root_file

    nohup ./do_hist_2d_no_data.sh $ii $root_prefix  > ${root_prefix}see_plot_${temp_root_file}.log &

done

echo done