export root_path="root_target"

export root_ext=".root"
E3_outana9=""
E4_outana9=""
E5_outana9=""
E6_outana9=""
E7_outana9=""

E3_outana7=""
E4_outana7=""
E5_outana7=""
E6_outana7=""
E7_outana7=""
#######################################################################################
mkexe.sh hist_nuSTORM_musig_GiBUU_v3_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

mkdir -p root_target
for ii in $root_path/*/*
do


    echo root file dir - $ii
    
    for jj in $ii/*.root
    do

        if [[ "$jj" == *"outAna9"* ]] || [[ "$jj" == *"outana9"* ]]
        then

            if [[ "$jj" == *"E3"* ]]
            then
                E3_outana9=$jj
            fi
            if [[ "$jj" == *"E4"* ]]
            then
                E4_outana9=$jj
            fi
            if [[ "$jj" == *"E5"* ]]
            then
                E5_outana9=$jj
            fi
            if [[ "$jj" == *"E6"* ]]
            then
                E6_outana9=$jj
            fi
            if [[ "$jj" == *"E7"* ]]
            then
                E7_outana9=$jj
            fi

        fi



        if [[ "$jj" == *"outAna7"* ]] || [[ "$jj" == *"outana7"* ]]
        then

            if [[ "$jj" == *"E3"* ]]
            then
                E3_outana7=$jj
            fi
            if [[ "$jj" == *"E4"* ]]
            then
                E4_outana7=$jj
            fi
            if [[ "$jj" == *"E5"* ]]
            then
                E5_outana7=$jj
            fi
            if [[ "$jj" == *"E6"* ]]
            then
                E6_outana7=$jj
            fi
            if [[ "$jj" == *"E7"* ]]
            then
                E7_outana7=$jj
            fi
        fi
        
    done

    pdf_output_path="$ii/"

    nohup ./hist_nuSTORM_musig_GiBUU_v3_arg $E3_outana7 $E4_outana7 $E5_outana7 $E6_outana7 $E7_outana7 $pdf_output_path > $ii/outAna7_see_plot_numu.log & #${root_prefix}
    echo outana7 input $E3_outana7 $E4_outana7 $E5_outana7 $E6_outana7 $E7_outana7 $ii

    nohup ./hist_nuSTORM_musig_GiBUU_v3_arg $E3_outana9 $E4_outana9 $E5_outana9 $E6_outana9 $E7_outana9 $pdf_output_path > $ii/outAna9_see_plot_numu.log & #${root_prefix}
    echo outana9 input $E3_outana9 $E4_outana9 $E5_outana9 $E6_outana9 $E7_outana9 $ii

done

echo done