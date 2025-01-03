export root_path="root_target"
export flux_path="nuSTORM_flux_txt_files"

export root_ext=".root"

E3_numu_flux=""
E3_nue_flux=""

E4_numu_flux=""
E4_nue_flux=""

E5_numu_flux=""
E5_nue_flux=""

E6_numu_flux=""
E6_nue_flux=""

E7_numu_flux=""
E7_nue_flux=""

#######################################################################################
mkexe.sh hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

for ii in $flux_path/*.txt
do

    if [[ "$ii" == *"Numu"* ]] || [[ "$ii" == *"numu"* ]] || [[ "$ii" == *"munu"* ]]
    then

        if [[ "$ii" == *"E3"* ]]
        then
            E3_numu_flux=$ii
        fi
        if [[ "$ii" == *"E4"* ]]
        then
            E4_numu_flux=$ii
        fi
        if [[ "$ii" == *"E5"* ]]
        then
            E5_numu_flux=$ii
        fi
        if [[ "$ii" == *"E6"* ]]
        then
            E6_numu_flux=$ii
        fi
        if [[ "$ii" == *"E7"* ]]
        then
            E7_numu_flux=$ii
        fi

    fi

    if [[ "$ii" == *"Nue"* ]] || [[ "$ii" == *"nue"* ]] || [[ "$ii" == *"enu"* ]]
    then

        if [[ "$ii" == *"E3"* ]]
        then
            E3_nue_flux=$ii
        fi
        if [[ "$ii" == *"E4"* ]]
        then
            E4_nue_flux=$ii
        fi
        if [[ "$ii" == *"E5"* ]]
        then
            E5_nue_flux=$ii
        fi
        if [[ "$ii" == *"E6"* ]]
        then
            E6_nue_flux=$ii
        fi
        if [[ "$ii" == *"E7"* ]]
        then
            E7_nue_flux=$ii
        fi

    fi

done

mkdir -p root_target
for ii in $root_path/*/*
do

    echo root file dir - $ii
    
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

    numu_flux="$E3_numu_flux $E4_numu_flux $E5_numu_flux $E6_numu_flux $E7_numu_flux"
    nue_flux="$E3_nue_flux $E4_nue_flux $E5_nue_flux $E6_nue_flux $E7_nue_flux"

    if [[ "$ii" == *"Numu"* ]] || [[ "$ii" == *"numu"* ]] || [[ "$ii" == *"munu"* ]]
    then

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $E3_outana7 $E4_outana7 $E5_outana7 $E6_outana7 $E7_outana7 $pdf_output_path $numu_flux > $ii/outAna7_see_plot_numu.log & #${root_prefix}
        echo numu outana7 input $E3_outana7 $E4_outana7 $E5_outana7 $E6_outana7 $E7_outana7 $ii $numu_flux

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $E3_outana9 $E4_outana9 $E5_outana9 $E6_outana9 $E7_outana9 $pdf_output_path $numu_flux > $ii/outAna9_see_plot_numu.log & #${root_prefix}
        echo numu outana9 input $E3_outana9 $E4_outana9 $E5_outana9 $E6_outana9 $E7_outana9 $ii $numu_flux 

    fi

    if [[ "$ii" == *"Nue"* ]] || [[ "$ii" == *"nue"* ]] || [[ "$ii" == *"enu"* ]]
    then

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $E3_outana7 $E4_outana7 $E5_outana7 $E6_outana7 $E7_outana7 $pdf_output_path $nue_flux > $ii/outAna7_see_plot_nue.log & #${root_prefix}
        echo nue outana7 input $E3_outana7 $E4_outana7 $E5_outana7 $E6_outana7 $E7_outana7 $ii $nue_flux

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $E3_outana9 $E4_outana9 $E5_outana9 $E6_outana9 $E7_outana9 $pdf_output_path $nue_flux > $ii/outAna9_see_plot_nue.log & #${root_prefix}
        echo nue outana9 input $E3_outana9 $E4_outana9 $E5_outana9 $E6_outana9 $E7_outana9 $ii $nue_flux

    fi

done

echo done