root_path="root_target"
flux_path="nuSTORM_flux_txt_files"

root_ext=".root"

nuFlux1p52_numu_flux=""
nuFlux1p52_nue_flux=""

nuFlux2p0_numu_flux=""
nuFlux2p0_nue_flux=""

nuFlux3p04_numu_flux=""
nuFlux3p04_nue_flux=""

nuFlux3p8_numu_flux=""
nuFlux3p8_nue_flux=""

nuFlux4p56_numu_flux=""
nuFlux4p56_nue_flux=""

nuFlux5p47_numu_flux=""
nuFlux5p47_nue_flux=""

#######################################################################################
mkexe.sh hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg -I${NUGENTKI}/include -lstyle -I${NUGENTKI}/style  -L${NUGENTKI}/style || exit

#nuFlux1p52, nuFlux2p0, nuFlux3p04, nuFlux3p8, nuFlux4p56, nuFlux5p47
#1p52, 2p0, 3p04, 3p8, 4p56, 5p47

for ii in $flux_path/*.txt
do

    if [[ "$ii" == *"Numu"* ]] || [[ "$ii" == *"numu"* ]] || [[ "$ii" == *"munu"* ]]
    then

        if [[ "$ii" == *"1p52"* ]]
        then
            nuFlux1p52_numu_flux=$ii
        fi
        if [[ "$ii" == *"2p0"* ]]
        then
            nuFlux2p0_numu_flux=$ii
        fi
        if [[ "$ii" == *"3p04"* ]]
        then
            nuFlux3p04_numu_flux=$ii
        fi
        if [[ "$ii" == *"3p8"* ]]
        then
            nuFlux3p8_numu_flux=$ii
        fi
        if [[ "$ii" == *"4p56"* ]]
        then
            nuFlux4p56_numu_flux=$ii
        fi
        if [[ "$ii" == *"5p47"* ]]
        then
            nuFlux5p47_numu_flux=$ii
        fi

    fi

    if [[ "$ii" == *"Nue"* ]] || [[ "$ii" == *"nue"* ]] || [[ "$ii" == *"enu"* ]]
    then

        if [[ "$ii" == *"1p52"* ]]
        then
            nuFlux1p52_nue_flux=$ii
        fi
        if [[ "$ii" == *"2p0"* ]]
        then
            nuFlux2p0_nue_flux=$ii
        fi
        if [[ "$ii" == *"3p04"* ]]
        then
            nuFlux3p04_nue_flux=$ii
        fi
        if [[ "$ii" == *"3p8"* ]]
        then
            nuFlux3p8_nue_flux=$ii
        fi
        if [[ "$ii" == *"4p56"* ]]
        then
            nuFlux4p56_nue_flux=$ii
        fi
        if [[ "$ii" == *"5p47"* ]]
        then
            nuFlux5p47_nue_flux=$ii
        fi

    fi

done
#nuFlux1p52, nuFlux2p0, nuFlux3p04, nuFlux3p8, nuFlux4p56, nuFlux5p47
#1p52, 2p0, 3p04, 3p8, 4p56, 5p47
mkdir -p root_target
for ii in $root_path/*/*
do

    echo root file dir - $ii
    
    nuFlux1p52_outana9=""
    nuFlux2p0_outana9=""
    nuFlux3p04_outana9=""
    nuFlux3p8_outana9=""
    nuFlux4p56_outana9=""
    nuFlux5p47_outana9=""

    nuFlux1p52_outana7=""
    nuFlux2p0_outana7=""
    nuFlux3p04_outana7=""
    nuFlux3p8_outana7=""
    nuFlux4p56_outana7=""
    nuFlux5p47_outana7=""

    for jj in $ii/*.root
    do

        if [[ "$jj" == *"outAna9"* ]] || [[ "$jj" == *"outana9"* ]]
        then

            if [[ "$jj" == *"1p52"* ]]
            then
                nuFlux1p52_outana9=$jj
            fi
            if [[ "$jj" == *"2p0"* ]]
            then
                nuFlux2p0_outana9=$jj
            fi
            if [[ "$jj" == *"3p04"* ]]
            then
                nuFlux3p04_outana9=$jj
            fi
            if [[ "$jj" == *"3p8"* ]]
            then
                nuFlux3p8_outana9=$jj
            fi
            if [[ "$jj" == *"4p56"* ]]
            then
                nuFlux4p56_outana9=$jj
            fi
            if [[ "$jj" == *"5p47"* ]]
            then
                nuFlux5p47_outana9=$jj
            fi

        fi


        if [[ "$jj" == *"outAna7"* ]] || [[ "$jj" == *"outana7"* ]]
        then

            if [[ "$jj" == *"1p52"* ]]
            then
                nuFlux1p52_outana7=$jj
            fi
            if [[ "$jj" == *"2p0"* ]]
            then
                nuFlux2p0_outana7=$jj
            fi
            if [[ "$jj" == *"3p04"* ]]
            then
                nuFlux3p04_outana7=$jj
            fi
            if [[ "$jj" == *"3p8"* ]]
            then
                nuFlux3p8_outana7=$jj
            fi
            if [[ "$jj" == *"4p56"* ]]
            then
                nuFlux4p56_outana7=$jj
            fi
            if [[ "$jj" == *"5p47"* ]]
            then
                nuFlux5p47_outana7=$jj
            fi
        fi
        
    done

    pdf_output_path="$ii/"

    numu_flux="$nuFlux1p52_numu_flux $nuFlux2p0_numu_flux $nuFlux3p04_numu_flux $nuFlux3p8_numu_flux $nuFlux4p56_numu_flux $nuFlux5p47_numu_flux"
    nue_flux="$nuFlux1p52_nue_flux $nuFlux2p0_nue_flux $nuFlux3p04_nue_flux $nuFlux3p8_nue_flux $nuFlux4p56_nue_flux $nuFlux5p47_nue_flux"

    if [[ "$ii" == *"Numu"* ]] || [[ "$ii" == *"numu"* ]] || [[ "$ii" == *"munu"* ]]
    then

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $nuFlux1p52_outana7 $nuFlux2p0_outana7 $nuFlux3p04_outana7 $nuFlux3p8_outana7 $nuFlux4p56_outana7 $nuFlux5p47_outana7 $pdf_output_path $numu_flux > $ii/outAna7_see_plot_numu.log & #${root_prefix}
        echo numu outana7 input $nuFlux1p52_outana7 $nuFlux2p0_outana7 $nuFlux3p04_outana7 $nuFlux3p8_outana7 $nuFlux4p56_outana7 $nuFlux5p47_outana7 $ii $numu_flux

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $nuFlux1p52_outana9 $nuFlux2p0_outana9 $nuFlux3p04_outana9 $nuFlux3p8_outana9 $nuFlux4p56_outana9 $nuFlux5p47_outana9 $pdf_output_path $numu_flux > $ii/outAna9_see_plot_numu.log & #${root_prefix}
        echo numu outana9 input $nuFlux1p52_outana9 $nuFlux2p0_outana9 $nuFlux3p04_outana9 $nuFlux3p8_outana9 $nuFlux4p56_outana9 $nuFlux5p47_outana9 $ii $numu_flux 

    fi

    if [[ "$ii" == *"Nue"* ]] || [[ "$ii" == *"nue"* ]] || [[ "$ii" == *"enu"* ]]
    then

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $nuFlux1p52_outana7 $nuFlux2p0_outana7 $nuFlux3p04_outana7 $nuFlux3p8_outana7 $nuFlux4p56_outana7 $nuFlux5p47_outana7 $pdf_output_path $nue_flux > $ii/outAna7_see_plot_nue.log & #${root_prefix}
        echo nue outana7 input $nuFlux1p52_outana7 $nuFlux2p0_outana7 $nuFlux3p04_outana7 $nuFlux3p8_outana7 $nuFlux4p56_outana7 $nuFlux5p47_outana77 $ii $nue_flux

        nohup ./do_hist_nuSTORM_musig_GiBUU_v4_evt_rate_arg.sh $nuFlux1p52_outana9 $nuFlux2p0_outana9 $nuFlux3p04_outana9 $nuFlux3p8_outana9 $nuFlux4p56_outana9 $nuFlux5p47_outana9 $pdf_output_path $nue_flux > $ii/outAna9_see_plot_nue.log & #${root_prefix}
        echo nue outana9 input $nuFlux1p52_outana9 $nuFlux2p0_outana9 $nuFlux3p04_outana9 $nuFlux3p8_outana9 $nuFlux4p56_outana9 $nuFlux5p47_outana9 $ii $nue_flux

    fi

done

echo done