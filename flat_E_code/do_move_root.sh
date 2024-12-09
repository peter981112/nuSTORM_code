cd root_target
mkdir -p nue numu anti_nue anti_numu

for ii in *
do

    #echo root file - $ii

    if [[ "$ii" == *"Numu"* ]] || [[ "$ii" == *"numu"* ]] || [[ "$ii" == *"munu"* ]]
    then
        if [[ "$ii" == *"anti"* ]]
        then
            echo anti numu root file $ii

            mv $ii anti_numu
        else
            echo numu root file $ii

            mv $ii numu
        fi
    fi

    if [[ "$ii" == *"Nue"* ]] ||[[ "$ii" == *"nue"* ]] ||[[ "$ii" == *"enu"* ]]
    then
        if [[ "$ii" == *"anti"* ]]
        then
            echo anti nue root file $ii

            mv $ii anti_nue
        else
            echo nue root file $ii

            mv $ii nue          
        fi
    fi

done
