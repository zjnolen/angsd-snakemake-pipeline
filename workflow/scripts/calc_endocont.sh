#!/usr/bin/env bash

(merged="${snakemake_input[merged]}"
paired="${snakemake_input[paired]}"

tot_merged=$(grep -E "^[0-9]+ \+ [0-9]+ in total" $merged | awk '{print $1}')
map_merged=$(grep -E "^[0-9]+ \+ [0-9]+ mapped" $merged | awk '{print $1}')
if [ $tot_merged == 0 ]; then
    perc_merged="No merged reads..."
else
    perc_merged=$(echo $tot_merged $map_merged | awk '{printf "%.2f", $2/$1*100}')
fi

if [ "$merged" == "$paired" ]; then
    tot_paired="NA"
    map_paired="NA"
    perc_paired="NA"
    tot_tot="$tot_merged"
    map_tot="$map_merged"
else
    tot_paired=$(grep -E "^[0-9]+ \+ [0-9]+ in total" $paired | awk '{print $1}')
    map_paired=$(grep -E "^[0-9]+ \+ [0-9]+ mapped" $paired | awk '{print $1}')
    if [ $tot_paired == 0 ]; then
        perc_paired="No paired reads..."
    else
        perc_paired=$(echo $tot_paired $map_paired | awk '{printf "%.2f", $2/$1*100}')
    fi
    tot_tot=$(echo $tot_merged"+"$tot_paired | bc)
    map_tot=$(echo $map_merged"+"$map_paired | bc)
fi

if [ $tot_tot == 0 ]; then
    perc_tot="No reads..."
else
    perc_tot=$(echo $tot_tot $map_tot | awk '{printf "%.2f", $2/$1*100}')
fi


echo -e "${snakemake_wildcards[sample]}\t"$perc_merged"\t"$perc_paired"\t"$perc_tot \
    > "${snakemake_output[endo]}") 2> "${snakemake_log[0]}"