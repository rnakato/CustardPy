for build in hg19 hg38 mm10 mm39 rn7 galGal5 galGal6 ce11 danRer11 dm6 xenLae2 sacCer3
do
    genome=/work/Database/UCSC/$build/genome.fa
    for enzyme in HindIII DpnII MboI Sau3AI Arima
    do
        python ../misc/generate_site_positions.py $enzyme $build $genome &
    done
done

#genome=/work/Database/Database_fromDocker/Ensembl-GRCg6a/genome.fa
#for enzyme in HindIII DpnII MboI; do
#    python ../misc/generate_site_positions.py $enzyme galGal6 $genome &
#done
