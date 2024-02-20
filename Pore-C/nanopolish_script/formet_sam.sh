#change the read name of porec generating bam file
samtools view -h -F 4 /public/home/lizw/task/pore_c/porec_1000_filter_mainchr_result/mapping/sortedDpnII_run06.all.bam|awk -F "\t" '{gsub(":[0-9]+:[0-9]+","",$1);print}' OFS="\t"|samtools view -b > formet_sorted_run06.all.bam
samtools view -h -F 4 /public/home/lizw/task/pore_c/porec_1000_filter_mainchr_result/mapping/sortedDpnII_run07.all.bam|awk -F "\t" '{gsub(":[0-9]+:[0-9]+","",$1);print}' OFS="\t"|samtools view -b > formet_sorted_run07.all.bam
samtools index formet_sorted_run06.all.bam
samtools index formet_sorted_run07.all.bam
