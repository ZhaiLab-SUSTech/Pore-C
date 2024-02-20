echo "cis" > cis_trans.txt
less MAPQ_1_sorted_dropdu.pairsam|grep -v "#"|awk -F "\t" '{if($2==$4) print}'|wc -l >> cis_trans.txt
echo "total" >> cis_trans.txt
less MAPQ_1_sorted_dropdu.pairsam|grep -v "#"|wc -l >> cis_trans.txt
#same method that used in ~/task/pore_c/porec_1000_filter_mainchr_result/pairs/cis_trans.sh and have been test with porec stat results.
