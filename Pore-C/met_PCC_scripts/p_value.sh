export R_LIBS_USER=/public/home/yuym/software/R_libs:$R_LIBS_USER
/public/home/yuym/software/miniconda3/bin/R --slave --args --infile 1k_filter.txt --outfile 1k_out.txt < /public/home/yuym/taiyi/data/Scripts/lib/fisherz.R 2>> tmp.log
/public/home/yuym/software/miniconda3/bin/R --slave --args --infile 100k_filter.txt --outfile 100k_out.txt < /public/home/yuym/taiyi/data/Scripts/lib/fisherz.R 2>> tmp.log
/public/home/yuym/software/miniconda3/bin/R --slave --args --infile 5k_filter.txt --outfile 5k_out.txt < /public/home/yuym/taiyi/data/Scripts/lib/fisherz.R 2>> tmp.log
/public/home/yuym/software/miniconda3/bin/R --slave --args --infile 25k_filter.txt --outfile 25k_out.txt < /public/home/yuym/taiyi/data/Scripts/lib/fisherz.R 2>> tmp.log
/public/home/yuym/software/miniconda3/bin/R --slave --args --infile 1M_filter.txt --outfile 1M_out.txt < /public/home/yuym/taiyi/data/Scripts/lib/fisherz.R 2>> tmp.log
/public/home/yuym/software/miniconda3/bin/R --slave --args --infile trans_filter.txt --outfile trans_out.txt < /public/home/yuym/taiyi/data/Scripts/lib/fisherz.R 2>> tmp.log
