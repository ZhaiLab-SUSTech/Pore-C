
library(optparse)
library(cocor)
library(vroom)
library(dplyr)
library(broom)
library(future.apply)


# input data format 
#r1      r2      n1      n2
#0.98    -0.04   100     100

# output data format
#r1      r2      n1      n2      statistic       distribution    p.value
#0.98    -0.04   100     100     16.279369841026206      z       0


parser <- OptionParser()
parser <- add_option(parser, c("-i", "--infile"), type="character", 
						help="infile")
parser <- add_option(parser, c("-o", "--outfile"), type="character", 
						help="Output File")
parser <- add_option(parser, c("-t", "--threads"), type="integer", default=10, 
						help="cpu number")
opt <- parse_args(parser)

infile <- opt$infile
outfile <- opt$outfile
threads <- opt$threads

#infile <- "/public/home/yuym/tmp.txt"
#outfile <- "/public/home/yuym/tmp_out.txt"
dat <- vroom(infile)

options(future.globals.maxSize=1e9)

plan(multicore, workers = threads)

stat <- future_lapply(1:nrow(dat), function(x){
    res <- cocor.indep.groups(dat$r1[x], dat$r2[x], dat$n1[x], dat$n2[x])@fisher1925
    res <- do.call(bind_cols,res)
    return(res)
})

stat <- do.call(bind_rows,stat)

dat <- bind_cols(dat,stat)

vroom_write(dat,path=outfile)