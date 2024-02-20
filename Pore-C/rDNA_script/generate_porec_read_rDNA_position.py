import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyarrow
import pyarrow.parquet as pq
from pathlib import Path
from pathlib import PurePath
import seaborn as sns
import pypairix

porec_dir = '/public/home/lizw/task/pore_c/porec_1000_filter_mainchr_result'
prefix = 'DpnII_run13'
mpq=1

align_dir = Path(PurePath(porec_dir,'align_table'))
align_dir_porec_generator = align_dir.glob(prefix + "*pore_c.parquet")
align_batch_porec_merge_df = pyarrow.concat_tables(list(map(pq.read_table,align_dir_porec_generator))).to_pandas()
align_batch_porec_merge_df_pass = align_batch_porec_merge_df.query('pass_filter==True').query('mapping_quality>=@mpq')
align_batch_porec_merge_df_pass_sort = align_batch_porec_merge_df_pass.sort_values(by=['read_name','read_start']).reindex(['chrom','start','end','strand','read_length','read_name','align_idx','read_start','read_end'],axis=1)
align_batch_porec_merge_df_pass_sort['align_name'] = align_batch_porec_merge_df_pass_sort['read_name'].str.cat(align_batch_porec_merge_df_pass_sort['align_idx'].astype('str'),sep="_")

from concurrent import futures
def remove_edge_overlap_thread(func,iter_object,max_workers):
    workers = min(max_workers,len(iter_object))
    with futures.ThreadPoolExecutor(workers) as excutors:
        res = excutors.map(func,iter_object)
    return list(res)

def remove_edge_overlap(group):
    last_read_end = 0
    group_df = group[1]
    for items in group_df.itertuples():
        index,chrom,start,end,strand,read_length,read_name,align_idx,read_start,read_end,align_name = items
        if read_start < last_read_end:
            cut_length = last_read_end - read_start
            group_df.loc[index,'read_start'] = last_read_end
            if strand == True:
                group_df.loc[index,'start'] = group_df.loc[index,'start'] + cut_length
            elif strand == False:
                group_df.loc[index,'end'] = group_df.loc[index,'end'] - cut_length
        last_read_end = read_end
    return group_df

group = align_batch_porec_merge_df_pass_sort.groupby(by='read_name')
align_batch_porec_merge_df_pass_sort_mod = pd.concat(remove_edge_overlap_thread(remove_edge_overlap,group,30),axis=0)

align_batch_porec_merge_df_pass_sort_mod.to_csv('/public/home/lizw/task/pore_c/rDNA3/align_batch_porec_merge_df_pass_sort_cutover.csv')
align_batch_porec_merge_df_pass_sort_mod.reindex(['read_name','start','end','align_name'],axis=1).to_csv('/public/home/lizw/task/pore_c/rDNA3/read_chrom_position_cutover.bed',index=None,header=None,sep="\t")
