import pandas as pd
sample_read_table = pd.read_csv('/public/home/lizw/task/pore_c/methylation/nanopolish/rep2/sample_read_table.csv')
sample_met_table = pd.read_csv('/public/home/lizw/task/pore_c/methylation/nanopolish/rep2/sample_met_table.csv')

#A positive value in the log_lik_ratio column indicates support for methylation
#all are CG, need not to be identified
#some are single nucleus acid resolution, but some are not. The judge are for all the CG inculded in the sequence.
#the context call across two framgment would be abandoned

cpu_count = 30
batch_size = len(sample_read_table) // cpu_count
batch_dict = dict(tuple(sample_read_table.groupby(np.arange(len(sample_read_table))//batch_size)))
#process or thread?
def multi_thread(func,iter_object,max_workers):
    workers = min(max_workers,len(iter_object))
    with futures.ThreadPoolExecutor(workers) as excutors:
        res = excutors.map(func,iter_object)
    return list(res)

def met_to_fragment_batch(key):
    df = batch_dict[key]
    met_level_list = []
    for items in df.itertuples(index=False):
        chrom,start,end,strand,read_name = items
        chrom_int = int(chrom)
        sub_df = sample_met_table.query('(read_name == @read_name) & (chromosome == @chrom_int) & (start >= @start) & (end < @end)')
        #according to Bing Ren: Only reads containing 2 or more CpGs on each end were included
        CG_count = sub_df['num_motifs'].sum()
        if CG_count >= 2:
            met_sub_df = sub_df.query('log_lik_ratio > 0')
            met_sub_value = met_sub_df['num_motifs'].sum()/CG_count
        else:
            met_sub_value = np.nan
        
        met_level_list.append(met_sub_value)
    
    df['met_level'] = met_level_list
    return df

met_and_read_df = pd.concat(multi_thread(met_to_fragment_batch,list(batch_dict.keys()),cpu_count),axis=0)
met_and_read_df.to_csv('/public/home/lizw/task/pore_c/methylation/met_structure/sample_fragment_met_levels.csv')
