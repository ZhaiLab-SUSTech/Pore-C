import pandas as pd
intersect_df = pd.read_table("/public/home/lizw/task/pore_c/rDNA3/intersect_0.1.tsv",header=None,names=['read_name1','start1','end1','align1','read_name2','start2','end2','align2'])
#keep only the read2 that located at the right side of the read1
intersect_df_single_dir = intersect_df.query('start1 <= start2')
intersect_df_sort = intersect_df_single_dir.sort_values(by=['read_name1','start1','end1'])
intersect_df_sort.to_csv('/public/home/lizw/task/pore_c/rDNA3/porec_filter_sorted_intersect_0.1.csv')
#must turn off the groupby sort to keep the original order !!!!
intersect_df_alignID = pd.DataFrame(intersect_df_sort.groupby('align1',sort=False).apply(lambda x:set(x['align2'])))
intersect_df_alignID.to_csv('/public/home/lizw/task/pore_c/rDNA3/porec_filter_intersect_df_alignID2.csv')
intersect_df_alignID['read_name'] = [x[0] for x in intersect_df_alignID.index.str.split("_")]

group = intersect_df_alignID.groupby(by='read_name',sort=False)

def filter(group):
    group_df = group[1]
    base_align = ''
    base_set = set()
    align_keep=[]
    singleton=[]
    for i in group_df.itertuples():
        align1,align_set,read_name = i
       #keeps align_set that contain interactions
        if len(align_set) < 2:
            singleton.append(align1)
        #no overlap, keep both
        elif len(base_set&align_set) == 0:
            align_keep.append(align1)
            base_set = align_set
            base_align = align1
        else:
            #if there are overlap, keeps only the largest set
            if len(align_set) > len(base_set):
                if base_align in align_keep:
                    align_keep.remove(base_align)
                align_keep.append(align1)
                base_set = align_set
                base_align = align1

            #the base set would be keep until it is smaller than other groups
            else:
                continue
    result_list = (align_keep,singleton)
    return result_list

from concurrent import futures
def filter_thread(func,iter_object,max_workers):
    workers = min(max_workers,len(iter_object))
    with futures.ThreadPoolExecutor(workers) as excutors:
        res = excutors.map(func,iter_object)
    return list(res)
    

from itertools import chain
align_both_list = filter_thread(filter,group,30)
#align_both_all = list(chain(*align_both_list))

interactions = []
singleton = []
for x in align_both_list:
    interactions.extend(x[0])
    singleton.extend(x[1])
    
intersect_df_filter = intersect_df_sort.query('align1 in @interactions')
intersect_df_filter.to_csv('/public/home/lizw/task/pore_c/rDNA3/porec_filter_intersect_df_filter2.csv')

all_align = []
all_align.extend(singleton)
all_align.extend(interactions)
all_filter = intersect_df_sort.query('align1 in @all_align')
all_filter.to_csv('porec_filter_all_filter.csv')
