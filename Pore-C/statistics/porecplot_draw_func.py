

#@Author: Zhuowen Li
#@LastEdit: 2021/10/5 下午9:44:25
#@Version: v2.13
#@Description: modify the line drawing


#@Author: Zhuowen Li
#@LastEdit: 2021/10/5 下午2:48:17
#@Version: v2.12 
#@Description: 

#@Author: Zhuowen Li
#@LastEdit: 2021/10/2 下午5:32:32
#@Version: v2.11
#@Description: generate matrix independently with heatmap draw,modify the plot specially for website


#@LastEdit: 2021/8/31 下午11:04:38
#@Version: v2.6
#@Description: 
#1.mark file columns: chrom start end name color
#2.provide optional sorting for the heatmap: maxily keep all reads in the plot

#@Author: Zhuowen Li
#@LastEdit: 2021/8/21 下午8:22:05
#@Version: v2.4
#@Description: 
#1.step sample according to the combinations of anchors, for cluster heatmap drawing
#2.anchor color ='grey'

#@Author: Zhuowen Li
#@LastEdit: 2021/8/28 下午4:57:53
#@Version: v2.5
#@Description: 
#1.change the anchor_detail file: add a new column to the anchor_filter file, retain the no anchor sites as normal sites in anchor_filter file, keep the complite reads information.
#2.correction: add brackets to all query and check the "in" judgments (must use list instead of pandas slice)
#3.correction: other['gene'] and other['mark'] would raise Keyerror


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from sys import exit
import os
import itertools
from collections import Counter

#calculating function

def chr_interval(df,region,column,binsize):
    pd.options.mode.chained_assignment = None
    df_region_filter_merge = pd.DataFrame([])
    chrinterval_all = pd.DataFrame([])
    df_bins_dict = {}
    for i in region.itertuples(index=False,name='Region'):
        region_chr,region_start,region_end,region_size,region_index = i
        df['chrom'] = df['chrom'].apply(lambda x:str(x))
        query_text = f'(chrom == @region_chr) and ({column} > @region_start) and ({column} <=  @region_end)'
        df_region_filter_temp = df.query(query_text)
        df_region_filter_temp['region_start'] = region_start
        df_region_filter_temp['region_end'] = region_end
        df_region_filter_temp['region_index'] = region_index
        df_cut,df_bins= pd.cut(df_region_filter_temp[column],bins=range(region_start,region_end+binsize-1,binsize),precision=0,retbins=True)
        df_region_filter_temp['interval']  = df_cut
        df_region_filter_merge = df_region_filter_merge.append(df_region_filter_temp)
        #add into intervals in region without target values calculated
        chrinterval_tmp = pd.DataFrame(pd.cut(np.arange(region_start+binsize,region_end+binsize,binsize),bins=df_bins,precision=0),columns = ['interval'])
        chrinterval_tmp['chrom'] = region_chr
        chrinterval_tmp['region_index'] = region_index
        chrinterval_all = chrinterval_all.append(chrinterval_tmp)
        df_bins_dict[region_index] = df_bins
    df_count = df_region_filter_merge.groupby(by=['region_index','interval']).count()[column]
    df_count.rename(f'{column}_count',inplace=True)
    return df_count,df_region_filter_merge,chrinterval_all,df_bins_dict

def region_reads(region):
    region_df = pd.read_csv(region,sep="\t",header=None,usecols=[0,1,2],names=['chrom','start','end'],converters={'chrom':str,'start':int,'end':int},comment="#")
    region_df['size'] = region_df['end']-region_df['start']
    region_df['region_index'] = region_df.index
    return region_df

def anchor_detail(anchor,df_merge_part,prefix):
    anchor_reads_all_df = pd.DataFrame([])
    for i in anchor.itertuples(index=False, name='Anchor'):
        #v2.3_9.1 input anchor files with name
        anchor_chr,anchor_start,anchor_end,anchor_name = i
        #filter out the alignment that related to the anchor sites
        criteria = '(chrom == @anchor_chr) and (pos <= @anchor_end) and (pos >= @anchor_start)'
        anchor_alignment = df_merge_part.query(criteria)
        anchor_alignment.to_csv(f'{prefix}_anchor_alignment.csv',sep='\t')
        
        #get the related reads 
        #filter out all alignments related to anchor sites related reads
        anchor_reads_list = list(anchor_alignment['read_name'])
        anchor_reads_df = df_merge_part.query('read_name in @anchor_reads_list')
        
        if len(anchor_reads_df) != 0:
            #mark the anchor alignment with anchor name,
            anchor_reads_df.loc[anchor_reads_df.eval(criteria),'anchor_name'] = anchor_name
            anchor_reads_all_df = anchor_reads_all_df.append(anchor_reads_df,ignore_index=True)
    
    anchor_reads_all_df.loc[:,'anchor_name'] = anchor_reads_all_df.loc[:,'anchor_name'].fillna('normal')
    anchor_reads_detail_df = anchor_reads_all_df.query('anchor_name != "normal"')
    
    anchor_reads_all_df.to_csv(f'{prefix}_anchor_reads_all.csv',sep='\t')
    anchor_reads_detail_df.to_csv(f'{prefix}_anchor_reads_detial.csv',sep='\t')
    
    return anchor_reads_detail_df


def anchor_statics(anchor_reads_detial_df,prefix):
    anchor_detial_site_nunique = anchor_reads_detial_df.groupby('read_name')['anchor_name'].nunique()
    anchor_detial_site_nunique.value_counts().to_csv(f'{prefix}_anchor_valuecount.csv',sep='\t')
    anchor_and2_list = anchor_detial_site_nunique[anchor_detial_site_nunique>=2].index.to_list()
    anchor_and2_detail = anchor_reads_detial_df.query('read_name in @anchor_and2_list')
    #only calculte the combination, would not calculate if the interaction caused by repeats of same anchor sites
    anchor_and2_detail_dropdu = anchor_and2_detail.drop_duplicates(['read_name','anchor_name'])
    
    for hubsize in [2,3,4]:
        anchor_detial_combinations = anchor_and2_detail_dropdu.groupby('read_name').apply(lambda x : list(itertools.combinations(x['anchor_name'],hubsize)))
        anchor_comb_list = list(itertools.chain(*anchor_detial_combinations))
        anchor_comb_df = pd.DataFrame.from_dict(dict(Counter(anchor_comb_list)),orient='index')
        fig, ax = plt.subplots(figsize=(10,6))
        if len(anchor_comb_df) >= 1:
            anchor_comb_df.sort_index().plot(kind='bar',ax=ax)
            ax.legend_.remove()
            plt.savefig(f'{prefix}_{hubsize}.png',format = 'png',dpi=300,bbox_inches = 'tight')
            anchor_comb_df.to_csv(f'{prefix}_{hubsize}.csv')
        else:
            continue

def anchor_reads(anchor,anchor_reads_detial_df,df_merge_part,region,column,binsize,anchor_mode,prefix,write,**others):
    anchor_detial_site_nunique = anchor_reads_detial_df.groupby('read_name')['anchor_name'].nunique()
    
    if anchor_mode[:3] == 'and':
        anchorcount = int(anchor_mode[3:])
    elif anchor_mode == 'all':
        anchorcount = len(anchor)
    elif anchor_mode == 'or':
        anchorcount = 1
        
    #step sample according to the combinations of anchors, for cluster heatmap drawing
    
    anchor_and2_list = anchor_detial_site_nunique[anchor_detial_site_nunique>=anchorcount].index.to_list()
    anchor_and2_detail = anchor_reads_detial_df.query('read_name in @anchor_and2_list')
    anchor_and2_detail_dropdu = anchor_and2_detail.drop_duplicates(['read_name','anchor_name'])
    df_comb = pd.DataFrame([])
    for i in range(anchorcount,len(anchor)+1):
        anchor_detial_combinations = pd.DataFrame(anchor_and2_detail_dropdu.groupby('read_name').apply(lambda x : list(itertools.combinations(x['anchor_name'],i))))
        #print(anchor_detial_combinations)
        anchor_detial_combinations_dropmu = anchor_detial_combinations[anchor_detial_combinations[0].apply(lambda x:len(x)==1)]
        df_comb = df_comb.append(anchor_detial_combinations_dropmu)
    df_comb['new'] = df_comb[0].apply(lambda x: x[0])
    
    if others['step_sample'] is not None:
        sample_frac = others['step_sample']
        df_comb_sample = df_comb.groupby(by='new').sample(frac=sample_frac)
        df_comb_sample.to_csv(f'{prefix}_comb_sample_frac{sample_frac}.csv')
        anchor_and_list = df_comb_sample.index.to_list()
    else:
        #print('why_no_comb?')
        df_comb.to_csv(f'{prefix}_comb_sample.csv')
        anchor_and_list = anchor_detial_site_nunique[anchor_detial_site_nunique>=anchorcount].index.to_list()
        
#include the alignment informaiton with anchored reads, not the samle of anchor_reads_detial_df
    anchor_merge = df_merge_part.query('read_name in @anchor_and_list')
    anchor_interval_count,anchor_filter_merge,*anchr_rest = chr_interval(anchor_merge,region,'pos',binsize)
    
    if (len(anchor_filter_merge) != 0) and (write == True):
        anchor_filter_merge.to_csv(f'{prefix}_anchor_filter_merge.csv',sep='\t')
        return anchor_interval_count,anchor_filter_merge,anchor_and_list
    
    else:
        with open (f'{prefix}.report','a+') as report:
            report.write('No reads match anchor')
            exit()
        
def anchor_calcu(merge_keepall_part,anchor_interval_count,region_df,binsize,prefix,write):        
    binnor_interval_count,*binnor_rest,chrinterval_all,df_bins_dict= chr_interval(merge_keepall_part,region_df,'site',binsize)
    anchor_binnor_merge = pd.merge(anchor_interval_count,binnor_interval_count,on = ['interval','region_index'],how='left')
    anchor_binnor_merge.eval('bin_nor=pos_count/site_count',inplace=True)
    anchor_binnor_merge_all = pd.merge(anchor_binnor_merge,chrinterval_all,on = ['interval','region_index'],how='right') \
                              .fillna({'pos_count':0,'site_count':0,'bin_nor':0})
    anchor_binnor_merge_all['bin_logy'] = np.log10(anchor_binnor_merge_all['pos_count']+1)
    anchor_binnor_merge_all_sort = anchor_binnor_merge_all.sort_values(by=['region_index','interval']).reset_index()
    if write == True:
        anchor_binnor_merge_all_sort.to_csv(f'{prefix}_anchor_binnor_merge_all.csv',index = None, sep='\t')
    
    return anchor_binnor_merge_all_sort

def heatmap_line(i,j,anchor_filter_merge):
    heat = anchor_filter_merge[anchor_filter_merge['region_index']==i].groupby(by=['read_name','interval']).count()['pos'].astype(np.int16).unstack().fillna(0)
    heat_T = heat.T
    heat_a = pd.merge(heat_T,j,on='interval',how ='outer').sort_values(by='interval')
    heat_b = heat_a.T.fillna(0)
    heat_b.drop(['interval'],axis=0,inplace=True)
    heat_c = heat_b.astype(float)
    heat_c.columns = np.arange(len(j))
    reads = anchor_filter_merge['read_name'].drop_duplicates().sort_values()
    heat_d = pd.merge(reads,heat_c,left_on='read_name',right_index=True,how ='left').fillna(0).sort_values(by='read_name')
    heat_d.index = np.arange(len(heat_d))
    return heat_d


def heat_df(anchor_binnor_merge_all_sort,anchor_filter_merge):
    anchor_binnor_merge_all_group = anchor_binnor_merge_all_sort.groupby('region_index')
    heat_concat = pd.DataFrame([])
    read_concat = pd.DataFrame([])
    for i,j in anchor_binnor_merge_all_group:
        heat_d = heatmap_line(i,j,anchor_filter_merge)
        #for cluster slice
        heat_concat = pd.concat([heat_concat,heat_d.iloc[:,1:]],axis=1)
        read_concat = pd.concat([read_concat,heat_d.iloc[:,0]],axis=1)
    return heat_concat,read_concat

#drawing functions        
def remove_draw(ax,**remove_dict):
    if remove_dict['tick'] == True:
        ax.tick_params(
        axis='x',          
        which='both',      
        bottom=False,     
        top=False,
        left=False,
        labelbottom=False,
        direction='inout') 
    if remove_dict['spine'] == True:
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    if remove_dict['locator'] == True:
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
    
def porec_draw(df,ax,calcute,plotkind,ylim):
    ax.set_xlabel('')
    remove_draw(ax,tick=True,spine=True,locator=False)
    df_len = len(df)
    useColDt = dict(raw='pos_count', logy="bin_logy", binnor="bin_nor")
    colorDt = dict(raw='darkgreen', logy="purple", binnor="#ED6D00")
    titleDt = dict(
        raw='Contact Count', 
        logy="log${_10}$(Contact Count + 1)", 
        binnor="Normalized Contact Count"
    )
    useCol = useColDt[calcute]
    plotDt = {}
    if plotkind == 'bar':
        plotDt['width'] = 0.9

    plotDt['kind'] = plotkind
    if ylim is not None:
        if calcute == 'logy':
            plotDt['ylim'] = (0,np.log10(float(ylim)))
        else:
            plotDt['ylim'] = (0,float(ylim))

    plotDt['xlim'] = (0,df_len+1)
    plotDt['color'] = colorDt[calcute]
        
    #parameters can be load in in dict way
    df[useCol].plot(**plotDt)
    ax.set_ylabel(titleDt[calcute])

def site_temp(s,region_df,binsize):
    site = pd.read_csv(s,header=None,sep='\t',names=['chrom','site'],converters={'chrom':str, 'site':int},comment="#")
    site_interval_count,_1,chrinterval_all,_2 = chr_interval(site,region_df,'site',binsize)
    site_interval_count_df = site_interval_count.reset_index()
    site_interval_count_merge = pd.merge(site_interval_count_df,chrinterval_all,how='right').sort_values(by=['chrom','interval']).fillna({'site_count':0})
    site_result = site_interval_count_merge['site_count']
    return site_result

    
def bg_temp(bg,region,binsize):
    bed_graph = pd.read_csv(bg,header=None,sep='\t',names=['chrom','start','end','values'],converters={'chrom':str,'start':int,'end':int,'values':float},comment="#")
    bed_graph.eval('result = (end-start)*values',inplace=True)

    _1,bg_region_filter_start,chrinterval_all,_2 = chr_interval(bed_graph,region,'start',binsize) 
    _3,bg_region_filter_end,*others = chr_interval(bed_graph,region,'end',binsize) 

    bg_region_filter_start.rename(columns={'interval':'start_interval'},inplace=True)
    bg_region_filter_end.rename(columns={'interval':'end_interval'},inplace=True)
    bg_region_filter_start_end = pd.merge(bg_region_filter_start,bg_region_filter_end,how='outer')

    #not_equeal
    bg_region_filter_start_end_notequal = bg_region_filter_start_end[bg_region_filter_start_end['start_interval']!=bg_region_filter_start_end['end_interval']]
    def middle_point(a,b):
        if type(a) == pd._libs.interval.Interval:
            return a.right
        elif type(b) == pd._libs.interval.Interval:
            return b.left
        else:
            pass


    bg_region_filter_start_end_notequal['middle_point'] = bg_region_filter_start_end_notequal.apply(lambda row: middle_point(row['start_interval'],row['end_interval']),axis=1).astype(int)
    bg_region_filter_start_end_notequal['start'] = bg_region_filter_start_end_notequal['start'].astype(int)
    bg_region_filter_start_end_notequal['end'] = bg_region_filter_start_end_notequal['end'].astype(int)
    bg_region_filter_start_end_notequal.eval('start_dist_part=result*(middle_point-start)/(end-start)',inplace=True,engine='python')
    bg_region_filter_start_end_notequal.eval('end_dist_part=result*(end-middle_point)/(end-start)',inplace=True,engine='python')

    bg_region_filter_notequal_a = bg_region_filter_start_end_notequal.loc[:,('chrom','start_interval','start_dist_part')].rename(columns={'start_interval':'interval','start_dist_part':'result'})
    bg_region_filter_notequal_b = bg_region_filter_start_end_notequal.loc[:,('chrom','end_interval','end_dist_part')].rename(columns={'end_interval':'interval','end_dist_part':'result'})

    #equal
    bg_region_filter_start_end_equal = bg_region_filter_start_end[bg_region_filter_start_end['start_interval']==bg_region_filter_start_end['end_interval']]
    bg_region_filter_start_end_equal_part = bg_region_filter_start_end_equal.loc[:,('chrom','start_interval','result')].rename(columns={'start_interval':'interval'})

    #all
    bg_region_filter_start_end_all = pd.concat([bg_region_filter_notequal_a,bg_region_filter_notequal_b,bg_region_filter_start_end_equal_part])
    bg_region_filter_start_end_all_drop_sort=bg_region_filter_start_end_all[~bg_region_filter_start_end_all['interval'].isin([0])].sort_values(by=['chrom','interval'])
    bg_region_filter_start_end_all_final = pd.merge(bg_region_filter_start_end_all_drop_sort,chrinterval_all,on=['chrom','interval'],how='right').fillna({'result':0})
    
    #sum
    bg_result = bg_region_filter_start_end_all_final.groupby(by=['chrom','interval']).sum()['result'].reset_index()['result']
    return bg_result

        
def track_draw(track_list,track_labels,track_ylims,mainplot,region_df,binsize,gs,color_order,kind):
    color_list = ['#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    track_order = 0
    for track in track_list: 
        ax_track = mainplot.add_subplot(gs[track_order])
        ax_track.set_xlabel('')
        ax_track.set_ylabel(track_labels[track_order])
        #ax_track.yaxis.set_label_coords(-0.05,0.6)
        remove_draw(ax_track,tick=True,spine=True,locator=False)
        if kind == 'bg':
            track_result = bg_temp(track,region_df,binsize)
        elif kind == 'site':
            track_result = site_temp(track,region_df,binsize)
        
        TrackDt = {}
        TrackDt['kind'] = 'bar'
        TrackDt['width'] = 0.9
        TrackDt['ax'] = ax_track
        TrackDt['color'] = color_list[track_order+color_order]
        
        if track_ylims:
            ylim_up = track_ylims[track_order]
            if ylim_up != None:                       
                TrackDt['ylim'] = (0,float(ylim_up))
        
        track_result.plot(**TrackDt)
        track_order += 1

def gene_model_draw(ax,gene_isoform,region_iter,yflen):
    region_chr,region_start,region_end,region_size,region_index,region_cum = region_iter
    gene_isoform['chrom'] = gene_isoform['chrom'].apply(lambda x:str(x))
    gene_isoform_sub = gene_isoform.query('(chrom == @region_chr) and (chromEnd>=@region_start) and (chromStart <= @region_end)')
    gene_isoform_sub['gene_id'] = gene_isoform['name'].map(lambda x: x.split('.')[0])
    gene_isoform_sub_part =  gene_isoform_sub.loc[:,('chrom','chromStart','chromEnd','strand','thickStart','thickEnd','blockCount','blockSizes','blockStarts','gene_id')]
    
    for gene in gene_isoform_sub_part.itertuples():
        _,gene_chr,gene_start,gene_end,gene_strand,gene_thickStart,gene_thickEnd,gene_blockCount,temp_blockSizes,temp_blockStarts,gene_id = gene
        gene_blockSizes = np.fromstring(temp_blockSizes, sep=',', dtype='int')
        gene_blockStarts = np.fromstring(temp_blockStarts, sep=',', dtype='int') + gene_start

        gene_size = gene_end - gene_start
        gene_color = 'k'
        arrowprops = dict(arrowstyle="-|>", connectionstyle="angle", color = gene_color)
        
        height = 0.075/yflen
        bottom = 0.60/yflen
        
        
        if gene_strand == '+':
            ax.annotate('', xy=(gene_start+gene_size/10, bottom + height*2), xytext=(gene_start, bottom), arrowprops=arrowprops)
        else:
            ax.annotate('', xy=(gene_end-gene_size/10, bottom + height*2), xytext=(gene_end, bottom), arrowprops=arrowprops)
        
        
        gene_block = mp.Rectangle((gene_start,bottom),gene_size,0.02/yflen,color=gene_color, linewidth=0)
        ax.add_patch(gene_block)
        
        #ax.plot([gene_start, gene_end], [0, 0], color = gene_color)
        for exonstart, size in zip(gene_blockStarts, gene_blockSizes):
            if (exonstart == gene_start) and (exonstart+size == gene_end):
                utr_size = gene_thickStart-gene_start
                utr = mp.Rectangle((exonstart, bottom + 0-height/2), utr_size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                utr_size = gene_end-gene_thickEnd
                utr = mp.Rectangle((gene_thickEnd, bottom + 0-height/2), utr_size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                exon = mp.Rectangle((gene_thickStart, bottom + 0-height), gene_thickEnd-gene_thickStart, height*2, color=gene_color, linewidth=0)
                ax.add_patch(exon)
            elif exonstart + size <= gene_thickStart:
                # only 5'/ 3'UTR
                utr = mp.Rectangle((exonstart, bottom + 0-height/2), size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)
            elif (exonstart < gene_thickStart) and (exonstart + size > gene_thickStart):
                # exon with 5' / 3' UTR 
                utr_size = gene_thickStart-exonstart
                utr = mp.Rectangle((exonstart, bottom + 0-height/2), utr_size, height, color=gene_color, linewidth=0)
                exon = mp.Rectangle((exonstart+utr_size, bottom + 0-height), size-utr_size, height*2, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                ax.add_patch(exon)
            elif (exonstart >= gene_thickStart) and (exonstart + size <= gene_thickEnd):
                # regular exon
                exon = mp.Rectangle((exonstart, bottom + 0-height), size, height*2, color=gene_color, linewidth=0)
                ax.add_patch(exon)
            elif (exonstart < gene_thickEnd) and (exonstart + size) > gene_thickEnd:
                # exon with 3' / 5' UTR
                utr_size = exonstart + size - gene_thickEnd
                utr = mp.Rectangle((gene_thickEnd, bottom + 0-height/2), utr_size, height, color=gene_color, linewidth=0)
                exon = mp.Rectangle((exonstart, bottom + 0-height), size-utr_size, height*2, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                ax.add_patch(exon)
            elif exonstart >= gene_thickEnd:
                # only 3'/ 5'UTR
                utr = mp.Rectangle((exonstart, bottom + 0-height/2), size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)

            ax.annotate(gene_id, xy=((gene_start+gene_end)/2, bottom+height*2.5), ha='center')
            
def rec_part(region_iter,df,kind,ax,binsize,yflen,**other):
    region_chr,region_start,region_end,region_size,region_index,region_cum = region_iter
    df['chrom'] = df['chrom'].apply(lambda x:str(x))
    df_sub = df.query('(chrom == @region_chr) and (end >= @region_start) and (start <= @region_end)')
    df_sub = df[df['chrom']==region_chr]
    for k in df_sub.itertuples(index=False):
        if kind == 'mark':
            df_chr,df_start,df_end,df_name,df_color = k
            bottom = 0.20/yflen
            ec_color = 'black'
        if kind == 'anchor':
            df_chr,df_start,df_end,df_name = k
            df_color = '#C00000'
            bottom = 0.40/yflen
            ec_color = '#C00000'
        
        left = df_start
        width = df_end - df_start + 1
        ax.add_patch(plt.Rectangle((left,bottom),width,0.075/yflen,color=df_color,ec = ec_color))
        remove_draw(ax,tick=False,spine = False,locator = False)
        ax.spines.top.set_visible(False)
        for axis in ['left','bottom','right']:
            ax.spines[axis].set_alpha(0.7)
            ax.spines[axis].set_color('#bdbebd')
        
        ax.set_xlabel(f'chr{region_chr}')
        ax.patch.set_alpha(0)
        ax.set_yticks([])
        ax.set_xlim(left=region_start,right=region_end)
        x_ticks = np.arange(region_start,region_end,10*binsize)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks,rotation=45,ha='right')
            
def rec_draw(region_df_cum,mainplot,anchor,binsize,gs,yflen,anchormode,**other):
    last_cum = 0
    for region_iter in region_df_cum.itertuples(index=False):
        region_cum = region_iter[-1]
        ax_rec = mainplot.add_subplot(gs[0,last_cum:region_cum])
        if other is not None:
            if 'gene' in other:
                gene_model_draw(ax_rec,other['gene'],region_iter,yflen)
            if 'mark' in other:
                rec_part(region_iter,other['mark'],'mark',ax_rec,binsize,yflen)
        rec_part(region_iter,anchor,'anchor',ax_rec,binsize,yflen)
        last_cum = region_cum
