import pandas as pd
import click 
import numpy as np

@click.command()
@click.option('--geno_file',help = 'geno bed')
@click.option('--met_file',help = 'file with single cytosine met level')
@click.option('--prefix',help='prefix')
def level_size(geno_file,met_file,prefix):
    geno_bed = pd.read_table(geno_file,header=None,names=['chrom','chr_start','chr_end'])
    met_df = pd.read_table(met_file,header=None, names=['chrom','start','end','strand','coverage','met_cout','met_percentage','kind'])
    binsize = 500000
    interval_Ser = pd.Series([])
    for i in met_df.groupby(by='chrom'):
        chrom,df = i
        chr_start = 0
        chr_end = geno_bed.loc[chrom-1,'chr_end']
        df_cut = pd.cut(df['start'],bins=np.arange(chr_start,chr_end+binsize,binsize),right=False)
        interval_Ser = interval_Ser.append(df_cut)

    interval_Ser.name = 'interval'
    met_df_interval = pd.merge(met_df,interval_Ser,left_index=True,right_index=True)
    met_df_interval_group = met_df_interval.groupby(by=['chrom','interval']).apply(lambda df: df.assign(met_level = df['met_cout'].sum()/df['coverage'].sum(),size = len(df)))
    met_df_interval_group.to_csv(f'/public/home/lizw/task/pore_c/methylation/{prefix}_level_size_detial.csv')
    met_df_interval_group.query('size>=4').reindex(['chrom','kind','interval','met_level','size'],axis=1).drop_duplicates().to_csv(f'/public/home/lizw/task/pore_c/methylation/{prefix}_500000_level_size_in_interval.csv')
if __name__ == '__main__':
    level_size()
