import calour as ca
import pandas as pd
import numpy as np
import scipy as sp

def import_nsf():
    '''Import the table of "good" and "bad" bacteria from Abass et al. 2021 metaanalysis paper.
    amnona/paper-metaanalysis.
    Returns
    -------
    pd.DataFrane

    '''
    nsd = pd.read_csv('https://raw.githubusercontent.com/amnona/paper-metaanalysis/main/ratios/nonspecific/nonspecific-down_feature.txt',sep='\t',index_col=0)
    nsd['dir'] = 'down'
    nsu = pd.read_csv('https://raw.githubusercontent.com/amnona/paper-metaanalysis/main/ratios/nonspecific/nonspecific-up_feature.txt',sep='\t',index_col=0)
    nsu['dir'] = 'up'
    ns = nsd.merge(nsu,how='outer')
    ns['dir'].value_counts()
    return ns

def dbi_binary(exp, nsf, thresh=0, outfile=None):
    '''Binary index:
    Parameters
    ----------
    exp: ca.Experiment
        Witht the samples to calculate the index
    nsf: pd.DataFrame
        with feature sequence as index, dir='up'/'down'
    
    Returns
    -------
    pd.DataFrane
        sample_id as index
        'score': dbi score
    '''
    res={}
    ca.set_log_level('ERROR')
    upf = nsf[nsf['dir']=='up']['_feature_id'].values
    downf = nsf[nsf['dir']=='down']['_feature_id'].values
    exp = exp.filter_ids(nsf._feature_id.values)
    exp.sparse = False
    exp.data = (exp.data > thresh)
    for cid, cexp in exp.iterate():
        tt = cexp.filter_ids(upf)
        nup = tt.data.sum(axis=1)[0]
        tt = cexp.filter_ids(downf)
        ndown = tt.data.sum(axis=1)[0]
        dbi =  np.log2((nup+0.1) / (ndown+0.1))
        res[cid] = dbi
    df=pd.DataFrame(res.items(), columns=['SampleID','Dysbiosis_index'])
    df=df.set_index('SampleID')
    if outfile is not None:
        df.to_csv(outfile, sep='\t')
    return df

def dbi_freqs(exp, nsf, thresh=0, outfile=None):
    '''Binary index:
    Parameters
    ----------
    exp: ca.Experiment
        Witht the samples to calculate the index
    nsf: pd.DataFrame
        with feature sequence as index, dir='up'/'down'
    
    Returns
    -------
    pd.DataFrane
        sample_id as index
        'score': dbi score
    '''
    res={}
    ca.set_log_level('ERROR')
    upf = nsf[nsf['dir']=='up']['_feature_id'].values
    downf = nsf[nsf['dir']=='down']['_feature_id'].values
    exp = exp.filter_ids(nsf._feature_id.values)
    exp.sparse = False
#     exp.data = (exp.data > thresh)
    for cid, cexp in exp.iterate():
        tt = cexp.filter_ids(upf)
        nup = tt.data.sum(axis=1)[0]
        tt = cexp.filter_ids(downf)
        ndown = tt.data.sum(axis=1)[0]
        dbi =  np.log2((nup+0.1) / (ndown+0.1))
        res[cid] = dbi
    df=pd.DataFrame(res.items(), columns=['SampleID','Dysbiosis_index'])
    df=df.set_index('SampleID')
    if outfile is not None:
        df.to_csv(outfile, sep='\t')
    return df

def dbi_ranks(exp, nsf, thresh = 0, outfile = None):
    '''Ranked index: 
    Parameters
    ----------
    exp: calour.Experiment
        Witht the samples to calculate the index
    nsf: pd.DataFrame
        with feature sequence as index, dir='up'/'down'
    
    Returns
    -------
    pd.DataFrane
        sample_id as index
        'score': dbi score
    '''
    exp.sparse = False
    zzz=exp.copy()
    zzz.data = sp.stats.rankdata(zzz.data, axis=0)
    df_rank=dbi_freqs(zzz, nsf)
    if outfile is not None:
        df_rank.to_csv(outfile, sep='\t')
    return df_rank
  
## example:
#nsf = import_nsf()
#exp = ca.read_amplicon("feature-table.biom", "metadata.tsv", min_reads=2000, normalize=10000)
#dbi = dbi_ranks(exp, nsf)