def dbi_binary(exp, nsf, thresh=0, outfile=None):
    '''Binary index:
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
    df = dbi_binary(exp, nsf)
    exp.sparse = False
    zzz=exp.copy()
    zzz.data = sp.stats.rankdata(zzz.data, axis=0)
    df_rank=dbi_freqs(zzz,ns,outfile='dysbiosis_index_ranks_1s.txt')
    if outfile is not None:
            df.to_csv(outfile, sep='\t')
    return df
  