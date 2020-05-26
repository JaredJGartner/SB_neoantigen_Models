import pandas as pd
import numpy as np
import os
from mhcflurry import Class1PresentationPredictor
import glob
import re
import subprocess
import itertools
from subprocess import Popen, PIPE, STDOUT
import shutil
import signal
from pathlib import Path
from random import randint
import multiprocessing as mp
from multiprocessing import Process
import pickle
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import BaggingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.calibration import CalibratedClassifierCV

os.environ['KMP_DUPLICATE_LIB_OK']='True'##this gets around an issue that might arise on Mac install  "OMP: Error #15: Initializing libiomp5.dylib, but found libiomp5.dylib already initialized."

base_path = Path(__file__).parent
####################################################################################################
##path to netStab 
netstab = '/Users/gartnerjj/Desktop/neoantigen_git/NetMHCStabpan1.0/netMHCstabpan-1.0/netMHCstabpan'
##path to immunogenicity
immuno_loc = '/Users/gartnerjj/Desktop/neoantigen_git/version2.0/src/immunogenicity/predict_immunogenicity.py3'
#####################################################################################################


def WT2MutpercentileRank(wtrank,mutrank):
    try:
        ratio = round(wtrank/mutrank, 4)
        return ratio
    except:
        return np.nan

def MHCflurry(df, alleles_list):
    predictor = Class1PresentationPredictor.load()
    #create dicts of index and nmer
    mut_dict = {}
    wt_dict = {}

    for i,r in df.loc[((df['Mut nmer']!='-')& (df['Mut nmer'].str.len() >=8))].iterrows():
        mut_dict[i] = r['Mut nmer']
    for i,r in df.loc[((df['Wt nmer']!='-')& (df['Wt nmer'].str.len() >=8))].iterrows():
        wt_dict[i] = r['Wt nmer']
        
    ##for whatever reason the result "all" isn't outputting everything 
    ##this works fine if I loop through each HLA and since MHCflurry is 
    ##quick this shouldn't be much of a problem
    MHCFlurryMut = pd.DataFrame()
    MHCFlurryWt = pd.DataFrame()
    
    for a in alleles_list:
        MHCFlurryMut = MHCFlurryMut.append(predictor.predict_sequences(sequences=mut_dict, result = 'all', alleles = [a], peptide_lengths=(8,9,10,11,12),  use_flanks=True, verbose = 0), ignore_index = True)
        MHCFlurryWt = MHCFlurryWt.append(predictor.predict_sequences(sequences=wt_dict, result = 'all', alleles = [a], peptide_lengths=(8,9,10,11,12),  use_flanks=True, verbose = 0), ignore_index = True)

    ##add lengths in 
    MHCFlurryMut['peptide length'] = MHCFlurryMut['peptide'].str.len()
    MHCFlurryWt['peptide length'] = MHCFlurryWt['peptide'].str.len()

    ##rename some columns
    MHCFlurryMut.rename(columns={'peptide': 'Mutant peptide', 'affinity': 'MHCflurry mutant affinity',\
                                'best_allele':'HLA', 'affinity_percentile':'MHCflurry affinity percentile rank mutant',\
                                'processing_score':'MHCflurry processing score mutant', 'presentation_score':'MHCflurry presentation score mutant'},\
                    inplace = True)
    MHCFlurryWt.rename(columns={'peptide': 'Wild type peptide', 'affinity': 'MHCflurry wild type affinity',\
                                'best_allele':'HLA', 'affinity_percentile':'MHCflurry affinity percentile rank wild type',\
                                'processing_score':'MHCflurry processing score wild type', 'presentation_score':'MHCflurry presentation score wild type'},\
                    inplace = True)


    ##merge mutant and wild type and remove wild type peptied in mutant peptide column
    merged = MHCFlurryMut.merge(MHCFlurryWt, how = 'left', left_on = ['sequence_name', 'pos', 'HLA', 'peptide length'],\
                    right_on = ['sequence_name', 'pos', 'HLA', 'peptide length'] )

    merged = merged.loc[~(merged['Mutant peptide']==merged['Wild type peptide'])][['sequence_name', 'HLA', 'peptide length','Mutant peptide','Wild type peptide',\
                                                                                    'MHCflurry mutant affinity','MHCflurry wild type affinity','MHCflurry affinity percentile rank mutant',\
                                                                                    'MHCflurry affinity percentile rank wild type','MHCflurry processing score mutant', 'MHCflurry processing score wild type',\
                                                                                    'MHCflurry presentation score mutant','MHCflurry presentation score wild type',]]

    ##merge with original dataframe
    out = df.merge(merged, how = 'left', left_index=True, right_on='sequence_name').drop(columns='sequence_name')
    out['MHCFlurry Wt:Mut rank'] = out.apply(lambda x: WT2MutpercentileRank(x['MHCflurry affinity percentile rank wild type'],x['MHCflurry affinity percentile rank mutant']), axis = 1)

    ##return out
    return out

def create_netSTAB_commands(df, loc=None):
    ##will generate fastas from all_netMHCI_no_Rank_cutoff.txt
    ## and create commands for netStab that can be run 
    ##in parallel.
    
    ##I only need  limited number of columns to generate the fastas for
    ##netSTABpan
    fasta_file = df[['HLA', 'Mutant peptide', 'Wild type peptide']]
    fasta_file['mut_len'] = fasta_file['Mutant peptide'].str.len()
    fasta_file['wt_len'] = fasta_file['Wild type peptide'].str.len()
    HLA = fasta_file.loc[((fasta_file['HLA']!='-')&(~fasta_file['HLA'].isnull()))]['HLA'].unique().tolist()
    fasta_file.drop_duplicates(subset= ['Mutant peptide', 'Wild type peptide'], inplace = True)
    
    ##generate mutant and wt fastas
    samples = ['mut', 'wt']
    for s in samples:
        for l in fasta_file[s+'_len'].unique():
            if l < 8:
                continue
            df_l = fasta_file.loc[fasta_file[s+'_len']==l]
            ##write fasta 
            with open(loc+'NetStab_{}_{}.fasta'.format(l,s),'w') as out:
                for i,r in df_l.iterrows():
                    out.write('>{}\n'.format(i))
                    if s == 'mut':
                        out.write('{}\n'.format(r['Mutant peptide']))
                    else:
                        out.write('{}\n'.format(r['Wild type peptide']))
                out.write('\n')
    
    ##generate commands for each fastas for each class I HLA
    commands = []
    for h in HLA:
        for f in glob.glob(loc+'NetStab*.fasta'):
            l = f.split('/')[-1].split('_')[1]
            out = f.replace('.fasta','_{}.txt'.format(h.replace(':','')))
            trash = f.replace('.fasta','_{}.trash'.format(h.replace(':','')))
            commands.append('{} -l {} -f {} -a {} -xls -xlsfile {} > {}'.format(netstab, l,f,h,out,trash))
    
    return commands

def init_worker():
    ##Courtesy of WR at biowulf
    '''This makes the worker processes ignore SIGINT. that means that any SIGING
        sent goes to the main process which can shut down the pool. this is required
        for scancel to terminate a multiprocessing script cleanly'''
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def sub_process(cmd):
    ##adapted from WR at biowulf 
    '''will accept a string or a tuple the tuple was needed
       so i could pass an environmental variable [TMPDIR]for
       each command of netMHCI/netCHOP'''
    my_env = os.environ.copy()
    if isinstance(cmd, tuple):
        cmd1 = cmd[0]
        tmp_dir = cmd[1]
        ##if somthing failed before this the tmp_dir may remain (unlikely due to randomint) I need to remove it
        if os.path.exists('tmp/{}'.format(tmp_dir)):
            shutil.rmtree('tmp/{}'.format(tmp_dir))
        path = Path('tmp/{}'.format(tmp_dir))
        path.mkdir(parents=True)
        ##set TMPDIR
        my_env['TMPDIR'] = 'tmp/{}'.format(tmp_dir)
        exitcode = subprocess.call(cmd1, shell = True, env = my_env)
        if exitcode != 0:
            print("{} failed with exit code {} ".format(cmd1, exitcode))
        # elif exitcode == 0:
        #     print('{} completed with no error.'.format(cmd1.split('>')[0]))
        return exitcode
    
    else:
        exitcode = subprocess.call(cmd, shell = True)
        if exitcode != 0:
            print("'{}' failed with exit code {} ".format(cmd, exitcode))
        # elif exitcode == 0:
        #     print('{} completed with no error.'.format(cmd1.split('>')[0]))
        return exitcode

def clean_up_netSTAB_output(df, loc = None):
    '''This will clean up the output of commands generated by
        create_netMHC_commands'''
    

    #if no location given execute in cwd
    if loc is None:
        loc = os.getcwd()+'/'

    df_netSTAB_mut = pd.DataFrame()
    df_netSTAB_wt = pd.DataFrame()

    for x in glob.glob(loc+'*NetStab*.txt'):
        f=open(x)
        for l in f:
            hla=l.strip()
            break
        a = pd.read_csv(x, skiprows=1, sep='\t')
        a['HLA']=hla
        a['ep_length']=a['Peptide'].str.len()
        a.rename(columns = {'Pred' : 'NetSTAB prediction'}, inplace = True)
        if re.search('mut', x.split('/')[-1]):
            df_netSTAB_mut = df_netSTAB_mut.append(a, ignore_index=True)
        else:
            df_netSTAB_wt = df_netSTAB_wt.append(a, ignore_index=True)
        os.remove(x)

    ##merge mutant and wt together
    netSTAB = df_netSTAB_mut[['ID', 'NetSTAB prediction', 'HLA', 'Peptide']].merge(df_netSTAB_wt[['ID','NetSTAB prediction','HLA','Peptide']],\
                                                                    how = 'left',  left_on = ['ID', 'HLA'],\
                                                                    right_on = ['ID', 'HLA'],\
                                                                    suffixes = (' Mutant', ' Wild Type'))
    netSTAB.drop_duplicates(subset = ['NetSTAB prediction Wild Type','HLA','Peptide Mutant',\
                                'Peptide Wild Type','NetSTAB prediction Mutant'], inplace = True)

    df = df.merge(netSTAB, how = 'left', left_on =['HLA','Mutant peptide'] , right_on=['HLA','Peptide Mutant'])
    df.drop(columns=['ID', 'Peptide Mutant', 'Peptide Wild Type'], inplace = True)

    return df

def run_immunogeniciy(df, loc):
    '''runs IEDB immunogenicity '''
    
    ##read in minimal peptides
    mers = df.copy()
    
    ##create Mut text to feed into immunogenicity prediction from IEDB
    run = ['Mutant', 'Wild type']
    for col in run:
        with open('{}immuno_peptides.txt'.format(loc), 'w') as out:
            for x in mers.loc[~mers['{} peptide'.format(col)].isnull()]['{} peptide'.format(col)].unique().tolist():
                if not re.search('-', x) and not re.search('X',x):
                    out.write('{}\n'.format(x))

        ##run Immunogenicity Mutant Peptides
        cmd_python = 'python {} {}immuno_peptides.txt '.format(immuno_loc, loc)
        p1 = subprocess.Popen(cmd_python, shell = True, stdout=PIPE)
        immuno_out = p1.stdout.read()

        #convert stdout to pd dataframe
        data = [x.split(',') for x in immuno_out.decode('ascii').split('\n')[4:]]
        header = immuno_out.decode('ascii').split('\n')[3].split(',')
        immuno_out = pd.DataFrame(data=data, columns = header)

        ##merge data back ino starting df
        mers = mers.merge(immuno_out[['peptide', 'score']], how = 'left', left_on='{} peptide'.format(col), right_on='peptide')
        mers.rename(columns = {'score':'IEDB {} immunogenicity score'.format(col)}, inplace = True)
        ##drop peptide columns
        mers = mers.drop('peptide', axis=1)
    
    ## add in change in immunogenicity score
    mers['IEDB Mutant immunogenicity score'] = mers['IEDB Mutant immunogenicity score'].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    mers['IEDB Wild type immunogenicity score'] = mers['IEDB Wild type immunogenicity score'].apply(lambda x: pd.to_numeric(x, errors='coerce'))

    os.remove('{}immuno_peptides.txt'.format(loc))

    return mers

def hydrophobicity_score(peptide_string):
    '''will give you the sum of the kyte doolittle score for the
    for the peptide given.'''
    
    aa_hydro = {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
                'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
                'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
                'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
    try:
        hydro_sum = 0
        for x in peptide_string[3:-1]:
            try:
                hydro_sum += aa_hydro[x]
            except:
                return np.nan
        return hydro_sum
    except:
        return np.nan

def MMP_Model_score(minimal_df):
    '''takes minimal df and produces MMP score for each minimal epitope'''
    
    ##load in models
    with open(str(base_path)+'/models/mmp_Model.pickle', 'rb') as f:
        imr_mmp, scstd_mmp, model_mmp = pickle.load(f)
    
    ##feature columns
    mmpfeats = ['Gene expression decile', 'Present in RNA-seq data', 'MHCflurry affinity percentile rank mutant', 'MHCFlurry Wt:Mut rank', 'NetSTAB prediction Mutant', 'IEDB Mutant immunogenicity score',  'Contact hydrophobicity mutant', 'Exome VAF decile']

    for x in mmpfeats:
        minimal_df[x] = minimal_df[x].apply(pd.to_numeric, errors='coerce')


    minimal_df.loc[minimal_df['Wild type peptide']=='-', 'MHCFlurry Wt:Mut rank'] = np.nan
    
    ##'-' in peptide mutant can't be ran through predictions model
    score_df = minimal_df.loc[minimal_df['Mutant peptide']!='-'].copy()

    ##fill NAs with mean of column
    _X = imr_mmp.transform(score_df[mmpfeats])

    ##Standardize all features
    _X = scstd_mmp.transform(_X)

    ##Score 
    score_df['MMP model score'] = [b for a,b in model_mmp.predict_proba(_X)]

    print('MMPs scored')
    out_df = minimal_df.merge(score_df['MMP model score'],\
                              left_index=True, right_index=True)

    out_df.drop_duplicates(inplace = True)
    out_df.reset_index(drop = True, inplace = True)
    
    return out_df

def NMER_model_score(nmers, mmpdf):
    '''Scores the NMERS with the NMER model'''
    ##Load models
    with open(str(base_path)+'/models/nmer_Model.pickle', 'rb') as f:
        nmer_imr, nmer_scstd, nmer_Model = pickle.load(f)
   
    # ##feature columns
    # feats = ['minimal 1', 'minimal 2', 'minimal 3','minimal 4', 'minimal 5', 'Exome VAF decile']
    ##generate input
    mmpScore = pd.DataFrame(columns = ['Unique identifier']+['mmp score {}'.format(x+1) for x in range(0,2)])

    for x in [(mmpdf,mmpScore)]:
        group = x[0].groupby('Unique identifier')
        for k,g in group:
            out = [k]
            ##add top netMHCrank
            top2 = g['MMP model score'].nlargest(2).tolist()
            top2.extend([np.nan]*(2-len(top2)))
            out.extend(top2)

            ##add to table
            x[1].loc[len(x[1])] = out
    
    ##nmers
    nmers = nmers.merge(mmpScore, how = 'outer', left_on='Unique identifier', right_on='Unique identifier', indicator=True)
  
    ##Train nmer model
    Nmer_features = ['mmp score 1', 'mmp score 2']
    imr = Imputer(missing_values='NaN', strategy = 'mean', axis=0)
    _X = nmer_imr.transform(nmers[Nmer_features])

    ##Standardize all features
    _X = nmer_scstd.transform(_X)

    ##score
    nmers['Nmer score'] = [b for a,b in nmer_Model.predict_proba(_X)]

    return nmers

def check4wt(mer, mer_list):
    check = mer_list
    if any(mer in string for string in check):
        return True
    else:
        return False
