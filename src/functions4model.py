import os
import pandas as pd
import sys
import glob
import re
import numpy as np
import itertools
from random import randint
import shutil
import signal
import subprocess
from pathlib import Path
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

base_path = Path(__file__).parent
###########################################################################
##ADJUST TOP LOCATIONS ON YOUR SYSTEM
##LOCATIONS OF NETMHC,NETCHOP, NETSTAB 
netmhcpan = '/data/CCRSB/apps/netMHCpan/netMHCpan-4.0/netMHCpan'
##path to netCHOP-3.1
netchop = '/data/CCRSB/apps/netCHOP/netchop-3.1/bin/netChop'
##path to netStab
netstab = '/data/CCRSB/apps/netMHCstabpan/netMHCstabpan-1.0/netMHCstabpan'
############################################################################

mhcI = (str(base_path)+'/mhc/MHC_list.txt')

def create_Fastas(df,name, loc = None):
    '''Takes df as input that contains the columns 'Mut nmer' & 'Wt nmer'
    and produces fasta files ready to be used for netMHC, netCHOP, and netStab. 
    A name is also taken to be used to name fastas'''
    
    ##set location if none given
    if loc is None:
        loc = os.getcwd()+'/'
    
    ######################
    ##Create fastas for ##
    ##netMHC            ##
    ######################
    if not 'Mut nmer' in df.columns or not 'Wt nmer' in df.columns:
        sys.exit('The Columns Mut nmer and Wt nmer are required.  One or both are missing.')
    
    ##remove any epitopes with stop codon in them
    df = df.loc[((df['Mut nmer']!=df['Wt nmer'])&(~df['Mut nmer'].str.contains('\*')))]

    ##get rid of any fastas if they are there
    if os.path.exists(loc+str(name)+'_EpitopesI_mut.fasta'):
        os.remove(loc+str(name)+'_EpitopesI_mut.fasta')
    if os.path.exists(loc+str(name)+'_EpitopesI_wt.fasta'):
        os.remove(loc+str(name)+'_EpitopesI_wt.fasta')

    for i,r in df.iterrows():
        sets = [('mut', 'Mut nmer'), ('wt', 'Wt nmer')]
        for tups in sets:
            if (len(r[tups[1]]) >=8):
                epitope_fasta=open(loc+str(name)+'_EpitopesI_{}.fasta'.format(tups[0]),'a')
                epitope_fasta.write('>'+str(i)+'\n'+r[tups[1]].strip()+'\n')
                epitope_fasta.close()

    for x in glob.glob(loc+'*Epitopes*.fasta'):
        y = open(x, 'a')
        y.write('\n')
        y.close()

def create_net_commands(MHC_list, name, loc = None ):
    '''This will create the netMHC commands that need to be run given a list of MHC
        a location of the fasta files (either EpitopesI_*.fasta for class I,
        these fastas can be produced with create_fastas function.'''

    #if no location given execute in cwd
    if loc is None:
        loc = os.getcwd()+'/'

    ##Extract list of usable HLAs
    ##CLASSI    
    MHCI_list = []
    for x in open(mhcI):
        MHCI_list.append(x.strip())
    
    ##Make commands to for netMHCpan, netMHCpanII, netCHOP
    netMHC_commands = []
    failed_HLAs= set()

    ##Class I commands
    for a in MHC_list:
        if a not in MHCI_list:
            failed_HLAs.add(a)
            continue
        for b in glob.glob(loc+'*EpitopesI_*.fasta'):
            fasta = b.split('/')[-1]
            if re.search('mut',fasta):
                ##binding
                netMHC_commands.append('{} {} -a {} -l 8,9,10,11,12 -BA -s -xls  -xlsfile {}{}_netMHCpan_{}_mut.txt > {}{}{}_mut.trash '.format(netmhcpan,b,a,loc,name,a.replace(':',''),loc,name,a.replace(':','')))
            else:
                netMHC_commands.append('{} {} -a {} -l 8,9,10,11,12 -BA -s -xls  -xlsfile {}{}_netMHCpan_{}_wt.txt > {}{}{}_wt.trash '.format(netmhcpan,b,a,loc,name,a.replace(':',''),loc,name,a.replace(':','')))

    ##netchop 
    for b in glob.glob(loc+'*EpitopesI_mut.fasta'):
        fasta = b.split('/')[-1]
        netMHC_commands.append('{} -v 0 -t 2 {} > {}{}_netCHOPCterm_mut.out'.format(netchop,b,loc,name))
        netMHC_commands.append('{} -v 1 -t 2 {} > {}{}_netCHOP20S_mut.out'.format(netchop,b,loc,name))
    for b in glob.glob(loc+'*EpitopesI_wt.fasta'):
        fasta = b.split('/')[-1]
        netMHC_commands.append('{} -v 0 -t 2 {} > {}{}_netCHOPCterm_wt.out'.format(netchop,b,loc,name))
        netMHC_commands.append('{} -v 1 -t 2 {} > {}{}_netCHOP20S_wt.out'.format(netchop,b,loc,name))
    
    ##Write out predicted HLAs that weren't in 
    ##netMHC lists
    if len(failed_HLAs) >0:
        o=open(loc+'failed_HLAs.txt', 'w')
        o.write(str('\n'.join(x for x in failed_HLAs)))
        o.close()
        
    ##return commands
    return netMHC_commands

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
        elif exitcode == 0:
             print('{} completed with no error.'.format(cmd1.split('>')[0]))
        return exitcode
    
    else:
        exitcode = subprocess.call(cmd, shell = True)
        if exitcode != 0:
            print("'{}' failed with exit code {} ".format(cmd, exitcode))
        elif exitcode == 0:
             print('{} completed with no error.'.format(cmd1.split('>')[0]))
        return exitcode

def anchor_cysteine_check(pep, AA):
    ##record if aa in anchor residues of
    ##peptide output list of positions (0-based)
    ##that need to be changed
    pep = pep.strip()
    anchor_residues = [1,2,len(pep)-1]
    cys_pos = [i for i, aa in enumerate(pep) if aa == AA]
    cys_anchors = [pos for pos in anchor_residues if pos in cys_pos]
    if len(cys_anchors)>0:
        return 1
    else:
        return 0
    
def is_anchor(mut,wt):
    mut_pos =[]
    for e,x in enumerate(itertools.zip_longest(mut,wt)):
        if x[0] != x[1]:
            mut_pos.append(e)

    if len([x for x in [0,1,2,len(mut)-1] if x in mut_pos])>0:
        return 1
    else:
        return 0
    
def mk_cutsite(a,b):
    '''Make cutsite function to help with 
        errors seen from netMHC where some int
        are reported as strings or where the 
        entry is corrupted not sure why this issue
        occurs'''
    if isinstance(a, int) and isinstance(b, int):
        return a+b
    elif (not isinstance(a, int)) or (not isinstance(b, int)) :
        return '-'

def clean_up_prediction_output(df, loc = None):
    '''This will clean up the output of commands generated by
        create_net_commands'''

    #if no location given execute in cwd
    if loc is None:
        loc = os.getcwd()+'/'


    ########################
    ## NETCHOP DATAFRAMES ##
    ########################

    print('making Cterm')
    Cterm_mut = pd.DataFrame()
    Cterm_wt = pd.DataFrame()

    if len(glob.glob(loc+'*netCHOPCterm_*.out')) <2:
        sys.exit('missing netCHOP.out file')

    for x in glob.glob(loc+'*netCHOPCterm_*.out'):
        df_Cterm=pd.DataFrame(columns=['pos','AA','C','score','ID'])
        f=open(x)
        for l in f:
            formatted_line=re.sub(r'  *','\t',l)
            c=formatted_line.strip().split('\t')
            if len(c) ==5:
                try:
                    float(c[0])
                except:
                    continue
                df_Cterm.loc[len(df_Cterm)]=np.array(c)
        df_Cterm[['ID', 'score', 'pos']]=df_Cterm[['ID', 'score', 'pos']].astype(float).fillna(np.nan)
        df_Cterm['look_up']=list(zip(df_Cterm['ID'],df_Cterm['pos']))
        df_Cterm=df_Cterm[['look_up','score']]
        if re.search('_mut.out', x):
            Cterm_mut = Cterm_mut.append(df_Cterm, ignore_index = True)
        else:
            Cterm_wt = Cterm_wt.append(df_Cterm, ignore_index = True)
        # os.remove(x)

    df_Cterm = Cterm_mut.merge(Cterm_wt, how = 'left', left_on ='look_up', right_on = 'look_up')

    print('making 20S')
    S20_mut = pd.DataFrame()
    S20_wt = pd.DataFrame()

    if len(glob.glob(loc+'*netCHOP20S_*.out')) <2:
        sys.exit('missing netCHOP20S.out file')
    
    for x in glob.glob(loc+'*netCHOP20S_*.out'):
        df_20S=pd.DataFrame(columns=['pos','AA','C','score','ID'])
        f=open(x)
        for l in f:
            formatted_line=re.sub(r'  *','\t',l)
            c=formatted_line.strip().split('\t')
            if len(c) ==5:
                try:
                    float(c[0])
                except:
                    continue
                df_20S.loc[len(df_20S)]=np.array(c)
        df_20S[['ID', 'score', 'pos']]=df_20S[['ID', 'score', 'pos']].astype(float).fillna(np.nan)
        df_20S['look_up']=list(zip(df_20S['ID'],df_20S['pos']))
        df_20S=df_20S[['look_up','score']]
        if re.search('_mut.out', x):
            S20_mut = S20_mut.append(df_20S, ignore_index = True)
        else:
            S20_wt = S20_wt.append(df_20S, ignore_index = True)
        # os.remove(x)


    df_20S = S20_mut.merge(S20_wt, how = 'left', left_on ='look_up', right_on = 'look_up')

    #adjust columns and add change
    df_Cterm.columns = ['look_up', 'Mut C-term netChop score', 'Wt C-term netChop score']
    df_Cterm['Cterm score change (mut - wt)'] = df_Cterm.apply(lambda x: x['Mut C-term netChop score'] - x['Wt C-term netChop score'], axis = 1)

    df_20S.columns = ['look_up', 'Mut 20S netChop score', 'Wt 20S netChop score']
    df_20S['20S score change (mut - wt)'] = df_20S.apply(lambda x: x['Mut 20S netChop score'] - x['Wt 20S netChop score'], axis = 1)

    netCHOP = df_Cterm.merge(df_20S, left_on ='look_up', right_on='look_up')

    ###############################
    ##create netMHCI final output##
    ###############################

    df_netMHCI_mut = pd.DataFrame()
    df_netMHCI_wt = pd.DataFrame()
    print('Formatting netMHcoutput')
    for x in glob.glob(loc+'*netMHCpan*.txt'):
        f=open(x)
        for l in f:
            hla=l.strip()
            break
        a = pd.read_csv(x, skiprows=1, sep='\t')
        a['HLA']=hla
        a['ep_length']=a['Peptide'].str.len()
        if re.search('mut', x.split('/')[-1]):
            df_netMHCI_mut = df_netMHCI_mut.append(a, ignore_index=True)
        else:
            df_netMHCI_wt = df_netMHCI_wt.append(a, ignore_index=True)
        # os.remove(x)

    ##merge mutant and wt together
    print('Merging mutant and Wild Type')
    netMHCI = df_netMHCI_mut.merge(df_netMHCI_wt, how = 'left',  left_on = ['Pos', 'ID', 'HLA', 'ep_length'], right_on = ['Pos', 'ID', 'HLA','ep_length'], suffixes = (' Mutant', ' Wild Type') )
    df_netMHCI_wt = None
    df_netMHCI_mut = None

    ##only retain the samples where the mutant is present
    netMHCI = netMHCI.loc[netMHCI['Peptide Mutant']!= netMHCI['Peptide Wild Type']]

    ##add in netCHOP scores
    netMHCI = netMHCI.drop_duplicates()
    netMHCI['Pos']=netMHCI['Pos'].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    netMHCI['ep_length']=netMHCI['ep_length'].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    netMHCI['cut site']=netMHCI.apply(lambda x: mk_cutsite(x['Pos'], x['ep_length']), axis=1)
    netMHCI['netChop_lookup']=list(zip(netMHCI['ID'],netMHCI['cut site']))
    netMHCI['netChop_lookup']=list(zip(netMHCI['ID'],netMHCI['cut site']))

    print('Merging netmhc with netCHOP')
    print(len(netCHOP), len(netMHCI))
    netMHCI = netMHCI.merge(netCHOP, left_on ='netChop_lookup', right_on = 'look_up', how='left').drop(columns=['netChop_lookup', 'look_up'])
    print('Dropping columns')
    netMHCI.drop(columns = ['Ave Mutant', 'NB Mutant', 'core Wild Type',\
                            'icore Wild Type', 'Ave Wild Type', \
                            'NB Wild Type', 'cut site'], inplace = True)
    
    print('Adding Wt:Mut Rank')
    netMHCI['Wt:Mut Rank'] = netMHCI.apply(lambda x: x['Rank Wild Type']/x['Rank Mutant'], axis = 1)
    netMHCI = netMHCI.fillna(value='-')
    print('Adding cys in anchor')
    netMHCI['cys in anchor'] = netMHCI.apply(lambda x: anchor_cysteine_check(x['Peptide Mutant'], 'C'), axis = 1)
    print('Adding is anchor')
    netMHCI['is anchor'] = netMHCI.apply(lambda x: is_anchor(x['Peptide Mutant'], x['Peptide Wild Type']), axis = 1)
    


    ##########################################
    ##combine results of netMHCI and netCHOP##
    ##with the variants to retain variant   ##
    ##information save all scores           ##
    ##########################################
    netMHCI = netMHCI[['ID','HLA', 'Peptide Wild Type','Peptide Mutant', 'core Mutant', 'icore Mutant',\
                      'Rank Mutant','Rank Wild Type', 'Wt:Mut Rank','nM Mutant',  'nM Wild Type', \
                      'Mut C-term netChop score', 'Wt C-term netChop score', 'Cterm score change (mut - wt)',\
                      'Mut 20S netChop score', 'Wt 20S netChop score','20S score change (mut - wt)','cys in anchor', 'is anchor']]

    print('merging in dataframe nmer info')
    netMHCI = df.merge(netMHCI, how = 'left', left_index= True, right_on='ID')
    netMHCI.drop(columns = 'ID', inplace = True)
    netMHCI.reset_index(drop=True, inplace = True)
    netMHCI.fillna(value='-', inplace=True)
    print('saving file')
    netMHCI.to_csv(loc+'all_netMHCI_no_Rank_cutoff.txt', sep = '\t', index = True)

    rm_trash='rm '+loc+'*.trash '+loc+'*.fasta '+loc+'*netMHCpan*.txt '+loc+'*netCHOP*'
    subprocess.call(rm_trash, shell = True)
    
def anchor_cysteine_check(pep, AA):
    ##record if aa in anchor residues of
    ##peptide output list of positions (0-based)
    ##that need to be changed
    pep = pep.strip()
    anchor_residues = [1,2,len(pep)-1]
    cys_pos = [i for i, aa in enumerate(pep) if aa == AA]
    cys_anchors = [pos for pos in anchor_residues if pos in cys_pos]
    if len(cys_anchors)>0:
        return 1
    else:
        return 0
    
def is_anchor(mut,wt):
    import itertools

    mut_pos =[]
    for e,x in enumerate(itertools.zip_longest(mut,wt)):
        if x[0] != x[1]:
            mut_pos.append(e)

    if len([x for x in [0,1,2,len(mut)-1] if x in mut_pos])>0:
        return 1
    else:
        return 0
    
def mk_cutsite(a,b):
    '''Make cutsite function to help with 
        errors seen from netMHC where some int
        are reported as strings or where the 
        entry is corrupted not sure why this issue
        occurs'''
    if isinstance(a, int) and isinstance(b, int):
        return a+b
    elif (not isinstance(a, int)) or (not isinstance(b, int)) :
        return '-'

def create_netSTAB_commands(loc=None):
    ##will generate fastas from all_netMHCI_no_Rank_cutoff.txt
    ## and create commands for netStab that can be run 
    ##in parallel.
    
    ##read in all class I epitopes
    if loc is None:
        loc = os.getcwd()+'/'
    
    if not os.path.exists(loc+'all_netMHCI_no_Rank_cutoff.txt'):
        sys.exit('No epitope file present.')
    
    df = pd.read_csv(loc+'all_netMHCI_no_Rank_cutoff.txt', sep = '\t', index_col=0)
    
    ##I only need  limited number of columns to generate the fastas for
    ##netSTABpan
    fasta_file = df[['HLA', 'Peptide Mutant', 'Peptide Wild Type']]
    fasta_file['mut_len'] = fasta_file['Peptide Mutant'].str.len()
    fasta_file['wt_len'] = fasta_file['Peptide Wild Type'].str.len()
    HLA = fasta_file.loc[fasta_file['HLA']!='-']['HLA'].unique().tolist()
    fasta_file.drop_duplicates(subset= ['Peptide Mutant', 'Peptide Wild Type'], inplace = True)
    
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
                        out.write('{}\n'.format(r['Peptide Mutant']))
                    else:
                        out.write('{}\n'.format(r['Peptide Wild Type']))
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

def divide_return(x,y):
    try:
        out = x/y
        return out
    except:
        return '-'

def clean_up_netSTAB_output(loc = None):
    '''This will clean up the output of commands generated by
        create_netMHC_commands'''
    

    #if no location given execute in cwd
    if loc is None:
        loc = os.getcwd()+'/'

    if not os.path.exists(loc+'all_netMHCI_no_Rank_cutoff.txt'):
        sys.exit('No file used to generate fasta not present.')

    df = pd.read_csv(loc+'all_netMHCI_no_Rank_cutoff.txt', sep = '\t', index_col=0)
    if 'ID' in df.columns:
        df.drop(columns = 'ID', inplace = True)

    ###############################
    ##create netSTAB final output##
    ###############################

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
    #add in STAB ratio column
    netSTAB['NetSTAB prediction Mutant']=netSTAB['NetSTAB prediction Mutant'].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    netSTAB['NetSTAB prediction Wild Type']=netSTAB['NetSTAB prediction Wild Type'].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    netSTAB['NetSTAB Mut:WT'] = netSTAB.apply(lambda x: divide_return(x['NetSTAB prediction Mutant'],x['NetSTAB prediction Wild Type']), axis = 1)

    df = df.merge(netSTAB, how = 'left', left_on =['HLA','Peptide Mutant'] , right_on=['HLA','Peptide Mutant'])
    df.drop(columns = 'ID', inplace = True)
    if 'Peptide Wild Type_x' in df.columns:
        df.drop(columns = 'Peptide Wild Type_y', inplace = True)
        df.rename(columns = {'Peptide Wild Type_x': 'Peptide Wild Type'}, inplace = True)

    ##convert rank mutant to number so I can filter
    df['Rank Mutant'] = df['Rank Mutant'].apply(lambda x: pd.to_numeric(x, errors='coerce'))
    df_out = df.loc[(df['Rank Mutant']<=20)]
    ##fill back in nas
    df.fillna(value='-', inplace = True)
    df_out.fillna(value='-', inplace = True)
    df.to_csv(loc+'all_netMHCI_no_Rank_cutoff.txt', sep = '\t', index = True)
    df_out.to_csv(loc+'Final_NetMHCIpan4.0.txt', index = True, sep ='\t')   

    rm_trash='rm '+loc+'*.trash '+loc+'*.fasta'
    subprocess.call(rm_trash, shell = True)

def TAP_transport_score(epitope, peptide):
    '''Using the Method from Peters et al J Immunol 2003:171:1741-1749
    Identifying MHC Class I Epitopes by predicting the tap transport
    efficiency of epitope precursors and adjusted in netCTLpan publication
    a TAP transport score is calculated by averaging the score of the precursor 
    N_terminal peptides, downweighting to .2 and adding in the C-term score
    using the consensus scoring matrix N1,N2,N3 and Cterm'''
    
    Consensus_scoring_matrix = {'A' : {0:-1.56,1:-0.25,2:-0.10,3:0.24,4:-0.10,5:0.17,6:0.27,7:0.00,8:0.55},\
                            'C' : {0:0.05,1:-0.01,2:-0.02,3:0.11,4:0.09,5:0.05,6:0.00,7:-0.13,8:0.00},\
                            'D' : {0:1.37,1:1.42,2:1.83,3:-0.23,4:0.33,5:0.32,6:1.07,7:0.32,8:1.83},\
                            'E' : {0:1.65,1:0.02,2:1.51,3:0.08,4:0.54,5:-0.13,6:0.64,7:0.44,8:1.58},\
                            'F' : {0:1.03,1:-0.45,2:-1.05,3:-0.50,4:-0.26,5:0.08,6:-0.50,7:0.17,8:-2.52},\
                            'G' : {0:0.28,1:1.14,2:1.70,3:0.45,4:0.66,5:0.12,6:1.41,7:-0.38,8:1.41},\
                            'H' : {0:0.21,1:0.33,2:-0.23,3:-0.21,4:-0.11,5:-0.06,6:-0.19,7:0.39,8:0.55},\
                            'I' : {0:-0.11,1:-0.49,2:-0.62,3:-0.09,4:-0.42,5:-0.75,6:-0.94,7:0.45,8:-0.52},\
                            'K' : {0:-1.03,1:-0.41,2:0.09,3:-0.23,4:-0.08,5:-0.26,6:0.44,7:0.12,8:-0.45},\
                            'L' : {0:-0.50,1:0.09,2:-0.11,3:0.11,4:-0.34,5:0.02,6:-0.73,7:0.01,8:-0.94},\
                            'M' : {0:-0.38,1:-0.46,2:-0.58,3:-0.35,4:-0.26,5:0.30,6:-0.64,7:-0.11,8:-0.29},\
                            'N' : {0:-1.43,1:0.69,2:1.01,3:0.38,4:0.49,5:-0.27,6:0.16,7:0.33,8:1.33},\
                            'P' : {0:1.43,1:3.00,2:0.22,3:-0.04,4:-0.72,5:-0.13,6:-0.84,7:0.03,8:-0.09},\
                            'Q' : {0:0.47,1:-0.97,2:0.39,3:0.15,4:0.15,5:-0.07,6:0.34,7:0.26,8:0.12},\
                            'R' : {0:-1.34,1:-1.47,2:-0.42,3:-0.27,4:-0.32,5:-0.75,6:-0.09,7:-0.42,8:-1.47},\
                            'S' : {0:-0.56,1:-0.34,2:0.11,3:0.27,4:0.45,5:0.31,6:0.87,7:-0.51,8:2.26},\
                            'T' : {0:-0.12,1:-0.04,2:0.43,3:0.23,4:0.43,5:0.49,6:0.39,7:-0.46,8:0.72},\
                            'V' : {0:-0.49,1:-0.50,2:-0.71,3:0.27,4:0.37,5:-0.02,6:-0.29,7:0.10,8:-0.30},\
                            'W' : {0:0.54,1:-0.64,2:-1.65,3:-0.18,4:-0.78,5:0.31,6:-0.50,7:-0.63,8:-0.87},\
                            'Y' : {0:0.50,1:-0.67,2:-1.80,3:-0.18,4:-0.13,5:0.28,6:-0.87,7:0.02,8:-2.91}}
    
    try:
        loc = epitope.find(peptide)
        if loc > 0:
            N_score = 0
            precursor_N = epitope[loc-1:loc+2]
            peptide_N = peptide[0:3]
            for N_term in [precursor_N, peptide_N]:
                for e,x in enumerate(N_term):
                    N_score += Consensus_scoring_matrix[x][e]*-1
            
            ##Average and downsample N_score
            N_score =  (N_score/2)*.2

            ##get C-score
            C_score = Consensus_scoring_matrix[peptide[-1]][8]*-1
            
            ##get Tap score
            tap_score = C_score+N_score
        
        elif loc == 0: ##peptide at begining of epitope
            N_score = 0
            for e,x in enumerate(peptide[0:3]):
                    N_score += Consensus_scoring_matrix[x][e]*-1
            
            ##downsample
            N_score = N_score*.2
            
            ##get C-score
            C_score = Consensus_scoring_matrix[peptide[-1]][8]*-1
            
            ##get Tap score
            tap_score = C_score+N_score
        
        else:##peptide not in eptiope
            tap_score = np.nan
            
        return tap_score
        
    except:
        return np.nan
    
def hydrophobicity_score(peptide_string):
    '''will give you the sum of the kyte doolittle score for the
    for the peptide given.'''
    
    aa_hydro = {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
                'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
                'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
                'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
    
    hydro_sum = 0
    for x in peptide_string:
        try:
            hydro_sum += aa_hydro[x]
        except:
            return np.nan
    return hydro_sum

def find_Tcell_contact_mutations(Wtpep, Mutpep):
        
    contact_pos = [3,4,5,6]#0 based = 4,5,6,7 1 based
    
    if Wtpep =='-':
        return np.nan
    
    ##some of the peptides have spaces for some reason
    Wtpep = Wtpep.strip()
    Mutpep = Mutpep.strip()
    
    ##mutant pep lengths vary I want to go up 
    ##to but not include Cterm position
    Cterm = len(Mutpep)
    if Cterm != 8:
        fill_C = list((range(7,Cterm-1,1)))
        contact_pos.extend(fill_C)
    mut_pos = []
    for e,z in enumerate(itertools.zip_longest(Wtpep, Mutpep)):
        if z[0]!= z[1]:
            mut_pos.append(e)

    ##check against contacts
    check = [x for x in mut_pos if x in contact_pos]

    if len(check)>=1:
        return 1
    else:
        return 0

def MMP_Model_score(minimal_df):
    '''takes minimal df and produces MMP score for each minimal epitope'''
    
    ##load in models
    with open(str(base_path)+'/models/MmpModel.pickle', 'rb') as f:
        imr_mmp, scstd_mmp, model_mmp = pickle.load(f)
    
    ##feature columns
    mmpfeats = ['Gene expression decile', 'Present in RNA-seq data', 'Rank Mutant', 'Wt:Mut Rank', 'Mut C-term netChop score', 'Cterm score change (mut - wt)', 'Mut 20S netChop score', '20S score change (mut - wt)', 'NetSTAB prediction Mutant', 'Peptide TAP score Mut', 'Contact Hydrophobicity Mutant']

    for x in mmpfeats:
        minimal_df[x] = minimal_df[x].apply(pd.to_numeric, errors='coerce')


    minimal_df.loc[minimal_df['Peptide Wild Type']=='-', ['Wt:Mut Rank', 'T-Cell contact']] = np.nan,np.nan
    
    ##some samples are coming up with '-' in peptide mutant they can't be ran through predictions model
    score_df = minimal_df.loc[minimal_df['Peptide Mutant']!='-'].copy()
    score_df.reset_index(drop = True, inplace = True)

    ##fill NAs with mean of column
    _X = imr_mmp.transform(score_df[mmpfeats])
    print('done with imputer')

    ##Standardize all features
    _X = scstd_mmp.transform(_X)
    print('done with standard scaler')

    ##Score 
    score_df['MMP model score'] = model_mmp.predict_proba(_X)[:,1]

    print('MMPs scored')
    out_df = minimal_df.merge(score_df[mmpfeats+['HLA','MMP model score','Mut nmer', 'Peptide Mutant']],\
                              how = 'left', on = mmpfeats+['HLA','Mut nmer', 'Peptide Mutant'])

    out_df.drop_duplicates(inplace = True)
    out_df.reset_index(drop = True, inplace = True)
    
    return out_df

def NMER_model_score(feature_df):
    '''Scores the NMERS with the NMER model'''
    ##Load models
    with open(str(base_path)+'/models/NmerModel.pickle', 'rb') as f:
        imr_nmer, scstd_nmer, model_nmer = pickle.load(f)
   
    ##feature columns
    feats = ['minimal 1', 'minimal 2', 'minimal 3','minimal 4', 'minimal 5', 'Exome VAF decile']
        
    ##fill NAs
    _Score = imr_nmer.transform(feature_df[feats])

    ##Standardize all features
    _Score = scstd_nmer.transform(_Score)
    _proba = model_nmer.predict_proba(_Score)
    
    feature_df['Nmer Model score'] = [b for a,b in _proba]
    
    ##return df with score
    return feature_df
