import argparse
import sys
sys.path.append('src/')
from functions4model import *

##Parse command line arguments##
parser=argparse.ArgumentParser(description='''Takes as input:
                                                    1) an excel with Unique Names, Mutant and Wild Type Nmers, the exome VAF decile,
                                                       the Gene Expression and Whether the Mutation is seen in RNA-seq.  
                                                    2) a list of up to 6 Class I MHC alleles formatted to match netMHCpan input example: HLA-A01:01''', prog='ScoreEpitopes.py', usage='ScoreEpitopes.py xlfile allelelist')

parser.add_argument('xlfile',metavar='xlfile', help='name of formatted input excel file')
parser.add_argument('alleles',metavar='alleles', nargs='+',type=str, help='list of up to 6 alleles')
parser.add_argument('-C','--cpus', default=2, type = int, help='number of cpus to use')

##assign arguments to variables
args=parser.parse_args()
xlsx = args.xlfile
al = args.alleles
cpu_count = args.cpus

##get prefix for output from excel
fileID = xlsx.split('/')[-1].replace('.xlsx', '')


##check that cpus are available
if cpu_count > os.cpu_count():
    sys.exit('{} cpus aren\'t available.  {} = available cpus on your system.'.format(cpu_count, os.cpu_count()))

##exit if more than 6 alleles given
if len(al) > 6:
    sys.exit('Only accept lists of up to 6 Class I alleles. {} listed: {}.'.format(len(al), ','.join(al)))


if __name__ == '__main__':
    ##get dataframe input
    df = pd.read_excel(xlsx)

    ##mk tmp dir
    if not os.path.exists('tmp/'):
        os.mkdir('tmp')
    ##mk tmp dir
    if not os.path.exists('fasta/'):
        os.mkdir('fasta')

    ##create fastas from nmers
    create_Fastas(df, fileID, 'fasta/')

    ##get netchop and netMHC commands that will be ran
    netCommands = create_net_commands(al, fileID, 'fasta/')

    rand_ = randint(1,50000)
    ##convert to tuples so I can include a TMPDIR variable for
    ##each subprocess. I'm going to spread this over the cpus
    ##and I want to keep things cleam with a temp dir
    c = list(zip(netCommands,["{}_{}_{}".format(fileID, rand_, a_) for a_ in range(len(netCommands))]))

    ##Launch Pools
    p = mp.Pool(cpu_count, initializer = init_worker)
    #WR using Pool.map instead of apply_async. map is synchronous
    #so it blocks until all the commands have been executed
    try:
        #WR run the processing pool
        exitcodes = p.map(sub_process, c)
    #WR catch SIGINT and shut down cluster
    except (KeyboardInterrupt, SystemExit):
        p.terminate()
        p.join()
        sys.exit(1)
    else:
        p.close()
        p.join()
    
    ##check to make sure all netmhc executed fine before cleaning up
    if 1 in exitcodes:
        sys.exit('netMHC/netCHOP command failed.') 

    ##cleanup into 
    clean_up_prediction_output(df, 'fasta/')

    ##check minimal epitopes for stability to MHC
    print('generating netSTAB commands')
    netSTAB_commands = create_netSTAB_commands(loc='fasta/')

    ##convert to tuples so I can include a TMPDIR variable for
    ##each subprocess.
    c = list(zip(netSTAB_commands,["{}_{}_{}".format(fileID, rand_, a_) for a_ in range(len(netSTAB_commands))]))
    ##Launch Pools
    p = mp.Pool(cpu_count, initializer = init_worker)
    #WR using Pool.map instead of apply_async. map is synchronous
    #   so it blocks until all the commands have been executed
    try:
        #WR run the processing pool
        exitcodes = p.map(sub_process, c)
    #WR catch SIGINT and shut down cluster
    except (KeyboardInterrupt, SystemExit):
        p.terminate()
        p.join()
        sys.exit(1)
    else:
        p.close()
        p.join()

    ##check to make sure all netmhc executed fine before cleaning up
    if 1 in exitcodes:
        sys.exit('netSTAB command failed.') 

    if os.path.exists('tmp'):
        shutil.rmtree(('tmp'))

    ##cleanup netSTAB results
    clean_up_netSTAB_output(loc = 'fasta/')

    ##add in TAP transport score, hydrophobicity for the T-cell contact residues, and whether the mutation 
    ##falls in the T-cell contact residue
    print('add in TAP transport score, hydrophobicity for the T-cell contact residues')
    all_MHCI = pd.read_csv('fasta/all_netMHCI_no_Rank_cutoff.txt', sep = '\t', index_col=0)
    all_MHCI['Peptide TAP score Mut'] = all_MHCI.apply(lambda x: TAP_transport_score(x['Mut nmer'], x['Peptide Mutant']), axis = 1 )
    all_MHCI['Peptide TAP score Wt'] = all_MHCI.apply(lambda x: TAP_transport_score(x['Wt nmer'], x['Peptide Wild Type']), axis = 1 )
    all_MHCI['Peptide TAP score change (Mut - Wt)'] = all_MHCI.apply(lambda x: x['Peptide TAP score Mut'] - x['Peptide TAP score Wt'], axis = 1)
    all_MHCI['Contact Hydrophobicity Mutant'] = all_MHCI.apply(lambda x: hydrophobicity_score(x['Peptide Mutant'][3:-1]), axis = 1)
    all_MHCI['Contact Hydrophobicity Wild Type'] = all_MHCI.apply(lambda x: hydrophobicity_score(x['Peptide Wild Type'][3:-1]), axis = 1)
    all_MHCI['Contact Hydrophobicity Change (Wt - Mut)'] =  all_MHCI.apply(lambda x: x['Contact Hydrophobicity Wild Type']-x['Contact Hydrophobicity Mutant'], axis = 1)
    all_MHCI.fillna(value='-', inplace = True)
    all_MHCI['T-Cell contact'] = all_MHCI.apply(lambda x: find_Tcell_contact_mutations(x['Peptide Wild Type'],x['Peptide Mutant']), axis = 1)       
    all_MHCI['Rank Mutant'] = all_MHCI['Rank Mutant'].apply(lambda x: pd.to_numeric(x, errors='coerce'))

    #score minimals
    print('score minimals')
    all_MHCI = MMP_Model_score(all_MHCI)

    ##fill back in nas
    all_MHCI.fillna(value='-', inplace = True)
    
    ##Score Nmers
    #Ensure these two columns are numeric
    all_MHCI['Rank Mutant'] = all_MHCI['Rank Mutant'].apply(pd.to_numeric, errors = 'coerce')
    all_MHCI['MMP model score'] = all_MHCI['MMP model score'].apply(pd.to_numeric, errors = 'coerce')
    all_MHCI.sort_values(by= 'MMP model score', ascending = False).to_csv('{}_mmps_scored.txt'.format(fileID), sep = '\t', index = False)

    ##create dataframe to be used for Nmer model input
    MMPScores = all_MHCI[['Unique identifier','Exome VAF decile','MMP model score']]
    NMERinput = pd.DataFrame(columns=['Unique identifier','Exome VAF decile','minimal 1','minimal 2','minimal 3','minimal 4','minimal 5'])

    group = MMPScores.groupby(['Unique identifier','Exome VAF decile'])
    for k,g in group:
        out = [k[0],k[1]]
        top5Minimals = g['MMP model score'].nlargest(5).tolist()
        top5Minimals.extend([np.nan]*(5-len(top5Minimals)))

        ##add top ranked minimal
        out.extend(top5Minimals)
        NMERinput.loc[len(NMERinput)] = out
        
    ##final score for 25mers
    NMER = NMER_model_score(NMERinput)
    NMER = df.merge(NMER[['Unique identifier', 'Exome VAF decile', 'Nmer Model score']], how = 'left', left_on = ['Unique identifier', 'Exome VAF decile'] , right_on = ['Unique identifier', 'Exome VAF decile'])

    NMER.sort_values(by='Nmer Model score', ascending = False).to_csv('{}_nmers_scored.txt'.format(fileID), sep = '\t', index = False)

    ##move appropriate files from tmp and delete directory
    rmvcmd = 'rm -rf fasta/'
    subprocess.call(rmvcmd, shell=True)

