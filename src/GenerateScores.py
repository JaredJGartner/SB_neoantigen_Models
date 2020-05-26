import argparse
import sys
sys.path.append('src/')
from prediction_modules import *

##Parse command line arguments##
parser=argparse.ArgumentParser(description='''Takes as input:
                                                    1) an excel with Unique Names, Mutant and Wild Type Nmers, the exome VAF decile,
                                                       the Gene Expression and Whether the Mutation is seen in RNA-seq.  
                                                    2) a list of up to 6 Class I MHC alleles formatted to match netMHCpan input example: HLA-A01:01''', prog='GenerateScores.py', usage='GenerateScores.py xlfile allelelist')

parser.add_argument('xlfile',metavar='xlfile', help='name of formatted input excel file')
parser.add_argument('alleles',metavar='alleles', nargs='+',type=str, help='list of up to 6 alleles')
parser.add_argument('-C','--cpus', default=2, type = int, help='number of cpus to use')

##assign arguments to variables
args=parser.parse_args()
xlsx = args.xlfile
al = args.alleles
cpu_count = args.cpus

print(al)

##get prefix for output from excel
fileID = xlsx.split('/')[-1].replace('.xlsx', '')

##check that cpus are available
if cpu_count > os.cpu_count():
    sys.exit('{} cpus aren\'t available.  {} = available cpus on your system.'.format(cpu_count, os.cpu_count()))

##exit if more than 6 alleles given
if len(al) > 6:
    sys.exit('Only accept lists of up to 6 Class I alleles. {} listed: {}.'.format(len(al), ','.join(al)))

if __name__ == '__main__':
    
    ##read input
    df = pd.read_excel(xlsx)##changed to mmpdf
    ##socre with MHCflurry
    print('running MHCflurry')
    mmpdf = MHCflurry(df, al)
    print('MHCflurry complete')

    #mk fasta dir
    if not os.path.exists('fasta/'):
        os.mkdir('fasta')

    ##check minimal epitopes for stability to MHC
    print('running netMHCSTABpan')
    netSTAB_commands = create_netSTAB_commands(mmpdf, loc='fasta/')
    
    ##convert to tuples so I can include a TMPDIR variable for
    ##each subprocess.
    rand_ = randint(1,50000)
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
    
    print('netMHCstabpan complete')
    
    ##cleanup netSTAB results
    mmpdf = clean_up_netSTAB_output(mmpdf, loc = 'fasta/')

    ##IEDB immunogenicity score
    ##mk tmp dir
    if not os.path.exists('tmp/'):
        os.mkdir('tmp')

    print('Running IEDB immunogenicity score')
    mmpdf = run_immunogeniciy(mmpdf, loc='tmp/')
    print('completed IEDB immunogenicity score')

    ##add in contact residue hydrophobicity score
    mmpdf['Contact hydrophobicity mutant'] = mmpdf.apply(lambda x: hydrophobicity_score(x['Mutant peptide']), axis = 1)
    mmpdf['Contact hydrophobicity wild type'] = mmpdf.apply(lambda x: hydrophobicity_score(x['Wild type peptide']), axis = 1)

    ##Score mmps
    print('Scoring MMP models')
    mmpdf = MMP_Model_score(mmpdf)

    mmpdf.to_csv('tmp/mmpdf.txt', sep = '\t', index = False)

    #Score Nmers
    print('Scoring NMers')
    nmers = NMER_model_score(df, mmpdf)

    ##I want to clean up out put where Mut nmer is completely within Wt nmer
    ##a score will be output because imputer fills the median of columns but 
    ##this could be confusing to users
    nmers.loc[nmers.apply(lambda x: check4wt(x['Mut nmer'], list(set(nmers['Wt nmer'].tolist()))), axis = 1), ['mmp score 1', 'mmp score 2', 'Nmer score']] = '-'
    nmers.drop(columns='_merge', inplace = True)
    mmpdf.loc[mmpdf.apply(lambda x: check4wt(x['Mut nmer'], list(set(mmpdf['Wt nmer'].tolist()))), axis = 1), 'MMP model score'] = '-'


    ##produce out excel
    writer = pd.ExcelWriter('{}_scored.xlsx'.format(fileID), engine='xlsxwriter')
    df.to_excel(writer, sheet_name='InputData', index = False)
    nmers.to_excel(writer, sheet_name='NmersScored', index = False)
    mmpdf.to_excel(writer, sheet_name='MmpsScored', index = False)
    writer.save()

    ##remove tmp and fasta
    if os.path.exists('tmp/'):
        shutil.rmtree('tmp')

    if os.path.exists('fasta'):
        shutil.rmtree('fasta')

