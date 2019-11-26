# SB Neoantigen scores

This i repository contains a python script for scoring mmps and nmers for their likelihood of recognition against a list of Class I alleles based upon the methods described in  
"Development of a model for ranking candidate HLA class I neoantigens based upon datasets of known neoepitopes" Gartner et al.

### Prerequisites

This code is written to run with python version 3.6 and requires NetMHCpan-4.0,NetChop3.1 and netMHCStab 1.0.  Links are provided below please accept user license and install on your system

```
NetMHCpan-4.0 (http://www.cbs.dtu.dk/services/NetMHCpan/)
NetMHCpan-2.8 (http://www.cbs.dtu.dk/services/NetMHCpan/)
NetChop 3.1 (http://www.cbs.dtu.dk/services/NetChop/)
NetMHCStabpan 1.0 (http://www.cbs.dtu.dk/services/NetMHCstabpan/)
**gawk is required

```

### Installing

After cloning repository on system alter NetMHCpan-4.0,NetChop 3.1 and NetMHCStab paths in src/functions4model.py to their locations on your system. 

Clone repository

```
git clone https://github.com/JaredJGartner/SB_neoantigen_Models.git
```

change directory and unzip MmpModel

```
cd SB_neoantigen_Models/src/models/
gunzip MmpModel.pickle.gz
```

Alter paths in SB_neoantigen_Models/src/functions4models.py and use your favorite text editor.

```
Using your favorite text editor edit the paths to netmhcpan4.0, netmhcstabpan, netChop to their location on your system and save changes.

```


## Testing the code

Example input file is provided in examples so that you can test that everything is working properly

### The script takes as input 

1) An excel sheet with the following header:
    Unique identifier,Wt nmer,Mut nmer,Exome VAF decile,Gene expression decile,Present in RNA-seq data  
Unique identifier = a unique ID to link all data back to; Wt nmer = Wild type amino acid sequence; Mut nmer = Mutant amino acid sequence; Exome VAF decile = The variant allele frequency for all mutations broken into decile for more information see "Development of a model for ranking candidate HLA class I neoantigens based upon datasets of known neoepitopes" Gartner et al.; Gene expression decile = The expression of the gene binned into deciles as based upon al genes expression for that sample; Present in RNA-seq data = 1 or 0  for whether mutation was found in RNA sequencing data 1 = yes, 0 = no.

2) a list of up to 6 Class I HLAs formatted to match netMHC input for example HLA-A01:01


Test to see outputs are the same as those given in examples

```
##starting from SB_neoantigen_Models
mkdir test
cp examples/example_input.xlsx test/
cd test
python ../src/ScoreEpitopes.py example_input.xlsx  HLA-A01:01 HLA-A02:01 HLA-B11:01 HLA-B08:01 HLA-C02:01 HLA-C08:01

```

This will produce 2 out files in this test directory same as seen in examples directory

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details
