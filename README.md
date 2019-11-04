# SB Neoantigen scores

This i repository contains a python script for scoring mmps and nmers for their likelihood of recognition against a list of Class I alleles based upon the methods described in  
"Development of a model for ranking candidate HLA class I neoantigens based upon datasets of known neoepitopes" Gartner et al.

### Prerequisites

This code is written to run with python version 3.6 and requires NetMHCpan-4.0,NetChop3.1 and netMHCStab 1.0.  Links are provided below please accept user license and install on your system

```
NetMHCpan-4.0 (http://www.cbs.dtu.dk/services/NetMHCpan/)
NetChop 3.1 (http://www.cbs.dtu.dk/services/NetChop/)
NetMHCStab 1.0 (http://www.cbs.dtu.dk/services/NetMHCstab/)
```

### Installing

After cloning repository on system alter NetMHCpan-4.0,NetChop 3.1 and NetMHCStab paths in src/functions4model.py to their locations on your system. 

Clone repository

```
git clone https://github.com/JaredJGartner/SB_neoantigen_Models.git
```

change directory

```
cd SB_neoantigen_Models/
```

Alter paths. You can open src/functions4models.py and use your favorite text editor or run the commands below.

```
sed -i  "s#/netMHCpan-4.0/netMHCpan#/path/to/your/netMHCpan-4.0/netMHCpan#" src/functions4model.py
sed -i  "s#/netchop-3.1/bin/netChop#/path/to/your/netchop-3.1/bin/netChop#" src/functions4model.py
sed -i  "s#/netMHCstabpan-1.0/netMHCstabpan#/path/to/your/netMHCstabpan-1.0/netMHCstabpan#" src/functions4model.py

```


## Testing the code

Example input file is provided in examples so that you can test that everything is working properly

### The script takes as input 

1) An excel sheet with the following header:
    Unique identifier,Wt nmer,Mut nmer,Exome VAF decile,Gene expression decile,Present in RNA-seq data\n    
    a) Unique identifier = a unique ID to link all data back to\n
    b) Wt nmer = Wild type amino acid sequence\n
    c) Mut nmer = Mutant aminoc acid sequence\n
    d) Exome VAF decile = The variant allele frequency for all mutations broken into decile for more information see "Development of a model for ranking candidate HLA class I neoantigens based upon datasets of known neoepitopes" Gartner et al.\n
    e) Gene expression decile = The expression of the gene binned into deciles as based upon al genes expression for that sample.\n
    f) Present in RNA-seq data = 1 or 0  for whether mutation was found in RNA sequencing data 1 = yes, 0 = no.\n

2) a list of up to 6 Class I HLAs formatted to match netMHC input for example HLA-A01:01


Test to see outputs are the same as those given in examples

```
mkdir test
cp examples/example_input.xlsx test/
cd test
../src/ScoreEpitopes.py example_input.xlsx  HLA-A01:01 HLA-A02:01 HLA-B11:01 HLA-B08:01 HLA-C02:01 HLA-C08:01

```

This will produce 2 out files in this test directory same as seen in examples directory

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details
