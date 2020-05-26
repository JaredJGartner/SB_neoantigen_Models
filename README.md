# SB Neoantigen scores

This repository contains a python script for scoring mmps and nmers for their likelihood of recognition against a list of Class I alleles based upon the methods described in  
"Development of a model for ranking candidate HLA class I neoantigens based upon datasets of known neoepitopes" Gartner et al.

### Prerequisites

This code is written to run with python version 3.6 and requires IEDB Immunogenicity tool, MHCflurry1.6 and netMHCStab 1.0.  
Links are provided below please accept user license and install on your system.
a conda environment .yml is distributed to create python environment mhcflurry-env to run model 

```
IEDB immunogenicity prediction Tool (http://tools.iedb.org/immunogenicity/download/)
NetMHCStabpan 1.0 (http://www.cbs.dtu.dk/services/NetMHCstabpan/)


```

### Installing

After cloning repository on system please read instructions full step by step description of set up found in instructions.txt. 
IEDB mmunogenicity needs to be converted to python 3. Instructions for doing so are in this file. 
After everything is installed alter IEDB immunogenicity and  NetMHCStab paths in src/prediction_modules.py to their locations on your system.


Clone repository

```
git clone https://github.com/JaredJGartner/SB_neoantigen_Models.git
```

change directory and unzip models

```
cd SB_neoantigen_Models/src/models/
gunzip mmp_Model.pickle.gz
gunzip nmer_Model.pickle.gz
```

## Alter paths in SB_neoantigen_Models/src/functions4models.py and use your favorite text editor.

```
Using your favorite text editor edit the paths to IEDB imuunogenicity, netmhcstabpan to their location on your system and save changes.
A full step by step description of setting model is found in instructions.txt

```


## Testing the code

Example input file is provided in examples so that you can test that everything is working properly

### The script takes as input 

1) An excel sheet with the following header:
    Unique identifier,Wt nmer,Mut nmer,Exome VAF decile,Gene expression decile,Present in RNA-seq data  
Unique identifier = a unique ID to link all data back to; Wt nmer = Wild type amino acid sequence; Mut nmer = Mutant amino acid sequence; Exome VAF decile = The variant allele frequency for all mutations broken into decile for more information see "Development of a model for ranking candidate HLA class I neoantigens based upon datasets of known neoepitopes" Gartner et al.; Gene expression decile = The expression of the gene binned into deciles as based upon al genes expression for that sample; Present in RNA-seq data = 1 or 0  for whether mutation was found in RNA sequencing data 1 = yes, 0 = no.

2) a list of up to 6 Class I HLAs formatted to match netMHC input for example HLA-A01:01
 
3) *optionally you can provide the flag --cpus and increase the number of cpus this will speed up the neMHCstab prediction portion by spreading the commands across cpus

Test to see outputs are the same as those given in examples


## Starting from SB_neoantigen_Models 

Ensure that conda environment is Active. directions to create environment are in instructions.txt
```
conda activate mhcflurry-env
```
Test model
```
mkdir test
cd test
python ../src/GenerateScores.py ../examples/nmer_test_input.xlsx HLA-A02:01 HLA-A03:01 HLA-B13:02 HLA-B15:01 HLA-C05:01 HLA-C06:02
```

This will produce a single output excel files in this test directory same as seen in examples directory

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details
