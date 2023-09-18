# Charm
**Ch**romosome re**ar**rangement **m**odeler.
Charm is python-based tool to simulate Hi-C-maps with the user-defined chromosomal rearrangements. This tool allows to model CNVs, inversions, translocation, and extra-chromosomal fragments.

## Requirements
1) Python >= 3.7 with numpy
2) [Juicer Tools](https://github.com/aidenlab/juicer) (for dumping the contacts from the existing .hic-files and/or creating the new .hic-files)
3) Java (the Juicer Tools requirement)
## Test dataset run
```
python3 charm.py -i testdataset/EXAMPLE.ini
```
As ending, the Charm creates in "testdataset" folder "out" containing the hi-c file with simulated rearrangement named "example.cnv-X.hic"

# Quick Start

## Example tasks
### 1) The establish of reference database without consequent simulations:
```
python3 charm.py -i testdataset/EXAMPLE.ini -S pre
```
or
```
python3 charm.py -i testdataset/example_PRE.ini
```
As ending, the Charm creates in "testdataset" the folders "pre/TEST" with: 

- the files "TEST.5000.stat", "TEST.5000.binCov" and the folder "TEST.5000" containing files "TEST.5000.{chr1}.{chr2}.allCon"

- the files "TEST.50000.stat", "TEST.50000.binCov" and the folder "TEST.50000" containing files "TEST.50000.{chr1}.{chr2}.allCon"

- the files "pab.TEST.{resolution_pab}.stat", "pab.TEST.{resolution_pab}.binCov" and the folder "pab.TEST.{resolution_pab}" containing files "Tpab.TEST.{resolution_pab}.{chr1}.{chr2}.allCon"

where {chr1} and {chr2} are the chromosome names and {resolution_pab} is the resolution from [global] section, "resolution_pab" key.

### 2) The establish of database of randomized wild-type contacts, IF the reference database was performed:
```
python3 charm.py -i testdataset/EXAMPLE.ini -S wt
```
or
```
python3 charm.py -i testdataset/example_WT.ini
```
As ending, the Charm creates in "testdataset" the folders "wt/TEST.cov_mult_f1/841160/0/" with files named like "TEST.cov_mult_f1.0.{chr1}.{chr2}.allCon" and the folders "wt/TEST.cov_mult_f1/841160/1/" with files named like "TEST.cov_mult_f1.1.{chr1}.{chr2}.allCon"

### 3) The simulation of *heterozygous* mutation, IF the reference database and the pseudoreplicas were performed:
```
python3 charm.py -i testdataset/example_HETEROZYGOUS.ini -S SVs+
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "heterozygous.del.hic"

### 4) The simulation of *homozygous* mutation, IF the reference database and the pseudoreplicas were performed:
```
python3 charm.py -i testdataset/example_HOMOZYGOUS.ini -S SVs+
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "homozygous.del.hic"

### 5) The simulation of *mutant* genome, IF the reference database:
```
python3 charm.py -i testdataset/example_MUTANT.ini -S SVs+
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "in_mut.cnv-X.hic".
This hi-c file will contain the new chromosome "1x".

### 6) The building of wild-type contact map, IF the reference database and the pseudoreplicas were performed:
```
python3 charm.py -i testdataset/example_REPLICAS.ini -S hic
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "replicas.TEST.cov_mult_f1.hic".

## The your task: step by step.
1) Create your file with the description of rearrangment (see [The SVs description](https://github.com/NuriddinovMA/Charm#the-svs-description))
2) Duplicate the any examplified ini-file accordingly your tasks (see [example tasks](https://github.com/NuriddinovMA/Charm#example-tasks) and change it:
  * [global] section:
    - "work_dir" - the path to the your work directory;
    - "chrom_sizes" - the path to the file with the chromosome sizes of reference genome;
    - "path_to_juicertools " - the path to the juicertools jar file;
    - "one_as_null" - "True" contacts == 1 are processed as 0 \("True" should be used for the whole genomic Hi-C, and "False" should be used for the enriched Hi-C, like promoter-capture\);
    - "simulation_id" - the preferred name of simulations.
  * [preprocessing] section
    - "path_to_hic" - the path to the hic file with the reference contact map;
  * [SVs] section
    - "path_to_svs_list" - the path to the your file with the description of rearrangments;
    - "rearrangment_id" - the unique id of simulated rearrangment from SVs list;
  * [simulation] section
    - "contact_count" - the summ of contacts on simulated hi-c map
    - "predict_null_contacts" - use or "cov_mult_f"/"cov_sq_f"/"cov_mult_f1"/"cov_sq_f1" for whole genomic Hi-C, "cov_mixed_f"/"cov_mixsq_f"/"cov_mixed_f1"/"cov_missq_f1" for enriched Hi-C
  * [hic]
    - "simulation_id" - the unique name of resulted simulation
    - "format" - "hic" for juicer tools hic-map, "pre" for the [pre-file] (https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format),
      "short" for the [extra short] (https://github.com/aidenlab/juicer/wiki/Pre#extra-short-format-dev) pre-file,
      "pre.gz" and "short.gz" - the gzipped output.
    - "hic_resolutions" - the list of Hi-C map bin sizes; the minimal bin size must be equal to [global] "resolution" or higher. 
3) python3 charm.py -i [ini-file] -S [step]

### The chromosome sizes file
This file contains chromosome sizes ([example](https://github.com/NuriddinovMA/Charm/blob/main/testdataset/data/test.chrom.sizes)). The chromosome names and chromosome sizes must correspond to the chromosome sizes and chromosome names in .hic-file. 
File format (see the example "test.chr.sizes")
```
<chromosome name> <chromosome size bp>
```

### The SVs description 
To simulate SVs, Charm requires the file with a description of rearrangement. It shows which fragments of reference chromosomes compose the rearranged chromosome and their order.

![grafical example](https://github.com/NuriddinovMA/Charm/blob/main/description.png)

The format of the SVs list file (also see the example "test.svs_list.txt" in the testdataset folder):
```
<reference genome id> <rearrangment id> <chromosome> <coordinate chromosome block start> <coordinate of chromosome block end> <indicator> <new chromosome> <copy number of locus on OLD position> <copy number of locus on NEW position>
```
The  \<reference genome id\> and  the \<rearrangment id\> are any names to description the reference genome and the modeled rearrangements. Every model must be named uniquelly. 

The \<chromosome\> is the name of the reference genome chromosome involved in the rearrangement.

The \<coordinate chromosome block start\> and the\<coordinate of chromosome block end\> are the coordinates of breakpoints in reference genome. The "+" should be used in \<coordinate of chromosome block end\> column as the symbol of the chromosome end.

The \<indicator\> variants:
  - Use "->" for the plain SVs, this indicator designs the start and the end of SVs description
  - Use "!>" for the start of the description of complex SVs
  - Use ">>" for the continuation of description SVs
  - Use ">!" for the end of the description of complex SVs, the all lines between "!>" and ">!" are processed by Charm as one SV.

All lines between the "!>" and "<!" indicators must have the same \<rearrangment id\>.

The \<new chromosome\> is the name of the simulated chromosome resulting from the rearrangement. This name can match with the \<chromosome\> or can be novel.

The values in \<copy number of the locus on OLD position\> must be 0 or 1. The 1 can be used only for CNV simulations.

The values in \<copy number of the locus on NEW position\> can be any; the negative values correspond to the inversion; the "0" corresponds to the deletion if \<copy number of the locus on OLD position\> is "0", too.

Examples:

*(1)* An general translocation; the moving of locus 1:1Mb-2Mb to 7Mb: 
```
test	trn	chr1	0	1000000	!>	chr1	0	1
test	trn	chr1	2000000	7000000	>>	chr1	0	1
test	trn	chr1	1000000	2000000	>>	chr1	0	1
test	trn	chr1	7000000	+	>!	chr1	0	1
```
*(2)* An tandem duplication of locus 1:1Mb-2Mb
```
test	dups	chr1	0	2000000	!>	chr1	0	1
test	dups	chr1	1000000	2000000	>>	chr1	1	1
test	dups	chr1	2000000	+	>!	chr1	0	1
```
  or
```
test	dups	chr1	0	1000000	!>	chr1	0	1
test	dups	chr1	1000000	2000000	>>	chr1	0	2
test	dups	chr1	2000000	+	>!	chr1	0	1
```
*(3)* An deletion of locus 1:1Mb-2Mb

```
test	del	chr1	0	1000000	!>	chr1	0	1
test	del	chr1	1000000	2000000	>>	chr1	0	0
test	del	chr1	2000000	+	>!	chr1	0	1
```
*(4)* An inversion of locus 1:1Mb-2Mb
```
test	inv	chr1	0	1000000	!>	chr1	0	1
test	inv	chr1	1000000	2000000	>>	chr1	0	-1
test	inv	chr1	2000000	+	>!	chr1	0	1
```
*(5)* An translocation of locus 1:1,000,000-2,000,000 to chromosome end with the locus saving on old position and cnv x10 on new position
```
test	trnx10	chr1	0	+	!>	chr1	0	1
test	trnx10	chr1	1000000	2000000	>!	chr1	1	10
```
*(6)* An interchomosome  translocation; the moving of locus 1:1Mb-2Mb to the chromosome 2:7Mb: 
```
test	trn	chr1	0	1000000	!>	chr1	0	1
test	trn	chr1	2000000	+	>>	chr1	0	1
test	trn	chr2	0	700000	>>	chr2	0	1
test	trn	chr1	1000000	2000000	>>	chr2	0	1
test	trn	chr2	7000000	+	>!	chr2	0	1
```

*(7)* An complex rearrangement: the translocation with cnv x3 of locus 1:1Mb-2Mb on a new chromosome, the translocation with inversion and cnv x3 of locus 1:3Mb-4Mb on the new chromosome,
the translocation with cnv x5 of locus 1:5Mb-7,5Mb on the new chromosome. The chromosomes 1 is intact.
```
test	compX	chr1	0	+	!> 1	0	1
test	compX chr1 1000000 2000000 >> chrNew + 1 3
test	compX chr1 3000000 4000000 >> chrNew + 1 -3
test	compX chr1 5000000 7500000 >! chrNew + 1 5
```
*(8)* Several rearrangements: a translocation from chromosome 1:1Mb-2Mb to 2:7Mb, tandem duplication of 1:7Mb-8Mb and deletion of 2:1Mb-2Mb
```
test	several	chr1	0	1000000	!>	chr1	0	1
test	several	chr1	2000000	7000000	>>	chr1	0	1
test	several	chr1	7000000	8000000	>>	chr1	0	2
test	several	chr1	8000000	+	>>	chr1	0	1
test	several	chr2	0	1000000	>>	chr2	0	1
test	several	chr2	1000000	2000000	>>	chr2	0	0
test	several	chr2	2000000	7000000	>>	chr2	0	0
test	several	chr1	2000000	7000000	>>	chr2	0	1
test	several	chr2	7000000	+	>!	chr2	0	1
```

## Advanced description
charm [-i ini_file] [-S stage] 
* [ini_file]: the path to ini-file containing paths to the working directory, hic-file, unique SV id(s), model paramaters, and others. See the full ini-file description in the [BIG_EXAMPLE.ini](https://github.com/NuriddinovMA/Charm/blob/main/BIG_EXAMPLE.ini). See the common description of ini-file in python module [confiparser](https://docs.python.org/3/library/configparser.html), class configparser.ExtendedInterpolation.
* [stage]: optional, must be one of "pre+","SVs+","sim+","lift+","wt+","hic" (default "pre+")
  - "pre+" is the default parameter, from hic-file with the wild-type data to hic-file with simulated SVs
  - Use "SVs+" when the database with contact statistics exists, but the file with rearrangement description has not yet been created 
  - Use "sim+" when the rearrangement description and database are done, but the contacts of the mutant genome are not yet been simulated
  - Use "lift+" when the contacts of the mutant genome are simulated ("in_mut.[simulation_id]" files created), but these contacts are not yet lifted on the reference genome
  - Use "wt+" when the simulation is fully processed ("in_ref.[simulation_id]" created), but wild-type replicas are not
  - Use "hic" when all previous stages are successfully completed, but .hic-file is not created
* advanced settings for [stage]:
  - Use "pre" to create the database with contact statistics only
  - Use "SVs" to create the rearrangement description files only; in this mode, Charm can process the file containing any number of independent SVs 
  - Use "sim" to only simulate contacts within *mutated* genome based on provided database
  - Use "lift" to liftover contacts from provided defined mutant genome
  - Use "wt" to simulate wild-type replicas
