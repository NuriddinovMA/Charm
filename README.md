# Charm
**Ch**romosome re**ar**rangement **m**odeler.
Charm is python-based tool to simulate Hi-C-maps with the user-defined chromosomal rearrangements. This tool allow to simulate CNVs, inversions, translocation, and extra-chromosomal fragments.

## Requirements
1) Python }= 3.7 with numpy
2) [Jucier Tools](https://github.com/aidenlab/juicer)(for dumping the contacts from the existing .hic-files and/or creating the new .hic-files)

## Test dataset run
```
charm.sh -i testdataset/EXAMPLE.ini
```
As ending, the Charm creates in "testdataset" folder "out" containing the hi-c file with simulated rearrangement named "example.cnv-X.hic"

# Quick Start

## Example tasks
### 1) The establish of reference database without consequent simulations:
```
charm -i testdataset/EXAMPLE.ini *-S pre*
```
or
```
charm -i testdataset/example_PRE.ini
```
As ending, the Charm creates in "testdataset" the folders "pre/TEST" with: 

- the files "TEST.5000.stat", "TEST.5000.binCov" and the folder "TEST.5000" containing files "TEST.5000.{chr1}.{chr2}.allCon"

- the files "TEST.50000.stat", "TEST.50000.binCov" and the folder "TEST.50000" containing files "TEST.50000.{chr1}.{chr2}.allCon"

- the files "pab.TEST.{resolution_pab}.stat", "pab.TEST.{resolution_pab}.binCov" and the folder "pab.TEST.{resolution_pab}" containing files "Tpab.TEST.{resolution_pab}.{chr1}.{chr2}.allCon"

where {chr1} and {chr2} are the chromosome names and {resolution_pab} is the resolution from [global] section, "resolution_pab" key.

### 2) The establish of database of randomized wild-type contacts, IF the reference database was performed:
```
charm -i testdataset/EXAMPLE.ini -S wt
```
or
```
charm -i testdataset/example_WT.ini
```
As ending, the Charm creates in "testdataset" the folders "wt/TEST.cov_mult_f1/841160/0/" with files named like "TEST.cov_mult_f1.0.{chr1}.{chr2}.allCon" and the folders "wt/TEST.cov_mult_f1/841160/1/" with files named like "TEST.cov_mult_f1.1.{chr1}.{chr2}.allCon"

### 3) The simulation of *heterozygous* mutation, IF the reference database and the pseudoreplicas were performed:
```
charm -i testdataset/example_HETEROZYGOUS.ini -S SVs+
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "heterozygous.del.hic"

### 4) The simulation of *homozygous* mutation, IF the reference database and the pseudoreplicas were performed:
```
charm -i testdataset/example_HOMOZYGOUS.ini -S SVs+
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "homozygous.del.hic"

### 5) The simulation of *mutant* genome, IF the reference database:
```
charm -i testdataset/example_MUTANT.ini -S SVs+
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "in_mut.cnv-X.hic".
This hi-c file will contain the new chromosome "1x".

### 6) The building of wild-type contact map, IF the reference database and the pseudoreplicas were performed:
```
charm -i testdataset/example_REPLICAS.ini -S hic
```
As ending, the Charm creates in "testdataset" the folder "out" containing the hi-c file with simulated rearrangement named "replicas.TEST.cov_mult_f1.hic".

## The your task: step by step.
1) Create your file with the description of rearrangment (see [The SVs description](https://github.com/NuriddinovMA/Charm#the-svs-description))
2) Duplicate the any examplified ini-file accordingly your tasks (see [example tasks](https://github.com/NuriddinovMA/Charm#example-tasks) and change it:
  * [global] section:
    - "work_dir" - the path to the your work directory;
    - "chrom_sizes" - the path to the file with the chromosome sizes of reference genome;
    - "path_to_juicertools " - the path to the juicertools jar file;
    - "noised" - "True" if the reference HI-C is whole genomic, "False" if the reference Hi-C is enriched like promoter-capture;
    - "simulation_id" - the preferred name of simulations.
  * [preprocessing] section
    - "path_to_hic" - the path to the hic file with the reference contact map;
  * [SVs] section
    - "path_to_svs_list" - the path to the your file with the description of rearrangments;
    - "rearrangment_id" - the unique id of simulated rearrangment;
  * [simulation] section
    - "contact_count" - the summ of contacts on simulated hi-c map
    - "predict_null_contacts" - use or "cov_mult_f"/"cov_sq_f"/"cov_mult_f1"/"cov_sq_f1" for whole genomic Hi-C, "cov_mixed_f"/"cov_mixsq_f"/"cov_mixed_f1"/"cov_missq_f1" for enriched Hi-C
  * [hic]
    - "simulation_id" - the unique name of resulted simulation
    - "format" - "hic" for juicer tools hic-map, "pre" for the [pre-file] (https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format),
      "short" for the [extra short] (https://github.com/aidenlab/juicer/wiki/Pre#extra-short-format-dev) pre-file,
      "pre.gz" and "short.gz" - the gzipped output.
    - "hic_resolutions" - the list of Hi-C map bin sizes; the minimal bin size must be equal to [global] "resolution" or higher. 
3) Run charm

### others OS
```
python3 scripts/charm_manager.py -i [ini_file] [-S stage]
```
### The chromosome sizes file
This file contains chromosome sizes ([example](https://github.com/NuriddinovMA/Charm/blob/main/testdataset/data/test.chrom.sizes)). The chromosome names and chromosome sizes must correspond to the chromosome sizes and chromosome names in .hic-file. 
File format (see the example "test.chr.sizes")
```
<chromosome name> <chromosome size bp>
```

### The SVs description 
To simulate SVs, Charm requires the description of rearrangement in the following format (also see the example "test.svs_list.txt" in the testdataset folder):
```
<reference genome id> <mutant genome uniq id> <chromosome> <coordinate chromosome block start> <coordinate of chromosome block end> <the indicator> <new chromosome> <copy number of locus on OLD position> <copy number of locus on NEW position>
```
The indicator variants:
  - Use "->" for the plain SVs, this indicator designs the start and the end of SVs description
  - Use "!>" for the start of the description of complex SVs
  - Use ">>" for the continuation of description SVs
  - Use ">!" for the end of the description of complex SVs, the all lines between "!>" and ">!" are processed by Charm as one SV.
 
Examples:

*(1)* An general translocation; the moving of locus 1:1Mb-2Mb to 7Mb: 
```
test	trn	1	0	1000000	!>	1	0	1
test	trn	1	2000000	7000000	>>	1	0	1
test	trn	1	1000000	2000000	>>	1	0	1
test	trn	1	7000000	+	>!	1	0	1
```
*(2)* An tandem duplication of locus 1:1Mb-2Mb
```
test	dups	1	0	2000000	!>	1	0	1
test	dups	1	1000000	2000000	>>	1	1	1
test	dups	1	7000000	+	>!	1	0	1
```
  or
```
test	dups	1	0	1000000	!>	1	0	1
test	dups	1	1000000	2000000	>>	1	0	2
test	dups	1	7000000	+	>!	1	0	1
```
*(3)* An deletion of locus 1:1Mb-2Mb

```
test	del	1	0	1000000	!>	1	0	1
test	del	1	1000000	2000000	>>	1	0	0
test	del	1	7000000	+	>!	1	0	1
```
*(4)* An inversion of locus 1:1Mb-2Mb
```
test	inv	1	0	1000000	!>	1	0	1
test	inv	1	1000000	2000000	>>	1	0	-1
test	inv	1	7000000	+	>!	1	0	1
```
*(5)* An translocation of locus 1:1,000,000-2,000,000 to chromosome end with the locus saving on old position and cnv x10 on new position
```
test	trnx10	1	0	+	!>	1	0	1
test	trnx10	1	1000000	2000000	>!	1	1	10
```
*(6)* An interchomosome  translocation; the moving of locus 1:1Mb-2Mb to the chromosome 2:7Mb: 
```
test	trn	1	0	1000000	!>	1	0	1
test	trn	1	7000000	+	>!	1	0	1
test	trn	2	0	700000	>>	2	0	1
test	trn	1	2000000	7000000	>>	2	0	1
test	trn	2	7000000	+	>!	2	0	1
```

*(7)* An complex rearrangement: the translocation with cnv x3 of locus 1:1Mb-2Mb on a new chromosome, the translocation with inversion and cnv x3 of locus 1:3Mb-4Mb on the new chromosome,
the translocation with cnv x5 of locus 1:5Mb-7,5Mb on the new chromosome. The chromosomes 1 is intact.
```
test	compX	1	0	+	!> 1	0	1
test	compX 1 1000000 2000000 >> chrNew + 1 3
test	compX 1 3000000 4000000 >> chrNew + 1 -3
test	compX 1 5000000 7500000 >! chrNew + 1 5
```
*(8)* Several rearrangements: a translocation from chromosome 1:1Mb-2Mb to 2:7Mb, tandem duplication of 1:7Mb-8Mb and deletion of 2:1Mb-2Mb
```
test	several	1	0	1000000	!>	1	0	1
test	several	1	2000000	7000000	>>	1	0	1
test	several	1	7000000	8000000	>>	1	0	2
test	several	1	8000000	+	>>	1	0	1
test	several	2	0	1000000	>>	2	0	1
test	several	2	1000000	2000000	>>	2	0	0
test	several	2	2000000	7000000	>>	2	0	0
test	several	1	2000000	7000000	>>	2	0	1
test	several	2	7000000	+	>!	2	0	1
```

## Advanced description
charm [-i ini_file] [-S stage] 
* [ini_file]: the path to ini-file containing paths to the working directory, hic-file, unique SV id(s), model paramaters, and others. See the full ini-file description in the [BIG_EXAMPLE.ini](https://github.com/NuriddinovMA/Charm/blob/main/BIG_EXAMPLE.ini)
The short useful example see in the [EXAMPLE.ini](https://github.com/NuriddinovMA/Charm/blob/main/EXAMPLE.ini)
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
