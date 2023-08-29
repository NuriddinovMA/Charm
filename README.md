# Charm
**Ch**romosome re**ar**rangement **m**odeler.
Charm is python-based tool to simulate Hi-C-maps with the user-defined chromosomal rearrangements. This tool allow to simulate CNVs, inversions, translocation, and extra-chromosomal fragments.

## Requirements
1) Python >= 3.7 with numpy
2) [Jucier Tools](https://github.com/aidenlab/juicer)(for dumping the contacts from the existing .hic-files and/or creating the new .hic-files)

## Quick Start
### Unix
charm.sh [-i ini_file] [-S stage] 
* [ini_file]: the path to ini-file containing paths to the working directory, hic-file, unique SV id(s), model paramaters, and others. See the full ini-file description in the [BIG_EXAMPLE.ini](https://github.com/genomech/Charm/blob/main/BIG_EXAMPLE.ini) The short usefull example see in the [TEST_EXAMPLE.ini](https://github.com/genomech/Charm/blob/main/TEST_EXAMPLE.ini)
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
### others OS
python3 scripts/charm_manager.py [-S stage] -i [ini_file]

## Test dataset
TODO: describe the minimal test here

## First run
Modify **TEST_EXAMPLE.ini**:
1) [global] section, "work_dir" key: provide the path to the working directory
TODO: refactor the rest of the manual as above.


2) provide the unique experiment ID in the [global] section, on the "simulation_id" key
3) write to ini-file the path to juicer tools .jar-file in the [global] section, on the "path_to_juicertools" key  
4.5 optional, if path to java not defined, write to ini-file twihe path to java directory in the [global] section, on the "path_to_java_dir" key 
5) write to ini-file the path to **you_chromosome_sizes_file** file in the [global] section, on the "chrom_sizes" key; the file format see below 
6) write to ini-file the path to **you_hic_map.hic** file in the [prerpoceesing] section, on the "path_to_hic" key
7) write to ini-file the path to **you_SVs_list_file** in the [SVs] section, on the "path_to_svs_list" key; the file format see below
8) change other parameters as you wish
9) run Charm

### The chromosome sizes file
This file contains chromosome sizes (see https://github.com/aidenlab/juicer/wiki/Pre). The chromosome names and chromosome sizes must correspond to the chromosome sizes and chromosome names in .hic-file. 
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
  1) An general translocation; the moving of locus 1,000,000-2,000,000 on the chromosome 1 to 7,000,000 position: 
  ```
test	trn	1	0	1000000	!>	1	0	1
test	trn	1	2000000	7000000	>>	1	0	1
test	trn	1	1000000	2000000	>>	1	0	1
test	trn	1	7000000	+	>!	1	0	1
  ```
  2) An tandem duplication of locus 1,000,000-2,000,000 on the chromosome 1
  ```
test	dups	1	0	*2000000*	!>	1	0	1
test	dups	1	1000000	*2000000*	>>	1	*1*	1
test	dups	1	7000000	+	>!	1	0	1

  ```
  or
  ```
test	dups	1	0	*1000000*	!>	1	0	1
test	dups	1	*1000000*	2000000	>>	1	0	*2*
test	dups	1	7000000	+	>!	1	0	1
  ```
  3) An deletion of locus 1,000,000-2,000,000 on the chromosome 1

  ```
test	del	1	0	1000000	!>	1	0	1
test	del	1	1000000	2000000	>>	1	0	0
test	del	1	7000000	+	>!	1	0	1

  ```
 4) An inversion of locus chr1:1,000,000-2,000,000 
  ```
test	inv	1	0	1000000	!>	1	0	1
test	inv	1	1000000	2000000	>>	1	0	*-1*
test	inv	1	7000000	+	>!	1	0	1
  ```
  5) An translocation of locus chr1:1,000,000-2,000,000 to chromosome end with the locus saving on old position and cnv x10 on new position
  ```
test	trnx10	1	0	+	!>	1	0	1
test	trnx10	1	1000000	2000000	>!	1	1	10
  ```
  6) An interchomosome  translocation; the moving of locus 1,000,000-2,000,000 on the chromosome 1 to 7,000,000 position on the chromosome 2: 
  ```
test	trn	1	0	1000000	!>	1	0	1
test	trn	1	7000000	+	>!	1	0	1
test	trn	2	0	700000	>>	2	0	1
test	trn	1	2000000	7000000	>>	*2*	0	1
test	trn	2	7000000	+	>!	2	0	1
  ```

  7) An complex rearrangement: the translocation with cnv x3 of locus 1,000,000-2,000,000 on a new chromosome, the translocation with inversion and cnv x3 of locus 3,000,000-4,000,000 on the new chromosome,  the translocation with cnv x5 of locus 5,000,000-7,500,000 on the new chromosome. The chromosomes 1 is intact.
  ```
 test	compX	1	0	+	!> 1	0	1
 test	compX 1 1000000 2000000 >> chrNew + 1 3
 test	compX 1 3000000 4000000 >> chrNew + 1 -3
 test	compX 1 5000000 7500000 >! chrNew + 1 5
  ```