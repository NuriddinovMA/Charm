# Charm
*Ch*romosome re*ar*rangement *m*odeler.
There ara python the based tool to simulate of Hi-C-map with the predetermened chromosomal rearrangments. This tools allow simulate such structural variations as CNVs, inversions, translocation and the extra-chromosome creation.

# Requirements
1) Python >= 3.7 with numpy module
2) Jucier Tools https://github.com/aidenlab/juicer (for dumping the contacts from the .hic-files and creation the new .hic-files)

# Quick Start

charm [-i ini_file] [-S stage] 
* [ini_file]: the path to ini-file contains all needed paramaters to simulate SVs, the ini-file description see in test.ini
* [stage]: optitional, must be one from "pre+","SVs+","sim+","lift+","wt+","hic" (default "pre+")
  - "pre+" is default parameter, from hic-file with wild_type to hic-file with SVs.
  - Use "SVs+" when database with contact statistics are done, but the file with rearrangment description has not yet been created 
  - Use "sim+" when the rearrangment description and database are done, but the contacts of mutant genome are not yet been simulated
  - Use "lift+" when the contacts of mutant genome are simulated, but this contacts are not yet lifted on the reference genome
  - Use "wt+" when the simulation is fully processed, but wild-type replicas are not done
  - Use "hic" when the all previously stages are successfully conplited, but .hic-file is not created
* advanced settings for [stage]:
  - Use "pre" to creation the database with contact statistics only
  - Use "SVs" to creation the rearrangment description files only; in this mode Charm can process the file containing any number of independent SVs 
  - Use "sim" to only simulate contacts within *mutanted* genome based on handly defined database
  - Use "lift" to lifover contacts from handly defined mutant genome
  - Use "wt" to simulated wild-type replicas

# The SVs description 
To simulated SVs, Charm requers the description of rearrangment given in the correct format
File format (see example "test_rear.txt" in the testdataset folder)
```
<reference genome id> <mutant genome uniq id> <chromosome> <coordinate of rearrangmnet locus start> <coordinate of rearrangmnet locus end> <the indicator> <new chromosome> <coordinate of new position of locus> <copy number of locus on OLD position> <copy number of locus on NEW position>
```
The indicator variants:
  - Use "->" for the plain SVs, this indicator designes the start and the end of SVs description
  - Use "!>" for the start of description of complicated SVs
  - Use ">>" for the continuation of description SVs
  - Use ">!" for the end of description of complicated SVs, the all lines between "!>" and ">!" are procesed by Charm as one SV.
 
Examples:
  1) the general translocation; the moving of locus chr1:1,000,000-2,000,000 to chr1:7,000,000 position: 
  ```
  hg19 my_translocation chr1 1000000 2000000 -> chr1 7000000 0 1
  ```
  2) the tandem duplication of locus chr1:1,000,000-2,000,000
  ```
  hg19 my_duplication1 chr1 1000000 2000000 -> chr1 2000000 1 1
  ```
  or
  ```
  hg19 my_duplication2 chr1 1000000 2000000 -> chr1 2000000 0 2
  ```
  3) the deletion of locus chr1:1,000,000-2,000,000 
  ```
  hg19 my_deletion chr1 1000000 2000000 -> chr1 1000000 0 0
  ```
  4) the cnv, the repeating locus chr1:1,000,000-2,000,000  10 times
  ```
  hg19 my_cnv10 chr1 1000000 2000000 -> chr1 1000000 0 10
  ```
  5) the inversion of locus chr1:1,000,000-2,000,000 
  ```
  hg19 my_inversion chr1 1000000 2000000 -> chr1 2000000 0 *-*1
  ```
  6) the translocation of locus chr1:1,000,000-2,000,000 to chromosome end with the locus saving on old position and cnv x10 on new position
  ```
  hg19 my_translocation2 chr1 1000000 2000000 -> chr1 + 1 10
  ```
  7) the complicated rearrangment: the translocation and inversion within one genome
  ```
  hg19 my_complicated2 chr1 1000000 2000000 !> chr1 7000000 0 1
  hg19 my_complicated2 chr1 6500000 7500000 >! chr1 6500000 0 -1
  ```
  8) the very complicated rearrangment: the translocation with cnv and creation of new chromosome
  ```
  hg19 my_complicated3 chr1 1000000 2000000 !> chrNew + 1 3
  hg19 my_complicated3 chr1 3000000 4000000 >> chrNew + 1 -3
  hg19 my_complicated3 chr1 5000000 7500000 >! chrNew + 1 -3
  ```
 
