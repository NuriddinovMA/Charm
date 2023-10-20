# The dataset
The simulated dataset of structural variations consist of [whole-genome Hi-C models](https://genedev.bionet.nsc.ru/ftp/by_Project/Charm/wgHi-C/) and [exome-capture Hi-C models](https://genedev.bionet.nsc.ru/ftp/by_Project/Charm/ecHi-C/).
Both the wgHi-C dataset and the ecHi-C dataset are built on the same rearrangements described in [simulation.txt](simulations.txt), where each line represents a unique SV
'''
<id>  <type SVs>  <chrm>  <start>  <end>  <note>
'''
- **id** is unique id number of rearrangment included in model file name;
- **type SVs** is type of structural variantion;
- **chrm**, **start**, **end** are genome coordianate of rearranged locus;
- **note** includes additional information for int**er**- and int**ra**chromosomal translocations: the genome coordinates of insertion and trancsloated locus orientation.

The every model in dataset is a gzipped tab separated file that contains, on each line
'''
<chrm1>  <pos1>  <chrm2>  <pos2>  <score>
'''
- **chrm1**, **pos1**, **chrm2**, **pos2** are genome coordinate of Hi-C contact;
- **score** is the contact count.

**Attention!** All SVs simalated as **hetero**zygous!
