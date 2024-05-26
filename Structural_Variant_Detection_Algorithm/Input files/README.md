# This folder contains the input files and a an example output.vcf file that how output will look like

### Sample Output

Here is an example of what the VCF output might look like:

```
##fileformat=VCFv4.2
#CHROM  POS     ID      REF     ALT             QUAL    FILTER  INFO
1       10000   .       N       <DEL:1500>      .       PASS    SVTYPE=DEL;END=11500;SVLEN=-1500
1       20000   .       N       <INV:500>       .       PASS    SVTYPE=INV;END=20500;SVLEN=500
2       30000   .       N       <DUP:300>       .       PASS    SVTYPE=DUP;END=30300;SVLEN=300
1       40000   .       N       <TRA>           .       PASS    SVTYPE=TRA;CHR2=2;POS2=50000