### Divergence

The nonCoveredRegion.py calculates the amount of sequence present/absent in the alignments and reports the total and percentage values. It needs a custom file, derived from maf file. 
For example, if the maf file has the name "pairwisesorangvsborangchr10.maf", and has two entries "sorang.chr10" and borang.chr10", to extract the custom file for the sorang, use the command:

`cat pairwisesorangvsborangchr10.maf | grep "s sorang.chr10" | python nonCoveredRegion.py - > out.dat`
