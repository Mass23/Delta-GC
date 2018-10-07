# Delta-GC

A python script that allows to calculate the sliding-window pairwise difference in GC-content between two aligned genomes in mauve (.xmfa) format. The result is shown on a manhattan plot and given in a text file.

Manhattan plot of delta-GC (the grey line represents the absolute value of the difference):

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_1.png)

The different color represents a scaffold of interest.

Identity (computed as the proportion of shared allele between the two genomes):

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_2.png)

Uncertainty (computed as the propotion of N's in the two sequences):

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_3.png)

Linear regression of GC1 on GC2:

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_4.png)
