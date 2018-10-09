# Delta-GC

# To Add:
- Use coordinates instead of scaffold number to define the region of interest

A python script that allows to calculate the sliding-window pairwise difference in GC-content between two aligned genomes in mauve (.xmfa) format. The result is shown on manhattan plots and given in a tab-separated file.

# Manhattan plot of delta-GC:
- Computed for each sliding window as: Genome 1 GC content - Genome 2 GC content
- The grey line represents the absolute value of the difference.
- The different colors represents a region of interest.

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_1.png)

# Identity 
- Computed as the proportion of shared allele between the two genomes.

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_2.png)

# Uncertainty
- Computed as the propotion of N's in the two sequences.

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_3.png)

# Linear regression of GC1 on GC2

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_4.png)

# Delta-GC distribution plots

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_4.png)

# Delta-GC ANOVA and means
