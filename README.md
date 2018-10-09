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

![alt text](https://github.com/Mass23/Delta-GC/blob/master/Figure_5.png)

# Delta-GC ANOVA and means

One-way ANOVA results:

|             | sum_sq      |df           |F          |PR(>F)       |
|-------------|-------------|-------------|-----------|-------------|
|of_interest  |0.000585     |1.0          |28.573769  |9.233950e-08 |
|Residual     |0.189652     |9271.0       |NaN        |NaN          |

Mean of Delta-GC: 
 - Scaffold of interest:  0.0009492675285102437
 - Genome:  -0.00011819850083390379
 
Mean of Identity: 
 - Scaffold of interest:  0.980309173841763
 - Genome:  0.9918252714834066
 
Mean of Uncertainty: 
 - Scaffold of interest:  0.048055825688073366
 - Genome:  0.013321694546287791
