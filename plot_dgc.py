# Author: Massimo Bourquin
# Contact: massimo.bourquin@unil.ch
# Arguments
# python3 Delta-GC.py"M_P_alignment.xmfa ['genom1', 'genome2'] window-size threshold-of-uncertainty-to-drop-window [region of : interest]

from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics

def plotting(file, genomename1, genomename2, window_size):

# Plotting:
    data = pd.read_csv('delta_gc.csv', delimiter="\t")

    mask_of_interest = data['of_interest'].eq(1)
    mask_genome = data['of_interest'].eq(0)

    of_interest = data[mask_of_interest].dropna()
    genome = data[mask_genome].dropna()

# Scaffolds coordinates
    scaff_coord = []
    actual_scaff = ''
    for i in range(0, data['window_number'].count()):
        if data['scaffold'][i] == actual_scaff:
            continue
        else:
            scaff_coord.append(data['window_number'][i])
            actual_scaff = data['scaffold'][i]
            continue

    print(scaff_coord)
# Plot delta-GC
    plt.figure(1)
    plt.title("Delta-GC")
    sns.set_style('white')
    sns.set_context("paper")

    # Vertical scaffold lines
    for i in range(0, len(scaff_coord)):
        plt.axvline(x = scaff_coord[i], color = "lightgrey")

    # Horizontal line at 0
    plt.plot((0, data['abs_dgc'].count()), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 0.8)

    # Plot absolute value of dgc
    plt.plot(data['window_number'],data['abs_dgc'], '-', color = 'darkgrey')

    # Plot the dgc    
    plt.plot(genome['window_number'], genome['delta_gc'], 'o', color = 'salmon')
    plt.plot(of_interest['window_number'], of_interest['delta_gc'], 'o', color = 'steelblue')

    # Axis labels
    plt.xlabel('Genome windows (size: ' + str(window_size) + 'bp.)')
    plt.ylabel('ΔGC: ' + genomename1 + ' - ' + genomename2)

    plt.show()
    plt.close()


# Plot Identity
    plt.figure(2)
    plt.title("Identity")
    sns.set_style('white')
    sns.set_context("paper")

    # Vertical scaffold lines
    for i in range(0, len(scaff_coord)):
        plt.axvline(x = scaff_coord[i], color = "lightgrey")

    #Plot an horizontal line at 0, -1 and 1
    plt.plot((0, data['identity'].count()), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 1)
    plt.plot((0, data['identity'].count()), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(genome['window_number'], genome['identity'],  'o', color = 'salmon')
    plt.plot(of_interest['window_number'], of_interest['identity'], 'o', color = 'steelblue')

    # Axis labels
    plt.ylabel('Identity', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)
    plt.show()
    plt.close()


# Plot Uncertainty
    plt.figure(3)
    plt.title("Uncertainty")
    sns.set_style('white')
    sns.set_context("paper")

    # Vertical scaffold lines
    for i in range(0, len(scaff_coord)):
        plt.axvline(x = scaff_coord[i], color = "lightgrey")

    #Plot an horizontal line at 0, -1 and 1
    plt.plot((1, data['uncertainty'].count()), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 1)
    plt.plot((1, data['uncertainty'].count()), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(genome['window_number'], genome['uncertainty'],  'o', color = 'salmon')
    plt.plot(of_interest['window_number'], of_interest['uncertainty'], 'o', color = 'steelblue')
    
    # Axis labels
    plt.ylabel('Uncertainty', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)
    
    # Output
    plt.show()
    plt.close()

    print("Output Manhattan plots Done!")

# Plot GC1/GC2
    plt.figure(4)
    plt.title("GC1 / GC2")
    sns.set_style('whitegrid')

    plt.subplot(211)
    sns.set_style('whitegrid')
    sns.set_context("paper")
    sns.kdeplot(genome['gc1'], genome['gc2'], cmap="Reds", shade=True, shade_lowest=False)
    plt.xlim(0.2,0.6)
    plt.ylim(0.2,0.6)

    plt.xlabel(genomename1 + " GC content", size = 15)
    plt.ylabel(genomename2 + " GC content", size = 15)

    plt.subplot(212)
    sns.set_style('whitegrid')
    sns.set_context("paper")
    sns.kdeplot(of_interest['gc1'], of_interest['gc2'], cmap="Blues", shade=True, shade_lowest=False)
    plt.ylim(0.2,0.6)
    plt.xlim(0.2,0.6)

    plt.xlabel(genomename1 + " GC content", size = 15)
    plt.ylabel(genomename2 + " GC content", size = 15)

    plt.show()
    plt.close()

# Plot distributions
    from statsmodels.formula.api import ols
    import statsmodels.api as sm

    plt.figure(5)
    plt.title("Density plot")
    sns.set_style('whitegrid')
    sns.set_context("paper")

    sns.distplot(genome['delta_gc'], color="salmon", hist=False)
    sns.distplot(of_interest['delta_gc'], color="steelblue", hist=False)

    plt.show()
    plt.close()

# ANOVA
    lm = ols('delta_gc ~ of_interest', data=data).fit()
    table = sm.stats.anova_lm(lm, typ=2) # Type 2 ANOVA DataFrame
    print("One-way ANOVA results:")
    print(table)

    print("Mean of Delta-GC: ")
    print(" - Scaffold of interest: ", of_interest['delta_gc'].mean())
    print(" - Genome: ", genome['delta_gc'].mean())

    print("Mean of Identity: ")
    print(" - Scaffold of interest: ", of_interest['identity'].mean())
    print(" - Genome: ", genome['identity'].mean())

    print("Mean of Uncertainty: ")
    print(" - Scaffold of interest: ", of_interest['uncertainty'].mean())
    print(" - Genome: ", genome['uncertainty'].mean())



plotting('delta_gc.csv', 'FselM', 'FselP', 50000)
