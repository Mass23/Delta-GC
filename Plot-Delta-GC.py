#!/usr/bin/env python3
# Author: Massimo Bourquin
# Contact: massimo.bourquin@unil.ch
# Arguments: list is for the region of interest (colored in blue)
# python3 Plot-Delta-GC.py Delta-GC.csv beginning:end Genome_name_1 Genome_name_2

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Parse arguments
args = sys.argv
data_file = args[1]
coord_list = args[2]
genomename1 = args[3]
genomename2 = args[4]
beginning = int(coord_list.split(':')[0])
end = int(coord_list.split(':')[1])

# Load data
data = pd.read_csv(data_file, delimiter="\t")
data = data.sort_values(by='Beginning')

data = data.dropna()

data = data[(data['Uncertainty1'] <= 0.3) & (data['Uncertainty2'] <= 0.3)]

of_interest = data[(data['Beginning'] >= beginning) & (data['End'] <= end)]
genome = data[(data['Beginning'] <= beginning) | (data['End'] >= end)]

##############################################################################
# Plot
plt.figure("Delta GC")
sns.set_style('white')
sns.set_context("paper")

# Vertical scaffold lines
scaff_coord = []
scaff_name = []
last_scaffold = ''
for index, row in data.iterrows():
    current_scaffold = row['Scaffold']
    if current_scaffold == last_scaffold:
        continue
    else:
        plt.axvline(row['Beginning'], color = "lightgrey", lw = 0.4)
        scaff_coord.append(row['Beginning'])
        scaff_name.append(int(row['Scaffold']))
        last_scaffold = current_scaffold

# Horizontal line at 0
plt.plot((0, data['Beginning'].max()), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 0.8)

# Plot absolute value of dgc
plt.plot(data['Beginning'], data['Absolute_Delta_GC'], '-', color = 'darkgrey')

# Plot the dgc
plt.plot(genome['Beginning'], genome['Delta_GC'], 'o', color = 'salmon')
plt.plot(of_interest['Beginning'], of_interest['Delta_GC'], 'o', color = 'steelblue')

# Axis labels
plt.tick_params(axis='x', which='major', pad=10)
plt.tick_params(axis='y', which='major', pad=10)
plt.tick_params(axis='both', which='major', pad=10)

plt.xticks(scaff_coord, scaff_name)
plt.xlabel('Scaffolds', size = 15)
plt.ylabel('ΔGC: ' + genomename1 + ' - ' + genomename2, size = 15)

plt.show()
plt.close()

# GC1 / GC2
plt.figure('GC1 / GC2')

plt.plot(genome['GC1'], genome['GC2'], 'o', color = 'salmon')
plt.plot(of_interest['GC1'], of_interest['GC2'], 'o', color = 'steelblue')

plt.xlabel('GC: ' + genomename1, size = 15)
plt.ylabel('GC: ' + genomename2, size = 15)

plt.show()
plt.close()

