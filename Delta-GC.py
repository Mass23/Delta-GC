# Author: Massimo Bourquin
# Contact: massimo.bourquin@unil.ch
# Arguments
# python3 Delta-GC.py"M_P_alignment.xmfa ['genom1', 'genome2'] window-size threshold-of-uncertainty-to-drop-window [region of : interest]

from scipy import stats
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
from Bio import SeqIO
import seaborn as sns
import statistics


# 1. Parse mauve .xmfa file
def absolute_value(value):
    try:
        if value < 0:
            return(value * -1)
        return(value)
    except:
        return(None)

def mauve_parser(alignment, genomes_list):
    ''' A function that takes as input a mauve .xmfa file and outputs two sorted scaffold class instances list, one per genome.'''

    genome1 = genomes_list[0]
    genome2 = genomes_list[1]

    genome1_seq = []
    genome2_seq = []
    
    print("Parsing alignment...")

    with open(alignment) as aln:
        data = AlignIO.parse(aln, 'mauve')
        current_genome = ''
    
        for record in data:

            for i in record:
                j = list(i.id.split("/"))
                
                if j[-2] == current_genome:
                    break

                elif j[-2] == genome1:
                    genome1_seq.extend([str(i.seq)])
                    current_genome = j[-2]
                    print(j[-2], j[-1])
                    continue
        
                elif j[-2] == genome2:
                    genome2_seq.extend([str(i.seq)])
                    current_genome = j[-2]
                    print(j[-2], j[-1])
                    continue
    
    genome1_seq = ''.join(genome1_seq)
    genome2_seq = ''.join(genome2_seq)
    print("Parsing done!")

    return([genome1, genome1_seq, genome2, genome2_seq])

# 2. Sliding-window function and results calculation
class result:
    def __init__(self, id_number, scaffold, diff_gc, gc1, gc2,identity, unknown, of_interest):
        self.id_number = int(id_number)
        self.diff_gc = diff_gc
        self.gc1 = gc1
        self.gc2 = gc2
        self.identity = identity
        self.unknown = unknown
        self.scaffold = scaffold
        self.of_interest = of_interest
        self.abs_dgc = absolute_value(diff_gc)

def sliding_window(genome1, genome2, size, scaffold_coord, scaffold_of_interest):
    ''' A function that takes two concatenated genomes and compares sliding-windows.'''

    # Create the output file
    windows_out = list()

    # Only positive numbers are allowed
    if size < 0:
        print("Negative value not allowed for the sliding-window size!")

    else:
        count_windows = 0

        # Line implemented to process on the length of the smallest sequence if size differs
        small_seq = min(len(genome1), len(genome2))

        # Sliding-window core part
        
        for i in range(0, small_seq, size):
            sw_begin = int(i)
            sw_end = int(i + size)
            window1 = str(genome1[sw_begin:sw_end])
            window2 = str(genome2[sw_begin:sw_end])
            count_windows += 1

            scaffold = 0
            for k in scaffold_coord:
                if k > sw_begin:
                    scaffold = int(scaffold_coord.index(k)) + 1
                    break
                else:
                    continue

            # counters initialisation
            acount1 = 0
            ccount1 = 0
            gcount1 = 0
            tcount1 = 0
            ncount1 = 0

            acount2 = 0
            ccount2 = 0
            gcount2 = 0
            tcount2 = 0
            ncount2 = 0

            identity_count = 0

            for nuc in range(min([len(window1), len(window2)])):
                # identity
                if window1[nuc] == window2[nuc]:
                    if window1[nuc] != 'N' or '?' or '-':
                        identity_count += 1
                
                # counters seq 1
                if window1[nuc] == 'A':
                    acount1 += 1
                elif window1[nuc] == 'C':
                    ccount1 += 1
                elif window1[nuc] == 'G':
                    gcount1 += 1
                elif window1[nuc] == 'T':
                    tcount1 +=1
                elif window1[nuc] == 'N' or '?' or '-':
                    ncount1 += 1
                else:
                    print("Char not allowed!")
                
                # counters seq 2
                if window2[nuc] == 'A':
                    acount2 += 1
                elif window2[nuc] == 'C':
                    ccount2 += 1
                elif window2[nuc] == 'G':
                    gcount2 += 1
                elif window2[nuc] == 'T':
                    tcount2 +=1
                elif window2[nuc] == 'N' or '?' or '-':
                    ncount2 += 1
                else:
                    print("Char not allowed!")
            
            unknown = float((ncount1 + ncount2) / (2 * size))

            nuc_count1 = acount1 + ccount1 + gcount1 + tcount1
            nuc_count2 = acount2 + ccount2 + gcount2 + tcount2
            
            try:
                gc1 = (ccount1 + gcount1)/nuc_count1
            except:
                gc1 = None
            
            try:
                gc2 = (ccount2 + gcount2)/nuc_count2
            except:
                gc2 = None
            
            try:
                diff_gc = gc1 - gc2
            except:
                diff_gc = None

            try:
                identity = identity_count / min([nuc_count1, nuc_count2])
            except:
                identity = 0

            if scaffold == scaffold_of_interest:
                results_i = result(count_windows, scaffold, diff_gc, gc1, gc2, identity, unknown, 1)
            else:
                results_i = result(count_windows, scaffold, diff_gc, gc1, gc2, identity, unknown, 0)

            windows_out.append(results_i)

    return windows_out

def main(alignment, genomes_list, window_size, threshold, interest_scaffold):

    # Parsing
    parsed_alignment = mauve_parser(alignment, genomes_list)
    genomename1 = parsed_alignment[0]
    genomename2 = parsed_alignment[2]

    scaffold_of_interest = interest_scaffold

    genome_ref = SeqIO.parse(open(genomename1), 'fasta')
    scaffold_coord = list()
    scaffold_coord.extend([0])

    for i in genome_ref:
        before_length = 0
        for k in range(0,len(scaffold_coord)):
            before_length += scaffold_coord[k]
        length = len(i.seq)
        actual_length = int(before_length + length)
        scaffold_coord.append(actual_length)

    # Genomes sequences
    genome1_seq = parsed_alignment[1]
    genome2_seq = parsed_alignment[3]
    
    print("Concatenation done!")
    print(genomename1, " size: ", len(genome1_seq))
    print(genomename2, " size: ", len(genome2_seq))

    # Sliding-window core function
    results = sliding_window(genome1_seq, genome2_seq, window_size, scaffold_coord, scaffold_of_interest)

    # Plot and file output
    print("Writing results...")

    out_file = open("delta_gc.csv", 'w')
    out_file.writelines("window_number\tscaffold\tdelta_gc\tgc1\tgc2\tidentity\tuncertainty\tof_interest\tabs_dgc\n")

    for i in results:
        if i.unknown > threshold:
            i.identity = np.nan
            i.gc2 = np.nan
            i.gc1 = np.nan
            i.diff_gc = np.nan
            i.uncertainty = np.nan
            i.abs_dgc = np.nan
            out_file.writelines(str(i.id_number) + '\t' + str(i.scaffold) + '\t' + str(i.diff_gc) + '\t' + str(i.gc1) + '\t' + str(i.gc2) + '\t' + str(i.identity) + '\t' + str(i.unknown) + '\t' + str(i.of_interest) + '\t' + str(i.abs_dgc) + "\n")

        else:
            out_file.writelines(str(i.id_number) + '\t' + str(i.scaffold) + '\t' + str(i.diff_gc) + '\t' + str(i.gc1) + '\t' + str(i.gc2) + '\t' + str(i.identity) + '\t' + str(i.unknown) + '\t' + str(i.of_interest) + '\t' + str(i.abs_dgc) + "\n")
    
    out_file.close() 


# Plotting:
    data = pd.read_csv('delta_gc.csv', delimiter="\t")

    mask_of_interest = data['of_interest'].eq(1)
    mask_genome = data['of_interest'].eq(0)

    of_interest = data[mask_of_interest].dropna()
    genome = data[mask_genome].dropna()

# Plot delta-GC
    plt.figure(1)
    plt.title("Delta-GC")
    sns.set_style('white')
    sns.set_context("paper")

    # Horizontal line at 0
    plt.plot((0, data['abs_dgc'].count()), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 0.8)

    # Plot the dgc    
    plt.plot(genome['window_number'], genome['delta_gc'], 'o', color = 'salmon')
    plt.plot(of_interest['window_number'], of_interest['delta_gc'], 'o', color = 'steelblue')

    # Plot absolute value of dgc
    plt.plot(data['window_number'],data['abs_dgc'], '-', color = 'lightgrey')

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
    sns.set_style('white')
    sns.set_context("paper")

    sns.kdeplot(genome['gc1'], genome['gc2'], cmap="Reds", shade=True, shade_lowest=False)
    sns.kdeplot(of_interest['gc1'], of_interest['gc2'], cmap="Blues", shade=True, shade_lowest=False)

    plt.xlabel(genomename1 + " GC content", size = 15)
    plt.ylabel(genomename2 + " GC content", size = 15)

    plt.show()
    plt.close()

# Plot distributions
    plt.figure(5)
    plt.title("Density plot")
    sns.set_style('white')
    sns.set_context("paper")

    sns.set_style('whitegrid')
    sns.distplot(genome['delta_gc'], color="salmon", hist=False)
    sns.distplot(of_interest['delta_gc'], color="steelblue", hist=False)

    plt.show()
    plt.close()

# F-test
    F_test = statistics.variance(of_interest['delta_gc']) / statistics.variance(genome['delta_gc'])
    df1 = of_interest['delta_gc'].count() - 1
    df2 = genome['delta_gc'].count() - 1

    p_value = stats.f.cdf(F_test, df1, df2)
    print("F-test: p.value = ", p_value)

main("M_P_alignment.xmfa", ['fsel_M.fasta', 'fsel_P.fasta'], 50000, 0.3, 3)
