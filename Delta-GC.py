# Author: Massimo Bourquin
# Contact: massimo.bourquin@unil.ch
# Arguments
# python3 Delta-GC.py"M_P_alignment.xmfa ['genom1', 'genome2'] window-size threshold-of-uncertainty-to-drop-window [region of : interest]

import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
from Bio import SeqIO

# 1. Parse mauve .xmfa file
def mauve_parser(alignment, genomes_list):
    ''' A function that takes as input a mauve .xmfa file and outputs two sorted scaffold class instances list, one per genome.'''

    genome1 = genomes_list[0]
    genome2 = genomes_list[1]

    genome1_seq = []
    genome2_seq = []
    
    print("Parsing alignment...")

    with open(alignment) as aln:
        data = AlignIO.parse(aln, 'mauve')
    
        for record in data:
            current_genome = ''

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
    def __init__(self, id_number, scaffold, diff_gc, gc1, gc2,identity, unknown):
        self.id_number = int(id_number)
        self.diff_gc = diff_gc
        self.gc1 = gc1
        self.gc2 = gc2
        self.identity = identity
        self.unknown = unknown
        self.scaffold = scaffold

def sliding_window(genome1, genome2, size, scaffold_coord):
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

            results_i = result(count_windows, scaffold, diff_gc, gc1, gc2, identity, unknown)

            windows_out.append(results_i)

    return windows_out

def main(alignment, genomes_list, window_size, threshold, interest_scaffold):
    
    # Parsing
    parsed_alignment = mauve_parser(alignment, genomes_list)
    genome1 = parsed_alignment[0]
    genome2 = parsed_alignment[2]

    scaffold_of_interest = interest_scaffold

    genome_ref = SeqIO.parse(open(genome1), 'fasta')
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
    print(parsed_alignment[0], " size: ", len(genome1_seq))
    print(parsed_alignment[2], " size: ", len(genome2_seq))

    # Sliding-window core function
    results = sliding_window(genome1_seq, genome2_seq, window_size, scaffold_coord)

    # Plot and file output
    print("Writing results...")

    out_file = open("delta_gc.csv", 'w')
    #out_file.writelines("window_number\tscaffold\tdelta_gc\tgc1\tgc2\tidentity\tuncertainty\n")

    window_number = []
    scaffold = []
    delta_gc = []
    gc1 = []
    gc2 = []
    identity = []
    uncertainty = []

    for i in results:
        window_number.append(i.id_number)
        scaffold.append(i.scaffold)

        if i.unknown > threshold or i.identity < 0.35:
            i.identity = np.nan
            identity.append(i.identity)
            i.gc2 = np.nan
            gc2.append(i.gc2)
            i.gc1 = np.nan
            gc1.append(i.gc1)
            i.diff_gc = np.nan
            delta_gc.append(i.diff_gc)
            i.uncertainty = np.nan
            uncertainty.append(i.unknown)
            out_file.writelines(str(i.id_number) + '\t' + str(i.scaffold) + '\t' + str(i.diff_gc) + '\t' + str(i.gc1) + '\t' + str(i.gc2) + '\t' + str(i.identity) + '\t' + str(i.unknown) + "\n")

        else:
            identity.append(i.identity)
            gc2.append(i.gc2)
            gc1.append(i.gc1)
            delta_gc.append(i.diff_gc)
            uncertainty.append(i.unknown)
            out_file.writelines(str(i.id_number) + '\t' + str(i.scaffold) + '\t' + str(i.diff_gc) + '\t' + str(i.gc1) + '\t' + str(i.gc2) + '\t' + str(i.identity) + '\t' + str(i.unknown) + "\n")
    
    out_file.close() 

# Output

    in_scaffold_of_interest = np.empty(len(window_number))
    for i in range(0, len(window_number)):
        if scaffold[i] == scaffold_of_interest:
            in_scaffold_of_interest[i] += 1
    
    abs_gc = []
    for i in delta_gc:
        try:
            i = abs(i)
            abs_gc.append(i)
        except:
            i = i
            abs_gc.append(i)

# Plot delta-GC
    #Plot an horizontal line at 0, -1 and 1
    plt.figure(1)

    #Plot the results
    plt.plot(window_number,abs_gc, '-', color = 'lightgrey')
    
    for i in range(0, len(window_number)):
        if scaffold[i] != scaffold_of_interest:
            plt.plot(window_number[i], delta_gc[i],  'o', color = 'salmon')
        else:
            plt.plot(window_number[i], delta_gc[i], 'o', color = 'steelblue')
    
    plt.plot((1, len(delta_gc)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 0.8)
    plt.plot((1, len(delta_gc)), (-1, -1), linestyle = 'solid', color = 'dimgrey', lw = 1)
    plt.plot((1, len(delta_gc)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    # Axis labels
    plt.axis([1, len(delta_gc), -0.2, 0.2])
    plt.ylabel('ΔGC (' + genome1 + ' - ' + genome2 + ')', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)

    plt.show()
    plt.close()


# Plot Identity
    plt.figure(2)
    #Plot an horizontal line at 0, -1 and 1
    plt.plot((1, len(identity)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 1)
    plt.plot((1, len(identity)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(window_number, identity, 'o', color= 'red')
    plt.axis([1, len(identity), 0, 1.1])

    # Axis labels
    plt.ylabel('Identity', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)
    plt.show()
    plt.close()


# Plot Uncertainty
    plt.figure(3)
    #Plot an horizontal line at 0, -1 and 1
    plt.plot((1, len(uncertainty)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 1)
    plt.plot((1, len(uncertainty)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(window_number, uncertainty, 'o', color= 'green')
    plt.axis([1, len(uncertainty), 0, 1.1])
    
    # Axis labels
    plt.ylabel('Uncertainty', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)
    
    # Output
    plt.show()
    plt.close()

    print("Output Manhattan plots Done!")

# Plot GC1/GC2
    plt.figure(4)

    for i in range(0, len(window_number)):
        if scaffold[i] == scaffold_of_interest:
            plt.scatter(gc1[i], gc2[i], color = 'steelblue')
        else:
            plt.scatter(gc1[i],gc2[i],color = 'salmon')
    plt.xlabel(genome1, size = 15)
    plt.ylabel(genome2, size = 15)

    plt.show()
    plt.close()

main("M_P_alignment.xmfa", ['fsel_M.fasta', 'fsel_P.fasta'], 20000, 0.2, 3)
