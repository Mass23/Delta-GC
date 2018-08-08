class scaffold:
    def __init__(self, genome_name, scaffold, sequence):
        self.genome_name = genome_name
        self.scaffold = scaffold
        self.sequence = sequence
        self.begin = scaffold.split("-")[0]
        self.end = scaffold.split("-")[1]

# 1. Parse mauve .xmfa file
def mauve_parser(alignment, genomes_list):
    ''' A function that takes as input a mauve .xmfa file and outputs two sorted scaffold class instances list, one per genome.'''
    from Bio import AlignIO
    import operator

    genome1 = genomes_list[0]
    genome2 = genomes_list[1]

    genome1_scaffolds = list()
    genome2_scaffolds = list()
    
    print("Parsing alignment...")

    data = AlignIO.parse(open(alignment), "mauve")
    
    for record in data:
        for i in record:
            j = list(i.id.split("/"))
            if j[-2] == genome1:
                scaffold_i = scaffold(genome1, j[-1], i.seq)
                genome1_scaffolds.append(scaffold_i)
            
            elif j[-2] == genome2:
                scaffold_i = scaffold(genome1, j[-1], i.seq)
                genome2_scaffolds.append(scaffold_i)
            
            print(j[-2], j[-1])
    print("Parsing done!")

    return([genome1, genome1_scaffolds, genome2, genome2_scaffolds])

# 2. Concatenate genome
def genome_concatenator(genome):
    ''' A function that takes as input a sorted scaffold class instance list and concatenates the sequence into a genome. Outputs the genome sequence.'''
    genome_seq = str()
    scaffold_list = list()

    for i in genome:
        genome_seq += i.sequence
        scaffold_list.extend(i.begin)

    return(genome_seq, scaffold_list)

# 3. Sliding-window function and results calculation
class result:
    def __init__(self, id_number, diff_gc, gc1, gc2,identity, unknown):
        self.id_number = int(id_number)
        self.diff_gc = diff_gc
        self.gc1 = gc1
        self.gc2 = gc2
        self.identity = identity
        self.unknown = unknown

def sliding_window(genome1, genome2, size):
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

            results_i = result(count_windows, diff_gc, gc1, gc2, identity, unknown)

            windows_out.append(results_i)

    return windows_out

def main(alignment, genomes_list, window_size, threshold):
    import matplotlib.pyplot as plt

    # Parsing
    parsed_alignment = mauve_parser(alignment, genomes_list)
    genome1 = parsed_alignment[0]
    genome2 = parsed_alignment[2]

    genome1_scaffolds = parsed_alignment[1]
    genome2_scaffolds = parsed_alignment[3]

    # Genomes concatenation

    genome1_seq = genome_concatenator(genome1_scaffolds)
    #scaffolds_separation_1 = genome1_seq[1] to use maybe for the output plot?
    genome1_seq = genome1_seq[0]
    
    genome2_seq = genome_concatenator(genome2_scaffolds)
    #scaffolds_separation_2 = genome2_seq[1]
    genome2_seq = genome2_seq[0]
    
    print("Concatenation done!")
    print(genome1, " size: ", len(genome1_seq))
    print(genome2, " size: ", len(genome2_seq))

    # Sliding-window core function
    results = sliding_window(genome1_seq, genome2_seq, window_size)

    # Plot and file output
    print("Writing results...")

    out_file = open("delta_gc.txt", 'w')
    out_file.writelines("Window_number \t delta_gc \t gc1 \t gc2 \t identity \t uncertainty \n")

    identity_list = []
    deltagc_list = []
    number_list = []
    uncertain_list = []

    for i in results:
        if i.unknown > threshold:
            i.identity = None
            i.gc2 = None
            i.gc1 = None
            i.diff_gc = None
            out_file.writelines(str(i.id_number) + '\t' + str(i.diff_gc) + '\t' + str(i.gc1) + '\t' + str(i.gc2) + '\t' + str(i.identity) + '\t' + str(i.unknown) + "\n")
            identity_list.append(i.identity)
            deltagc_list.append(i.diff_gc)
            number_list.append(i.id_number)
            uncertain_list.append(i.unknown)
        else:
            out_file.writelines(str(i.id_number) + '\t' + str(i.diff_gc) + '\t' + str(i.gc1) + '\t' + str(i.gc2) + '\t' + str(i.identity) + '\t' + str(i.unknown) + "\n")
            identity_list.append(i.identity)
            deltagc_list.append(i.diff_gc)
            number_list.append(i.id_number)
            uncertain_list.append(i.unknown)
    
    out_file.close()

# Plot delta-GC
    #Plot an horizontal line at 0, -1 and 1
    plt.figure(1)

    #Plot the results
    plt.plot(number_list, deltagc_list,  'o', color= 'royalblue')
    
    plt.plot((1, len(deltagc_list)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 0.8)
    plt.plot((1, len(deltagc_list)), (-1, -1), linestyle = 'solid', color = 'dimgrey', lw = 1)
    plt.plot((1, len(deltagc_list)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    # Axis labels
    plt.axis([1, len(deltagc_list), -1.1, 1.1])
    plt.ylabel('ΔGC (' + genome1 + ' - ' + genome2 + ')', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)

    plt.show()
    plt.close()


# Plot Identity
    plt.figure(2)
    #Plot an horizontal line at 0, -1 and 1
    plt.plot((1, len(identity_list)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 1)
    plt.plot((1, len(identity_list)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(number_list, identity_list, 'o', color= 'red')
    plt.axis([1, len(identity_list), 0, 1.1])

    # Axis labels
    plt.ylabel('Identity', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)
    plt.show()
    plt.close()


# Plot Uncertainty
    plt.figure(3)
    #Plot an horizontal line at 0, -1 and 1
    plt.plot((1, len(uncertain_list)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 1)
    plt.plot((1, len(uncertain_list)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(number_list, uncertain_list, 'o', color= 'green')
    plt.axis([1, len(uncertain_list), 0, 1.1])
    
    # Axis labels
    plt.ylabel('Uncertainty', size = 15)
    plt.xlabel('Genome windows (size:' + str(window_size) + 'bp.)', size = 15)
    
    # Output
    plt.show()
    plt.close()

    print("Output Manhattan plots Done!")


main("M_P_alignment.xmfa", ['fsel_M.fasta', 'fsel_P.fasta'], 5000, 0.1)
