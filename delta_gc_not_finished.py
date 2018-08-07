""" The aim of this program is to take as an input two aligned genomes in a mauve .xmfa file and assesss a differential GC content analysis using a sliding-window approach. """

##########################################################################################################################################################################
##### 1. Classes creation and Mauve .xmfa file parsing
class results:
    def __init__(self, window_n, gc1, gc2, dgc, window_begin, window_end):
        self.window_n = window_n
        self.window_begin = window_begin
        self.window_end = window_end
        self.gc1 = gc1
        self.gc2 = gc2
        self.dgc = dgc

class scaffold:
    def __init__(self, genome_name, scaffold, sequence, results):
        self.genome_name = genome_name
        self.scaffold = scaffold
        self.sequence = sequence
        self.begin = scaffold.split("-")[0]
        self.end = scaffold.split("-")[1]
        self.results = results

    def get_sequence(self, genome_name, scaffold):
        return self.sequence

def parse_mauve(in_file):
    import operator
    from Bio import AlignIO
    
    genomes_list = list()
    
    data = AlignIO.parse(open(in_file), "mauve")
    parsed_data = list()
    
    for record in data:
        for i in record:
            j = list(i.id.split("/"))
            scaffold_i = scaffold(j[-2], j[-1], i.seq)
            print(j[-2], j[-1], i.seq)
            
            if j[-2] in genomes_list:
                continue
            else:
                genomes_list.append(j[-2])
            parsed_data.append(scaffold_i)



##########################################################################################################################################################################
##### 2. GC content function

def gc_content(seq):
    """ Simple function to calculate the GC content given a sequence.
    If the number of unknown nucleotides exceeds a given threshold, the function returns -1"""

    # Counters
    acount = seq.count("A")
    ccount = seq.count("C")
    gcount = seq.count("G")
    tcount = seq.count("T")
    
    ncount = seq.count("N" or "-" or "—" or "?")
    
    nuc_count = acount + ccount + gcount + tcount

    # Count the number of unknown nucleotides in order to discard too uncertain regions
    try:
        n_cont = ncount / (nuc_count + ncount)
    
    # If division by zero, fails --> n_cont == 1
    except:
        n_cont = 1
    threshold_undefined = 0.1

    # Try statement to test whether the window percentage of unknown nucleotides is higher or lower than the threshold defined
    try:
        if n_cont < threshold_undefined:
		    # GC content - placed here to avoid division by zero if the sequence is fully undetermined
            gc_cont = (gcount + ccount) / nuc_count
            return(float(gc_cont))
    # If too many unknown nucleotides, the program returns -1
    except:
        return(int(-1))


##########################################################################################################################################################################
##### 3. Sliding-window function and LCB iteration

def sliding_window(seq1, seq2, size):
    """ Function that returns the sliding window number, the coordinates in the sequence, the difference in GC content (Also plotted) when two sequences are given. 
    The window size can be adjusted."""

    import matplotlib.pyplot as plt

    # Create the output file
    windows_out = list()

    # Only positive numbers are allowed
    if size < 0:
        print("Negative value not allowed for the sliding-window size!")

    else:
        count_windows = 0

        # Line implemented to process on the length of the smallest sequence if size differs
        small_seq = min(len(seq1), len(seq2))

        # Sliding-window core part
        for i in range(0, small_seq, size):
            sw_begin = int(i)
            sw_end = int(i + size)
            window1 = seq1[sw_begin:sw_end]
            window2 = seq2[sw_begin:sw_end]
            count_windows += 1
            
            GC1 = gc_content(window1)
            GC2 = gc_content(window2)

            # If too many unknown nucleotides, GC function returns -1, if one of the two GC results is equal to -1, the window result is = None
            try:
                diff_GC = gc1 - gc2
            except:
                diff_GC = None

            result_i = results(count_windows, GC1, GC2, diff_GC, sw_begin + 1, sw_end)
            gc_out.append(windows_out)

    return windows_out


##########################################################################################################################################################################
##### 4. Output Manhattan plot and main function

def Main():
    parse_mauve("M_P_Alignment.xmfa")
    
    Genome_1 = genomes_list[0]
    Genome_2 = genomes_list[1]

    #
    gc_out = open("gc_out.txt")
    for scaffold in scaffold_list:
        for sequence in scaffold.genome_name:

        sequence_1 = sequence.genome_name.sequence
        sequence_2 = scaffold.genome_name.sequence
        result_scaffold = sliding_window(sequence_1, sequence_2, 5000)
        scaffold(genome_name, scaffold, )
        
        gc_out.writelines("Scaffold: ", scaffold, "(",scaffold.begin, "-", scaffold.end, ")")
        for window in scaffold.results:
            gc_out.writelines("Window ", window.window_n, " (Pos.",window.window_begin , " - Pos.", window.window_end ,"): \t", window.dgc ," (", ,": ", gc1, " - ", ,": ",gc2,")" )
    
    #Plot an horizontal line at 0, -1 and 1
    plt.plot((1, len(x)), (0, 0), linestyle = 'dashed', color = 'darkgrey', lw = 0.8)
    plt.plot((1, len(x)), (-1, -1), linestyle = 'solid', color = 'dimgrey', lw = 1)
    plt.plot((1, len(x)), (1, 1), linestyle = 'solid', color = 'dimgrey', lw = 1)

    #Plot the results
    plt.plot(x, y, 'o', color= 'royalblue')
    plt.axis([1, len(x), -1.1, 1.1])

    # Axis labels
    plt.ylabel('ΔGC', size = 15)
    plt.xlabel('Genome windows', size = 15)

    # Output
    plt.show()
    plt.close()

    print("Output Manhattan plot Done!")

