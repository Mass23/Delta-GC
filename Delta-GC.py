#!/usr/bin/env python3
# Author: Massimo Bourquin
# Contact: massimo.bourquin@unil.ch
# Arguments
# python3 Delta-GC.py M_P_alignment.xmfa ['genome1', 'genome2'] window-size threshold-of-uncertainty-to-drop-window [region of : interest]

class Scaffold:
    def __init__(self, beginning, end, scaffold):
        self.beginning = beginning
        self.end = end
        self.scaffold = scaffold

class AlignmentScaffold:
    def __init__(self, aln_scaffold_number, beginning, end, ref_sequence, alt_sequence):
        self.aln_scaffold_number = aln_scaffold_number
        self.beginning = beginning
        self.end = end
        self.ref_sequence = ref_sequence
        self.alt_sequence = alt_sequence

def AbsoluteValue(value):
    try:
        if value < 0:
            return(value * -1)
        else:
            return(value)
    except:
        return('nan')

def GetGenomeScaffoldsSize(reference_genome):
    from Bio import SeqIO

    scaffold_dict = dict()
    scaffold_list = list()

    with open(reference_genome) as genome:
        total_length = 0
        for record in SeqIO.parse(genome, 'fasta'):
            beginning = total_length
            length = len(record.seq)
            end = beginning + length

            total_length += len(record.seq)

            scaffold_dict[record.id] = Scaffold(beginning = beginning, end = end, scaffold = record.id)
            scaffold_list.append(record.id)

    return(scaffold_dict, scaffold_list)

def ParseMauve(alignment_file):
    from Bio import AlignIO

    # Separate pairs of sequences
    with open(alignment_file) as file:
        alignment = AlignIO.parse(file, 'mauve')
        aln_scaffold_count = 0

        aln_scaffold_dict = dict()
        aln_scaffold_list = list()

        for aln_scaffold in alignment:

            if len(aln_scaffold) == 2:
                aln_scaffold_count += 1

                ref_scaffold = aln_scaffold[0]
                alt_scaffold = aln_scaffold[1]

                ref_id = ref_scaffold.id.split('/')
                ref_beginning = int(ref_id[-1].split('-')[0])
                ref_end = int(ref_id[-1].split('-')[-1])

                ref_sequence = ref_scaffold.seq
                alt_sequence = alt_scaffold.seq

                aln_scaffold_list.append(aln_scaffold_count)
                aln_scaffold_dict[aln_scaffold_count] = AlignmentScaffold(aln_scaffold_count, ref_beginning, ref_end, ref_sequence, alt_sequence)

            else:
                break

        return(aln_scaffold_dict, aln_scaffold_list)

def GetIdentity(ref_sequence, alt_sequence):
    identity_count = 0
    seq_length = len(ref_sequence)

    for i in range(0, len(ref_sequence)):
        if ref_sequence[i] == alt_sequence[i]:
            identity_count += 1
        else:
            continue

    try:
        identity = identity_count / seq_length
    except:
        identity = 'nan'

    return(identity)

def GetGC(sequence):
    acount = sequence.count('A')
    ccount = sequence.count('C')
    gcount = sequence.count('G')
    tcount = sequence.count('T')

    try:
        gc_content = (gcount + ccount) / (acount + ccount + gcount + tcount)
    except:
        gc_content = 'nan'

    return(gc_content)

def GetUncertainty(sequence):
    ncount = sequence.count('-') + sequence.count('?') + sequence.count('N')
    length = len(sequence)

    try:
        uncertainty = ncount / length
    except:
        uncertainty = 'nan'
    return(uncertainty)

def sliding_window(sequence1, sequence2, window_size, length):

    windows_list = list()
    identity_list = list()
    uncertainty1_list = list()
    uncertainty2_list = list()
    gc1_list = list()
    gc2_list = list()
    delta_gc_list = list()
    window_count = 0

    # Sliding-window core part
    for window in range(0, length, window_size):
        sw_begin = window
        sw_end = window + window_size
        window1 = sequence1[sw_begin:sw_end]
        window2 = sequence2[sw_begin:sw_end]
        window_count += 1

        windows_list.append(window_count)
        identity_list.append(GetIdentity(window1, window2))
        uncertainty1_list.append(GetUncertainty(window1))
        uncertainty2_list.append(GetUncertainty(window2))
        gc1 = GetGC(window1)
        gc2 = GetGC(window2)
        gc1_list.append(gc1)
        gc2_list.append(gc2)
        try:
            delta_gc = gc1 - gc2
        except:
            delta_gc = 'nan'
        delta_gc_list.append(delta_gc)

    return(windows_list, identity_list, uncertainty1_list, uncertainty2_list, gc1_list, gc2_list, delta_gc_list)

def Main(reference_genome, alignment_file, window_size):

    print("Parsing reference genome...")
    parse_reference = GetGenomeScaffoldsSize(reference_genome)
    scaffold_dict = parse_reference[0]
    scaffold_list = parse_reference[1]
    print("Done!")

    print("Parsing alignment...")
    parse_mauve = ParseMauve(alignment_file)
    aln_scaffold_dict = parse_mauve[0]
    aln_scaffold_list = parse_mauve[1]
    print("Done!")

    print("Writing results...")
    # Create output file
    output_file = open('Delta_GC.csv', 'w')
    output_file.write(''.join(['Scaffold', '\t',
                            'Alignment_scaffold', '\t',
                            'Window_number','\t',
                            'Beginning', '\t',
                            'End', '\t',
                            'Delta_GC', '\t',
                            'Absolute_Delta_GC', '\t',
                            'GC1', '\t',
                            'GC2', '\t',
                            'Identity', '\t',
                            'Uncertainty', '\t',
                            'Uncertainty1', '\t',
                            'Uncertainty2', '\n']))

    # Iterate through each scaffold
    for scaffold in scaffold_list:

        scaffold_beginning = scaffold_dict[scaffold].beginning
        scaffold_end = scaffold_dict[scaffold].end


        for aln_scaffold in aln_scaffold_list:

            aln_scaffold_count = aln_scaffold_dict[aln_scaffold].aln_scaffold_number
            aln_scaffold_beginning = aln_scaffold_dict[aln_scaffold_count].beginning
            aln_scaffold_end = aln_scaffold_dict[aln_scaffold_count].end
            aln_scaffold_length = aln_scaffold_end - aln_scaffold_beginning

            aln_scaffold_ref_seq = aln_scaffold_dict[aln_scaffold_count].ref_sequence
            aln_scaffold_alt_seq = aln_scaffold_dict[aln_scaffold_count].alt_sequence


            # If the alignment scaffold is in the scaffold: process
            if aln_scaffold_beginning > scaffold_beginning and aln_scaffold_end < scaffold_end:

                # Process the alignment scaffold
                aln_scaffold_results = sliding_window(aln_scaffold_ref_seq, aln_scaffold_alt_seq, window_size, aln_scaffold_length)

                # Create lists to store the results
                windows_list = aln_scaffold_results[0]
                identity_list = aln_scaffold_results[1]
                uncertainty1_list = aln_scaffold_results[2]
                uncertainty2_list = aln_scaffold_results[3]
                gc1_list = aln_scaffold_results[4]
                gc2_list = aln_scaffold_results[5]
                delta_gc_list = aln_scaffold_results[6]

                # Iterate through the result lists
                for i in range(0, len(windows_list)):

                    # Add the results to the file
                    beginning = aln_scaffold_beginning + (i * window_size)
                    end = aln_scaffold_beginning + (i * window_size) + window_size
                    output_file.write(''.join([str(scaffold.replace('Scaffold', '')).replace('Contig', '10'), '\t',
                                            str(aln_scaffold_count), '\t',
                                            str(windows_list[i]),'\t',
                                            str(beginning), '\t',
                                            str(end), '\t',
                                            str(delta_gc_list[i]), '\t',
                                            str(AbsoluteValue(delta_gc_list[i])), '\t',
                                            str(gc1_list[i]), '\t',
                                            str(gc2_list[i]), '\t',
                                            str(identity_list[i]), '\t',
                                            str((uncertainty1_list[i] + uncertainty2_list[i]) / 2), '\t',
                                            str(uncertainty1_list[i]), '\t',
                                            str(uncertainty2_list[i]), '\n']))
                # Remove the aln_scaffold processed:
                aln_scaffold_list.remove(aln_scaffold)

            else:
                continue
    print("Done!")

import sys
args = sys.argv
reference_genome = args[1]
alignment_file = args[2]
window_size = int(args[3])

Main(reference_genome, alignment_file, window_size)
