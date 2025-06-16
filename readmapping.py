# -*- coding: utf-8 -*-
"""
INDU- Read mapping

@author: Javier Hidalgo
"""
import matplotlib.pyplot as plt
import argparse
import os

def check_dna(seq):
    """
    

    Checks if a DNA has invalid characters.
    
    Parameters
    ----------
    seq : string
        Sequence of DNA.

    Raises
    ------
    ValueError
        The sequence contains a value different from A,T,G,C.

    Returns
    -------
    None.

    """
    #Checks if there are characters in a DNA sequence that are not ATGC
    if len(set(seq) - set({'T', 'A', 'G', 'C'})) > 0:
        raise ValueError("DNA sequence {0} contains at least one invalid character".format(seq))

def read_fasta(filename):
    """
    Reads a fasta file and store the information in a dictionary.

    Parameters
    ----------
    filename : str
        Name of the fasta file (must include .fa or .fasta).

    Returns
    -------
    out : dict
        A dictionary where the keys are the sequence annotations and the values
        are the sequences of DNA.
    """
    out={}   
    try:
        h= open(filename, "r")
        if filename[-3:] != '.fa' and filename[-5:] != '.fasta':
            raise ValueError('The file must be .fa or .fasta')
        for line in h:
            line= line.strip() # We take away the \n to use check_dna() afterwards
            if '>' in line:
                key= line[1:] #Keys are introduced by >
                out[key]=0 # Create an element in the dict with key=line and no value
            elif out[key]==0: #If there is no value, we add the line
                check_dna(line)
                out[key]= line
            else: # If there was already a line, we add the following (sequences that are >1 line long)
                check_dna(line)
                out[key]+=line 
    except FileNotFoundError:
        print(filename, "was not found on the working directory")
        out= None
        
    return out


             
def get_kmers(seq, k):
    """
    Given a sequence, it returns a dictionary with all different k-mers and
    the indexes in which they are located in the sequence

    Parameters
    ----------
    seq : str
        A sequence of DNA.
    k : int
        Length of the k-mers.

    Returns
    -------
    kmers : dict
        A dictionary in which the keys are the sequences of the k-mers and the
        values are a list with the indexes in which these k-mers are located.

    """
    kmers = {}
    if k > len(seq):
        raise ValueError('The length of the k-mers cannot exceed the length of the sequence')
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k] #Get the kmers
        if kmer not in kmers:
            kmers[kmer] = []
        kmers[kmer].append(i) #Append the index
    return kmers


def kmers_alignment(seq_read, dict_ref, k):
    """
    It breaks the read and the reference sequences into their kmers and store the
    kmers from the read that align with the kmers in the reference and their position
    in which the homology starts in the reference.

    Parameters
    ----------
    seq_read : str
        A sequence of DNA of the read.
    dict_ref : dict
        Dictionary of the reference sequences. It can be obtained using read_fasta()
    k : int
        Length of the k-mers.

    Returns
    -------
    output : dict
        A dictionary where the keys are the reference annonations of the reference sequences
        that align with the read and the values are lists of the indexes where the homology
        starts.

    """
    reads_kmers = get_kmers(seq_read, k)  # Get the kmers from the read
    output = {}

    for key, ref_seq in dict_ref.items():  
        ref_kmers = get_kmers(ref_seq, k)  # Get the kmers from the ref
        # Look for coincidences between kmers in the read and the ref
        output[key] = []
        for kmer, indexes in reads_kmers.items():
            for j in indexes:
                if kmer in ref_kmers:
                    # Get all the positions for that kmer
                    for pos in ref_kmers[kmer]:
                        adjusted_pos = pos - j # Adjust the position of the k-mer by the index of the read k-mer
                 # Positions must be unique and in between (0, len(ref_seq)-len(seq_read)) to calculate the h_distance
                        if adjusted_pos not in output[key] and 0 <= adjusted_pos <= (len(ref_seq)-len(seq_read)):
                            output[key].append(adjusted_pos) 
    return output



def hamming_distance(seq1, seq2):
    """
    It calculates the number of mismatches between two sequences

    Parameters
    ----------
    seq1 : str
        A DNA sequence.
    seq2 : str
        A DNA sequence.

    Returns
    -------
    distance : int
        Hamming distance: number of mismatches between two sequences.

    """
    if len(seq1) != len(seq2):
        raise ValueError('Both sequences must have the same length')
    distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
    return distance


def best_mapping_location(seq_read, dict_ref, k):
    """
    It takes a sequence and a dictionary with reference sequence and returns
    a dictionary with the aligments with lowest hamming distance.

    Parameters
    ----------
    seq_read : str
        A sequence of DNA of the read.
    dict_ref : dict
        Dictionary of the reference sequences. It can be obtained using read_fasta()
    k : int
        Length of the k-mers.

    Returns
    -------
    best_mapping : dict
        A dictionary with the best mapping locations. Keys are the reference annotations
        and the values are lists that contain tuples with the index of the mapping location
        and the hamming distance.

    """
    if k > len(seq_read):
        raise ValueError("The length of the k-mer cannot exceed the length of the sequence")
    kmers_indexes = kmers_alignment(seq_read, dict_ref, k)
    best_mapping = {}
    lowest_hamming_d = len(seq_read)

    for ref, index in kmers_indexes.items():
        for i in index:
            hamming_d = hamming_distance(seq_read, dict_ref[ref][i: i+len(seq_read)])
            if hamming_d < lowest_hamming_d:
                # If we find a lower hamming distance, we restore the dict to empty
                best_mapping = {}
                # And store the new positions
                best_mapping[ref] = [(i, hamming_d)]
                lowest_hamming_d = hamming_d
            elif hamming_d == lowest_hamming_d:
                if ref not in best_mapping:
                    best_mapping[ref] = []
                best_mapping[ref].append((i, hamming_d))
    # best_mapping (for read_i) = {'ref_i': []}
    return best_mapping

def get_data_for_output(reads, refs, k):
    """
    It takes the dictionary from reads and references sequences and returns the
    data needed to obtain the output files

    Parameters
    ----------
    reads : dict
        A dictionary where keys are the anotations and values are the sequences.
    refs : dict
        A dictionary where keys are the anotations and values are the sequences.
    k : int, optional
        Length of the kmers.

    Returns
    -------
    file1_dict : dict
        A dictionary where keys are the read anotations and the values are a dictionary
        with reference anotations as keys and a list with tuples (index of best match location, errors).
    file2_dict : dict
        A dictionary with reference anotations as keys and a list with tuples (read anotation, length of read
        sequence, length of reference sequence) of the best mapping locations.
    plot_dict : dict
        A dictionary where keys are the read anotations and the values are the number of best locations
        for each read (defined as m_i in coverage calculation).

    """
    file1_dict = {}
    file2_dict = {key: [] for key in refs}
    plot_dict = {}
    
    for read_annotation in reads:
        # We get the best mapping location for each read
        mapping = best_mapping_location(reads[read_annotation], refs, k)
        # File for output one
        file1_dict[read_annotation] = mapping
        # For each read we store the number of best aligments (needed for plots)
        counts = sum(len(mapping[ref]) for ref in mapping)
        plot_dict[read_annotation] = counts
        
        for ref_annotation in mapping:
            for i in range(len(mapping[ref_annotation])):
                # Store information for output 2
                file2_dict[ref_annotation].append((read_annotation, 
                                                   len(reads[read_annotation]), 
                                                   len(refs[ref_annotation])))
                
    return file1_dict, file2_dict, plot_dict
        
def write_output_file_1(file1_dict, output_filename):
    """
    It takes the output one of get_data_for_output() and write a .txt file with
    'read anotation', 'reference anotation', 'index', 'errors'.

    Parameters
    ----------
    file1_dict : dict
        A dictionary where keys are the read anotations and the values are a dictionary
        with reference anotations as keys and a list with tuples (index of best match location, errors).
    output_filename : str
        Name of the output file.

    Returns
    -------
    None.

    """
    with open(output_filename + '_1.txt', 'w') as h_out_1:
        for read_annotation in file1_dict:
            mapping = file1_dict[read_annotation]
            for ref_annotation in mapping:
                for i in range(len(mapping[ref_annotation])):
                    # Write the file for output 1
                    h_out_1.write(
                        f"{read_annotation},"
                        f"{ref_annotation},"
                        f"{mapping[ref_annotation][i][0]},"
                        f"{mapping[ref_annotation][i][1]}\n"
                    )

def calculate_coverage( file2_dict, m_i_counts, output_filename):
    """
    It take the data from the second dictionary returned by get_data_for_output(),
    and the m_i_counts from the third dictionary returned by the same function,
    calculates the coverage (Cj = sum(r_i(m_i_j/m_i)/l_j)) and write a .txt file
    with the reference anotation and its coverage

    Parameters
    ----------
    file2_dict : dict
        A dictionary with reference anotations as keys and a list with tuples (read anotation, length of read
        sequence, length of reference sequence) of the best mapping locations.
    m_i_counts : dict
        A dictionary where keys are the read anotations and the values are the number of best locations
        for each read (defined as m_i in coverage calculation).
    output_filename : str
        Name of the output file.

    Returns
    -------
    None.

    """
    with open(output_filename + '_2.txt', 'w') as h_out_2:
        for ref_annotation in file2_dict:
            coverage = 0
            # Tuples are for each reference annotation are:
                # (read_annotation, read_length, ref_length)
            coverage_tuples = file2_dict[ref_annotation]
            # Calculate coverage
            for tup in coverage_tuples:
                read_name, read_length, ref_length = tup
                # m_ij is hardcoded to 1. If m_ij was > 1, it will already be
                # accounted with two tuples with the same read_name
                m_ij = 1
                m_i = m_i_counts[read_name]
                try: 
                    coverage += (read_length * (m_ij / m_i))/ref_length
            
                except:
                    # We can have errors if ref_length is not defined
                    # That is the case of when a reference sequence does not align with any read
                    # Therefore, we set coverage to 0 manually
                    coverage = 0
            h_out_2.write(f"{ref_annotation},{coverage}\n")
            
def get_plot(plot_dict, output_filename):
    """
    Takes the third dictionary from get_data_for_output() and makes an histogram of
    the frequency of the number of best mapping location that reads have.

    Parameters
    ----------
    plot_dict : dict
        A dictionary where keys are the read anotations and the values are the number of best locations
        for each read (defined as m_i in coverage calculation).
    output_filename : str
        Name of the output file.

    Returns
    -------
    None.

    """
    # We get a list with all the values to build the histogram
    histogram_data = list(plot_dict.values())
    # We get the max value to adjust the x-axis
    max_val = max(histogram_data)
    plt.hist(histogram_data, bins= range(0, max_val+2), align= 'left')
    plt.xlabel('Nr. best locations')
    plt.ylabel('Counts')
    plt.savefig(output_filename + '.pdf')
    
def main():
    parser = argparse.ArgumentParser('The program takes a fasta file with reads and a fasta file with reference\
                                     sequences and produces three files: 1. All best locations to a file\
                                     2. The reference coverage for each reference 3. A histogram plot displaying the\
                                     number of mapping locations per read')
    parser.add_argument('--k', type = int, default= 10,
                        help = "Length of the k-mers. It should be chosen by the user. The default is 10\
                        because it the maximum k that can detect 4 errors in a read that is 50 pb long (Pigeon principle)")
    parser.add_argument('ref_filename', type = str,
                        help = 'Filename of the fasta with reference sequences')
    parser.add_argument('read_filename', type = str,
                        help = 'Filename of the fasta with read sequences')
    parser.add_argument('output_folder', type = str,
                        help = 'Folder where the output files go')
    
    args = parser.parse_args()
    k = args.k
    ref_filename = args.ref_filename
    read_filename = args.read_filename
    output_folder = args.output_folder
    
    
    ref = read_fasta(ref_filename)
    read = read_fasta(read_filename)
    os.makedirs(output_folder, exist_ok=True)
    output_filename = read_filename + '_' + ref_filename
    output_filename = os.path.join(output_folder, output_filename)
    
    file1_d, file2_d, plot_d = get_data_for_output(read, ref, k)
    write_output_file_1(file1_d, output_filename)
    # plot_d also have the m_i_counts
    calculate_coverage(file2_d, plot_d, output_filename)
    get_plot(plot_d, output_filename)
    
if __name__ == '__main__':
    main()
    

    
