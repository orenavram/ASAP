import logging
import os
import subprocess

import numpy as np
from plots_generator import plot_venn, plot_correlation, generate_alignment_report_pie_chart

from text_handler import read_table_to_dict, write_dict_to_file, write_mapping_file

logger = logging.getLogger('main')


def join_runs_analyses(number_of_runs, run_output_paths, joint_parsed_mixcr_output_path, chains, sequence_annotation_file_suffix, mutations_file_suffix, minimal_overlap = 1):

    logger.info('Number of replicates is: ' + str(number_of_runs))
    for chain in chains:
        sequence_to_entry_dict = {}
        aa_seq_to_counts_in_each_run_dict = {}
        shared_dna_reads = None

        for i in range(len(run_output_paths[:number_of_runs])): #without joint
            run_output_path = run_output_paths[i]

            annotation_path = os.path.join(run_output_path, 'parsed_mixcr_output', chain + sequence_annotation_file_suffix)
            if not os.path.exists(annotation_path):
                logger.info(f'No annotation path for run {i+1} chain {chain} (path does not exist: {annotation_path}). Continuing for the next run (if any).')
                continue
            add_annotation_to_seqs_dict(sequence_to_entry_dict, aa_seq_to_counts_in_each_run_dict, annotation_path, i, number_of_runs)

            core_dna_to_to_mutations_info_dict_path = os.path.join(run_output_path, 'parsed_mixcr_output', chain + mutations_file_suffix)
            current_run_core_dna_to_mutations_info_dict_path = read_table_to_dict(core_dna_to_to_mutations_info_dict_path, value_type=list, skip_rows=1)
            if shared_dna_reads == None:
                shared_dna_reads = set(current_run_core_dna_to_mutations_info_dict_path.keys())
            else:
                shared_dna_reads = shared_dna_reads.intersection(set(current_run_core_dna_to_mutations_info_dict_path.keys()))

        joint_core_dna_to_mutations_info_dict = dict((key, current_run_core_dna_to_mutations_info_dict_path[key]) for key in shared_dna_reads)


        #remove entries that do not pass the minimal count threshold in all runs
        #after that, these dictionaries keys will contain only joint aa sequences
        aa_seq_to_counts_with_less_than_minimal_threshold = filter_seqs_below_overlap_minimal_threshold(sequence_to_entry_dict, aa_seq_to_counts_in_each_run_dict, minimal_overlap)

        #generate joint files and plot for each chain
        if sequence_to_entry_dict != {}:
            joint_annotation_path = os.path.join(joint_parsed_mixcr_output_path, chain + sequence_annotation_file_suffix)
            isotypes_count_dict = {} #:{str:int}
            with open(joint_annotation_path, 'w') as f:
                for aa_seq in sequence_to_entry_dict:
                    entry = sequence_to_entry_dict[aa_seq]
                    f.write('\t'.join(entry) + '\n')
                    isotype = entry[1]
                    isotypes_count_dict[isotype] = isotypes_count_dict.get(isotype, 0) + 1

            if chain == 'IGH':
                outfile_pie_chart = os.path.join(joint_parsed_mixcr_output_path, 'alignment_report.png')
                generate_alignment_report_pie_chart(outfile_pie_chart, isotypes_count_dict, 'Joint')

            final_fasta_path = os.path.join(run_output_paths[-1], f'V{chain[-1]}_AA_sequences.fasta')
            generate_final_fasta(final_fasta_path, sequence_to_entry_dict, aa_seq_to_counts_in_each_run_dict)

            generate_aa_seq_to_dna_and_accession_mapping_file(chain, sequence_to_entry_dict, run_output_paths, number_of_runs)

            plot_runs_correlation(aa_seq_to_counts_in_each_run_dict, aa_seq_to_counts_with_less_than_minimal_threshold, number_of_runs, run_output_paths[-1], chain)

            generate_intersection_plot(number_of_runs, joint_annotation_path, sequence_annotation_file_suffix)

        #generate joint mutation counts file for each chain
        # if joint_mutation_counts_dict != {}:
        #     mutations_file = os.path.join(joint_parsed_mixcr_output_path, chain + mutations_file_suffix)
        #     write_dict_to_file(mutations_file, joint_mutation_counts_dict, sort_by=int)

        #generate joint mutation counts file for each chain
        if joint_core_dna_to_mutations_info_dict != {}:
            mutations_file = os.path.join(joint_parsed_mixcr_output_path, chain + mutations_file_suffix)
            write_dict_to_file(mutations_file, joint_core_dna_to_mutations_info_dict, value_type=list, header='\t'.join(['dna', 'Ka_per_codon', 'Ks_per_codon', 'number_of_baspair_mutations']))


def filter_seqs_below_overlap_minimal_threshold(sequence_to_entry_dict, aa_seq_to_counts_in_each_run, minimal_overlap):
    keys_to_remove = []
    for aa_seq in aa_seq_to_counts_in_each_run:
        if not passes_overlap_minimal_threshold_in_all_runs(aa_seq_to_counts_in_each_run[aa_seq], minimal_overlap):
            keys_to_remove.append(aa_seq)

    aa_seq_to_counts_with_less_than_minimal_threshold = {} #:{str:int}
    #two loops are necessary because mutating a dict while iterating on it raises an error
    for aa_seq in keys_to_remove:
        logger.debug('Omitting {} {} from joint analysis.'.format(aa_seq, aa_seq_to_counts_in_each_run[aa_seq]))
        sequence_to_entry_dict.pop(aa_seq)
        aa_seq_to_counts_with_less_than_minimal_threshold[aa_seq] = aa_seq_to_counts_in_each_run.pop(aa_seq)

    return aa_seq_to_counts_with_less_than_minimal_threshold


def plot_runs_correlation(aa_seq_to_counts_in_each_run, aa_seq_to_counts_with_less_than_minimal_threshold, number_of_runs, joint_path, chain):
    counts_vector_files = [0] * number_of_runs

    for i in range(number_of_runs):
        counts_vector_file = os.path.join(joint_path, 'run' + str(i+1) + '_'+ chain + '_counts_vector.txt')
        counts_vector_files[i] = open(counts_vector_file, 'w')

    # write counts vector of each run to a file
    for aa_seq in aa_seq_to_counts_in_each_run:
        counts = aa_seq_to_counts_in_each_run[aa_seq]
        for i in range(number_of_runs):
            counts_vector_files[i].write(str(counts[i]) + '\n')

    # write counts that are less than the minimal overlapping threshold
    for aa_seq in aa_seq_to_counts_with_less_than_minimal_threshold:
        counts = aa_seq_to_counts_with_less_than_minimal_threshold[aa_seq]
        for i in range(number_of_runs):
            counts_vector_files[i].write(str(counts[i]) + '\n')

    for i in range(number_of_runs):
        counts_vector_files[i].close()

    #for each runs pair compute correlation and plot it in the joint dir
    for i in range(number_of_runs-1):
        run_i_path = os.path.join(joint_path, 'run' + str(i + 1) + '_' + chain + '_counts_vector.txt')
        vector_i = read_vector_from_file(run_i_path)
        for j in range(i+1, number_of_runs):
            run_j_path = os.path.join(joint_path, 'run' + str(j + 1) + '_' + chain + '_counts_vector.txt')
            vector_j = read_vector_from_file(run_j_path)
            plot_path = os.path.join(joint_path, 'run' + str(i + 1) + '_' + 'run' + str(j + 1) + '_' + chain + '_correlation.png')
            plot_correlation(vector_i, vector_j, i+1, j+1, plot_path)


def read_vector_from_file(path):
    with open(path) as f:
        vector = f.read().split()
    return np.array([int(x) for x in vector])


def generate_final_fasta(final_fasta_path, sequence_to_entry_dict, aa_seq_to_counts_in_each_run=None):
    with open(final_fasta_path, 'w') as f:
        for aa_seq in sequence_to_entry_dict:
            chain, isotype, dna_read, aa_read, cdr3, v_type, d_type, j_type, read_frequency = sequence_to_entry_dict[aa_seq]
            if chain != 'IGH':
                d_type = 'NA'
            f.write('|'.join(['>' + chain, isotype, cdr3, v_type, d_type, j_type]))
            if aa_seq_to_counts_in_each_run:
                # add counts in each run only for the joint fasta!
                f.write('|' + str(sum(x for x in aa_seq_to_counts_in_each_run[aa_seq])))
                f.write('|' + ','.join(str(x) for x in aa_seq_to_counts_in_each_run[aa_seq]))
            else:
                f.write('|' + read_frequency)
            #from aa_sequences import mass_spec_seq
            #f.write('\n' + aa_seq + (mass_spec_seq if chain == 'IGH' else '') + '\n') # add mass_spec only to IGH!!
            f.write('\n' + aa_seq + '\n')


def add_counts_to_mutations_dict(joint_mutation_counts_frequency, mutation_counts_frequency_path):
    if not os.path.exists(mutation_counts_frequency_path):
        logger.info('Can\'t find mutation frequency path:' + mutation_counts_frequency_path)
        return
    current_run_frequency = read_table_to_dict(mutation_counts_frequency_path, value_type=int, skip_rows=1)
    for key in current_run_frequency:
        joint_mutation_counts_frequency[key] = joint_mutation_counts_frequency.get(key, 0) + current_run_frequency[key]


def generate_aa_seq_to_dna_and_accession_mapping_file(chain, sequence_to_entry_dict, run_output_paths, number_of_runs):
    joint_core_aa_to_dna_reads_and_accession_numbers = {}
    for i in range(len(run_output_paths[:number_of_runs])): #without joint
        run_output_path = run_output_paths[i]
        core_aa_to_dna_reads_and_accession_numbers = os.path.join(run_output_path, 'parsed_mixcr_output', f'{chain}_AA_to_DNA_reads.fasta')
        logger.info(f'Parsing {core_aa_to_dna_reads_and_accession_numbers} and adding to the joint file')
        with open(core_aa_to_dna_reads_and_accession_numbers) as f:
            ''' For example:
            >EVQLLESGAEVKKPGESLKISCKGSGYSFTSYWIGWVRQMPGKGLEWMGIIYPGDSDTRYSPSFQGQVTISADKSISTAYLQWSSLKASDTAMYYCARSRGGSYGTDYFDYWGQGTLVTVSS
            <M04473:15:000000000-B8Y9M:1:1101:21471:2444 1:N:0:1
            GAGGTGCAGCTGTTGGAGTCTGGAGCAGAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCTGGATACAGCTTTACCAGCTACTGGATCGGCTGGGTGCGCCAGATGCCCGGGAAAGGCCTGGAGTGGATGGGGATCATCTATCCTGGTGACTCTGATACCAGATACAGCCCGTCCTTCCAAGGCCAGGTCACCATCTCAGCCGACAAGTCCATCAGCACCGCCTACCTGCAGTGGAGCAGCCTGAAGGCCTCGGACACCGCCATGTATTACTGTGCGAGGTCCAGAGGCGGGAGCTACGGAACGGACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA
            <M04473:15:000000000-B8Y9M:1:2110:20973:16619 1:N:0:1
            GAGGTGCAGCTGTTGGAGTCTGGAGCAGAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCTGGATACAGCTTTACCAGCTACTGGATCGGCTGGGTGCGCCAGATGCCCGGGAAAGGCCTGGAGTGGATGGGGATCATCTATCCTGGTGACTCTGATACCAGATACAGCCCGTCCTTCCAAGGCCAGGTCACCATCTCAGCCGACAAGTCCATCAGCACCGCCTACCTGCAGTGGAGCAGCCTGAAGGCCTCGGACACCGCCATGTATTACTGTGCGAGGTCCAGAGGCGGGAGCTACGGAACGGACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA
            >QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARDLGPYDSSGYFQHWGQGTLVTVSS
            <M04473:15:000000000-B8Y9M:1:1101:22849:2507 1:N:0:1
            CAGGTCCAGCTTGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCATCTGGATACACCTTCACCAGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAATAATCAACCCTAGTGGTGGTAGCACAAGCTACGCACAGAAGTTCCAGGGCAGAGTCACCATGACCAGGGACACGTCCACGAGCACAGTCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGACCTGGGCCCCTATGATAGTAGTGGTTACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCA
            <M04473:15:000000000-B8Y9M:1:1102:15422:10777 1:N:0:1
            CAGGTCCAGCTTGTACAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCATCTGGATACACCTTCACCAGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAATAATCAACCCTAGTGGTGGTAGCACAAGCTACGCACAGAAGTTCCAGGGCAGAGTCACCATGACCAGGGACACGTCCACGAGCACAGTCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGACCTGGGCCCCTATGATAGTAGTGGTTACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCA
            '''
            aa_seq = ''
            while True:
                line = f.readline().rstrip()
                if not line:
                    break
                if line.startswith('>'):
                    aa_seq = line[1:]
                    if aa_seq not in sequence_to_entry_dict:
                        aa_seq = ''
                    continue
                if aa_seq:  # skip dna and accession for aa_seq that doesn't belong to the joint run
                    # reached here? aa_seq belongs to the joint run
                    accession = line[1:]  # line startswith '<'
                    dna_seq = f.readline().rstrip()
                    joint_core_aa_to_dna_reads_and_accession_numbers[aa_seq] = joint_core_aa_to_dna_reads_and_accession_numbers.get(aa_seq, []) + [(dna_seq, accession)]

    if joint_core_aa_to_dna_reads_and_accession_numbers != {}:
        joint_core_aa_to_read_and_accession_path = os.path.join(run_output_paths[-1], 'parsed_mixcr_output', f'{chain}_AA_to_DNA_reads.fasta')
        write_mapping_file(joint_core_aa_to_dna_reads_and_accession_numbers, joint_core_aa_to_read_and_accession_path)


def add_annotation_to_seqs_dict(sequence_to_entry_dict, aa_seq_to_counts_in_each_run_dict, annotation_path, i, number_of_runs):
    with open(annotation_path) as f:
        for line in f:
            '''
            QVQLQESGPGLVKPSETLSLTCTVSGGSISSYYWSWIRQPAGKGLEWIGRIYTSGSTNYNPSLKSRVTMSVDTSKNQFSLKLSSVTAADTAVYYCARTYSGSYYGRFDYWGQGTLVTVSS	IGH	CARTYSGSYYGRFDYW	IGHV4	IGHD1	IGHJ4*02	CCCAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCGCCGGGAAGGGACTGGAGTGGATTGGGCGTATCTATACCAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATGTCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGAACTTATAGTGGGAGCTACTACGGGCGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCCTCTCTCTCTCTCCTCCAGA	13
            '''
            tokens = line.rstrip().split('\t')
            if tokens[0] == 'chain':
                continue  # skip header
            aa_seq = tokens[3]
            count = int(tokens[-1])
            if aa_seq not in aa_seq_to_counts_in_each_run_dict:
                aa_seq_to_counts_in_each_run_dict[aa_seq] = [0] * number_of_runs  # initialize on first time
            aa_seq_to_counts_in_each_run_dict[aa_seq][i] = count
            if aa_seq not in sequence_to_entry_dict:
                sequence_to_entry_dict[aa_seq] = tokens
            else:
                # just update the tokens counter (last item)
                sequence_to_entry_dict[aa_seq][-1] = str(int(sequence_to_entry_dict[aa_seq][-1]) + count)


def passes_overlap_minimal_threshold_in_all_runs(l, threshold):
    for item in l:
        if item < threshold:
            return False
    return True


def generate_intersection_plot(number_of_runs, joint_annotation_path, sequence_annotation_file_suffix):
    runs_annotations_sets = []
    runs = []
    percent_of_intersected_annotations_per_run = []
    num_of_joint_annotations = int(subprocess.check_output(['wc', '-l', joint_annotation_path]).split()[0]) #number of lines in annotation file
    for i in range(number_of_runs):  # without joint
        run = 'run' + str(i + 1)
        runs.append(run)
        run_annotation_path = joint_annotation_path.replace('joint', run)
        with open(run_annotation_path) as f:
            lines = f.readlines()
        num_of_run_annotations = len(lines)
        logger.info('{} has {} annotations'.format(run, num_of_run_annotations))
        aa_seqs = [line.split()[3] for line in lines]
        runs_annotations_sets.append(set(aa_seqs))
        percent_of_intersected_annotations_per_run.append(num_of_joint_annotations / num_of_run_annotations * 100)

    logger.info(percent_of_intersected_annotations_per_run)
    raw_data_path = joint_annotation_path.replace(sequence_annotation_file_suffix, '_runs_intersections.txt')
    write_dict_to_file(raw_data_path, {(i+1):percent_of_intersected_annotations_per_run[i] for i in range(number_of_runs)}, header='\t'.join(['run_number', 'percent_of_intersected_annotations_with_the_joint']))

    out_path = raw_data_path.replace('txt', 'png')
    plot_venn(out_path, runs_annotations_sets, runs)


'''
def generate_intersection_plot(number_of_runs, joint_annotation_path, sequence_annotation_file_suffix):
    percent_of_intersected_annotations_per_run = [0]*number_of_runs
    runs = []
    num_of_joint_annotations = int(subprocess.check_output(['wc', '-l', joint_annotation_path]).split()[0]) #number of lines in annotation file
    logger.debug('joint' + str(num_of_joint_annotations))
    for i in range(number_of_runs): #without joint
        run = 'run' + str(i+1)
        runs.append(run)
        if num_of_joint_annotations != 0:
            run_annotation_path = joint_annotation_path.replace('joint', run)
            num_of_run_annotations = int(subprocess.check_output(['wc', '-l', run_annotation_path]).split()[0]) #number of lines in annotation file
            logger.info(run)
            logger.info(num_of_run_annotations)
            percent_of_intersected_annotations_per_run[i] = num_of_joint_annotations/num_of_run_annotations*100
    logger.info(percent_of_intersected_annotations_per_run)

    raw_data_path = joint_annotation_path.replace(sequence_annotation_file_suffix, '_runs_intersections.txt')
    write_dict_to_file(raw_data_path, {i:percent_of_intersected_annotations_per_run[i] for i in range(number_of_runs)})
    plot_path = raw_data_path.replace('txt', 'png')
    plot_intersection_histogram(runs, percent_of_intersected_annotations_per_run, plot_path)
'''

