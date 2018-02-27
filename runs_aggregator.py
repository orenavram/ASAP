import os, subprocess
from text_handler import read_table_to_dict, write_dict_to_file
from plots_generator import plot_intersection_histogram
import logging
logger = logging.getLogger('main')


def join_runs_analyses(number_of_runs, run_output_paths, joint_parsed_mixcr_output_path, chains, sequence_annotation_file_suffix, mutation_count_file_suffix, minimal_overlap = 1):

    logger.info('run_number ' + str(number_of_runs))
    for chain in chains:
        sequence_to_entry_dict = {}
        aa_seq_intersection = {}
        mutation_counts_frequency ={}

        for i in range(len(run_output_paths[:-1])): #without joint
            run_output_path = run_output_paths[i]

            annotation_path = os.path.join(run_output_path, 'parsed_mixcr_output', chain + sequence_annotation_file_suffix)
            add_annotation_to_seqs_dict(sequence_to_entry_dict, aa_seq_intersection, annotation_path, i, number_of_runs)

            mutation_counts_frequency_path = os.path.join(run_output_path, 'parsed_mixcr_output', chain + mutation_count_file_suffix)
            add_counts_to_mutations_dict(mutation_counts_frequency, mutation_counts_frequency_path)


        #generate joint files and plot for each chain
        if sequence_to_entry_dict != {}:
            joint_annotation_path = os.path.join(joint_parsed_mixcr_output_path, chain + sequence_annotation_file_suffix)
            with open(joint_annotation_path, 'w') as f:
                for aa_seq in sequence_to_entry_dict:
                    if passes_overlap_minimal_threshold_in_all_runs(aa_seq_intersection[aa_seq], minimal_overlap):
                        f.write('\t'.join(sequence_to_entry_dict[aa_seq]) + '\n')
                    '''
                    else: # for debugging
                        #if passes_overlap_minimal_threshold_in_all_runs(aa_seq_intersection[aa_seq], 2):
                        logger.debug('Omitting', aa_seq, aa_seq_intersection[aa_seq], 'from joint annotation file')
                    '''
            final_fasta_path = os.path.join(run_output_paths[-1], chain + '_final.fasta')
            generate_final_fasta(aa_seq_intersection, final_fasta_path, minimal_overlap, sequence_to_entry_dict)

            generate_intersection_plot(number_of_runs, joint_annotation_path, sequence_annotation_file_suffix)

        #generate joint mutation counts file for each chain
        if mutation_counts_frequency != {}:
            with open(os.path.join(joint_parsed_mixcr_output_path, chain + mutation_count_file_suffix), 'w') as f:
                for mutation_count in mutation_counts_frequency:
                    f.write(str(mutation_count) + '\t' + str(int(mutation_counts_frequency[mutation_count]/len(run_output_paths))) + '\n')


def generate_final_fasta(aa_seq_intersection, final_fasta_path, minimal_overlap, sequence_to_entry_dict):
    with open(final_fasta_path, 'w') as f:
        i = 0
        for aa_seq in sequence_to_entry_dict:
            if passes_overlap_minimal_threshold_in_all_runs(aa_seq_intersection[aa_seq], minimal_overlap):
                i += 1
                aa_read, chain, cdr3, v_type, d_type, j_type, dna_read, isotype, read_frequency = sequence_to_entry_dict[aa_seq]
                f.write('|'.join(['> {}'.format(i), '{}({})'.format(chain, isotype), cdr3, v_type, d_type, j_type,
                                  ','.join(str(x) for x in aa_seq_intersection[aa_seq])]))
                f.write('\n' + aa_seq + '\n')


def add_counts_to_mutations_dict(mutation_counts_frequency, mutation_counts_frequency_path):
    if not os.path.exists(mutation_counts_frequency_path):
        logger.info('Can\'t find mutation frequency path:' + mutation_counts_frequency_path)
        return
    current_run_frequency = read_table_to_dict(mutation_counts_frequency_path)
    for key in current_run_frequency:
        mutation_counts_frequency[key] = mutation_counts_frequency.get(key, 0) + int(current_run_frequency[key])


def add_annotation_to_seqs_dict(sequence_to_entry_dict, aa_seq_intersection, annotation_path, i, number_of_runs):
    if not os.path.exists(annotation_path):
        logger.info('Can\'t find annotation path:' + annotation_path)
        return
    with open(annotation_path) as f:
        for line in f:
            '''
            QVQLQESGPGLVKPSETLSLTCTVSGGSISSYYWSWIRQPAGKGLEWIGRIYTSGSTNYNPSLKSRVTMSVDTSKNQFSLKLSSVTAADTAVYYCARTYSGSYYGRFDYWGQGTLVTVSS	IGH	CARTYSGSYYGRFDYW	IGHV4	IGHD1	IGHJ4*02	CCCAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCGCCGGGAAGGGACTGGAGTGGATTGGGCGTATCTATACCAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATGTCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGAACTTATAGTGGGAGCTACTACGGGCGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCCTCTCTCTCTCTCCTCCAGA	13
            '''
            tokens = line.rstrip().split('\t')
            aa_seq = tokens[0]
            count = int(tokens[-1])
            if aa_seq not in aa_seq_intersection:
                aa_seq_intersection[aa_seq] = [0] * number_of_runs  # initialize on first time
            aa_seq_intersection[aa_seq][i] = count
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
    percent_of_intersected_annotations_per_run = [0]*number_of_runs
    runs = []
    num_of_joint_annotations = int(subprocess.check_output(['wc', '-l', joint_annotation_path]).split()[0]) #number of lines in annotation file
    logger.debug('joint' + str(num_of_joint_annotations))
    for i in range(number_of_runs): #without joint
        run = 'run' + str(i+1)
        runs.append(run)
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



