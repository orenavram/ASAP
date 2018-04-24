from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align.Generic import Alignment
from Bio.Alphabet import Gapped, IUPAC
from msa_parser import remove_sparse_columns
from weblogo_generator import generate_weblogo
from plots_generator import plot_barplot, generate_clonal_expansion_histogram
from text_handler import string_similarity, read_table_to_dict
from subprocess import check_output
from directory_creator import create_dir
from logging import getLogger
logger = getLogger('main')

import os

def analyze_cdr3(gp):
    plot_cdr3_length_distributions(gp)
    analyze_top_clones(gp)



def plot_cdr3_length_distributions(gp):
    for cdr3_analysis_dir in gp.cdr3_analysis_paths:
        for chain in gp.chains:
            cdr3_len_file = os.path.join(cdr3_analysis_dir, chain + '_cdr3_len_counts.' + gp.raw_data_file_suffix)
            if not os.path.exists(cdr3_len_file):
                logger.info('Skipping: {} (file does not exists)'.format(cdr3_len_file))
                continue
            plot_path = cdr3_len_file.replace(gp.raw_data_file_suffix, 'png')
            logger.info('Plotting cdr3_len_file: ' + cdr3_len_file)
            d = read_table_to_dict(cdr3_len_file, key_type=int, value_type=int, skip_rows=1)
            plot_barplot(d, plot_path, x_label='\nLength of CDR3')


def analyze_top_clones(gp):
    #top clones analysis
    k = gp.top_cdr3_clones_to_further_analyze
    for i in range(gp.number_of_runs+1): # +1 for the joint analysis
        for chain in gp.chains:
            annotations_path = os.path.join(gp.parsed_mixcr_output_paths[i], chain + gp.sequence_annotation_file_suffix)
            if not os.path.exists(annotations_path):
                logger.info('No such file: {}. Skipping its analysis...'.format(annotations_path))
                continue

            cdr3_to_aa_reads, cdr3_to_counts = parse_sequence_annotations_file(annotations_path, skip_rows=1)
            cdr3_to_num_of_different_reads = dict((cdr3, len(cdr3_to_aa_reads[cdr3])) for cdr3 in cdr3_to_aa_reads)

            most_common_cdr3 = sorted(cdr3_to_counts, key=lambda cdr3:(cdr3_to_counts.get(cdr3), cdr3_to_num_of_different_reads.get(cdr3)), reverse=True)

            cdr3_annotations_path = os.path.join(gp.cdr3_analysis_paths[i], chain + gp.cdr3_annotation_file_suffix)
            write_cdr3_counts_file(cdr3_annotations_path, most_common_cdr3, cdr3_to_counts, cdr3_to_num_of_different_reads)

            clonal_expansion_histogram_path = cdr3_annotations_path.replace(gp.raw_data_file_suffix, 'png')
            generate_clonal_expansion_histogram(cdr3_annotations_path, clonal_expansion_histogram_path, gp.top_cdr3_clones_to_clonal_expansion_graph)

            analyze_most_k_common_cdr3(annotations_path, gp.cdr3_analysis_paths[i], most_common_cdr3, cdr3_to_aa_reads, cdr3_to_counts, chain, gp, k)


def write_cdr3_counts_file(cdr3_annotations_path, most_common_cdr3, cdr3_to_counts, cdr3_to_num_of_different_reads):
    sum_of_cdr3_counts = sum(cdr3_to_counts.values())
    sum_of_cdr3_different_reads = sum(cdr3_to_num_of_different_reads.values())
    with open(cdr3_annotations_path, 'w') as f:
        for cdr3 in most_common_cdr3:
            cdr3_counts = cdr3_to_counts[cdr3]
            cdr3_different_reads = cdr3_to_num_of_different_reads[cdr3]
            f.write('\t'.join([cdr3, str(cdr3_counts), str(cdr3_different_reads), str(cdr3_counts / sum_of_cdr3_counts),
                               str(cdr3_different_reads / sum_of_cdr3_different_reads)]) + '\n')


def analyze_most_k_common_cdr3(annotations_path, cdr3_analysis_dir, most_common_cdr3, cdr3_to_aa_reads, cdr3_to_counts, chain, gp, k):
    most_k_common_cdr3_to_entry = {} #:{str:[str,str]}
    most_k_common_cdr3 = most_common_cdr3[:k]
    most_k_common_cdr3_to_counts = {cdr3: cdr3_to_counts[cdr3] for cdr3 in most_k_common_cdr3}
    most_k_common_cdr3_to_aa_reads = {cdr3: cdr3_to_aa_reads[cdr3] for cdr3 in most_k_common_cdr3}
    top_clones_path = os.path.join(cdr3_analysis_dir, chain + '_top_{}_clones'.format(k))
    create_dir(top_clones_path)
    if most_k_common_cdr3 != []:
        for i in range(gp.top_cdr3_clones_to_further_analyze):
            cdr3 = most_k_common_cdr3[i]
            multiple_sequences = most_k_common_cdr3_to_aa_reads[cdr3]
            cluster_prefix = os.path.join(top_clones_path, 'cluster_' + str(i))
            ms_path = cluster_prefix + '_ms.fasta'
            msa_path = cluster_prefix + '_msa.aln'

            try:
                aln, msa = align_sequences(multiple_sequences, ms_path, msa_path, i)
            except:
                print(multiple_sequences)
                raise
            summary_align = AlignInfo.SummaryInfo(aln)

            # find consensus sequence by majority rule
            consensus = summary_align.dumb_consensus(threshold=0, ambiguous='#')
            consensus = str(consensus)

            # find most similar sequence
            aa_most_similar_to_consesnsus, similarity_rate = find_most_similar_sequence(msa, consensus)

            #try:
            dna_most_similar_to_consesnsus = find_correspnding_dna(aa_most_similar_to_consesnsus, annotations_path)
            #except:
            #    pass

            # add entry to dict
            most_k_common_cdr3_to_entry[cdr3] = [cdr3, str(most_k_common_cdr3_to_counts[cdr3]),
                                                 str(len(most_k_common_cdr3_to_aa_reads[cdr3])), consensus,
                                                 aa_most_similar_to_consesnsus, dna_most_similar_to_consesnsus,
                                                 str(similarity_rate)]

            # generate web logo
            weblogo_file_path = cluster_prefix + '_weblogo.pdf'
            try:
                generate_weblogo(msa_path, weblogo_file_path)
            except Exception as e:
                logger.warning('Couldn\'t generate weblogo {} due to the following reason: {}'.format(msa_path, e.args))

        # write cdr3_to_entry to file
        top_cdr3_extended_annotations_path = os.path.join(cdr3_analysis_dir, chain + gp.top_cdr3_annotation_file_suffix)
        with open(top_cdr3_extended_annotations_path, 'w') as f:
            #TODO: CHECK THAT THESE HEADERS AND DATA ARE CORRECT!!
            f.write('\t'.join(['cdr3', 'cdr3_counts',  'aa_counts', 'consensus', 'aa_most_similar_to_consesnsus', 'dna_most_similar_to_consesnsus', 'similarity_rate']))
            for cdr3 in most_k_common_cdr3:
                f.write('\t'.join(most_k_common_cdr3_to_entry[cdr3]) + '\n')
                # write_dict_to_file(cdr3_annotations_path, most_k_common_cdr3_to_entry, sort_by=most_k_common_cdr3_to_counts.get, reverse=True)


def parse_sequence_annotations_file(annotations_path, skip_rows=0):
    cdr3_to_counts = {} #:{str:int}
    cdr3_to_aa_reads = {} #:{str:[str]}

    with open(annotations_path) as f:
        for i in range(skip_rows):
            f.readline()
        for line in f:
            chain, isotype, dna_read, aa_read, cdr3, v_type, d_type, j_type, read_frequency = line.split('\t')
            cdr3_to_counts[cdr3] = cdr3_to_counts.get(cdr3, 0) + int(read_frequency)
            cdr3_to_aa_reads[cdr3] = cdr3_to_aa_reads.get(cdr3, []) + [aa_read]

    return cdr3_to_aa_reads, cdr3_to_counts


def find_most_similar_sequence(msa, consensus):
    max_similarity = 0
    max_seq = ''
    for seq in msa:
        similarity = round(string_similarity(seq, consensus), 2)
        if similarity > max_similarity:
            max_similarity = similarity
            max_seq = seq
    return max_seq, max_similarity


def find_correspnding_dna(aa_most_similar_to_consesnsus, annotations_path):
    # assuming that the dna sequence appears at the 2nd position in the annotation_file
    '''
    b'IGH\tM\tGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t301\nIGH\tM\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t61\nIGH\tM\tGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tEVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t33\nIGH\tM\tGAGGTGCAGCTGGTGGAGACTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tEVQLVETGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t31\nIGH\tM\tGAGGTGCAGCTGTTGGAGCCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tEVQLLEPGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t4\nIGH\tM\tGAGGTGCAGCTGTTGGAGACTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tEVQLLETGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t24\nIGH\tM\tGAGGTGCAGCTGGTGGAGCCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAACTATTAGTGGTAGTGGTGGTAGCACATTCTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATATATTACTGTGTCAAGGTGTCCGGGGGTTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\tEVQLVEPGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSTISGSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCVKVSGGFDYWGQGTLVTVSS\tCVKVSGGFDYW\tIGHV3\tIGHD3\tIGHJ4\t9\n'
    '''
    raw_dna_most_similar_to_consesnsus = check_output(['grep', aa_most_similar_to_consesnsus, annotations_path]).decode('utf-8').split()[2] #check_output returns bytes class object so it needs to be decoded
    dna_most_similar_to_consesnsus = str(raw_dna_most_similar_to_consesnsus)[2:-1] #without ^"b'" an $"'"
    logger.info('dna_most_similar_to_consesnsus is {}'.format(dna_most_similar_to_consesnsus))
    return str(dna_most_similar_to_consesnsus)


def align_sequences(multiple_sequences, ms_path, msa_path, cluster_number):
    # create alignment object from a list of multiple sequences and write it to a file

    aln = Alignment(Gapped(IUPAC.protein, '-'))

    # write ms to file
    with open(ms_path, 'w') as f:
        f.write('\n'.join('>'+str(i) + '\n' + multiple_sequences[i] for i in range(len(multiple_sequences))))

    if len(multiple_sequences) > 1:
        # align file
        logger.info('running mafft:')
        mafft_cmd = 'mafft ' + ms_path + ' > ' + msa_path
        logger.info(mafft_cmd)
        os.system(mafft_cmd)
        # remove columns with more than 20% gaps
        remove_sparse_columns(msa_path, msa_path, 0.2)
    else:
        logger.info('Only 1 sequence in {}. Skipping mafft and copying ms file as is.'.format(cluster_number))
        with open(msa_path, 'w') as f:
            f.write('\n'.join('>' + str(i) + '\n' + multiple_sequences[i] for i in range(len(multiple_sequences))))

    for record in SeqIO.parse(msa_path, 'fasta'):
        aln.add_sequence(record.id, str(record.seq))

    logger.debug('Alignment is:')
    msa_str = ''
    for record in aln._records:
        msa_str += record.seq._data + '\n'
    logger.debug(msa_str)

    return aln, msa_str.split()