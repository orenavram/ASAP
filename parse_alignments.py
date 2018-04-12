import re,regex, os
#import warnings
#warnings.filterwarnings('error', category=Bio.Seq.BiopythonWarning)

import matplotlib     # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay https://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay
import time
import Bio.Seq
from plots_generator import generate_alignment_report_pie_chart
import aa_sequences as aa
from text_handler import write_dict_to_file, logger
import traceback


def parse_alignment_file(mixcr_output_path, parsed_mixcr_output_path, sequence_annotation_file_suffix, mutation_count_file_suffix, len_threshold, qlty_threshold, allowed_chain_types = ['IGH', 'IGK', 'IGL']):
    '''parse alignment procedure'''
    # input: alignments file, path for output files, length threshold, quality thresholds of total sequence and of CDR3 region
    # output: none. creates output files as specified in "notes" file

    # column indices of the relevant data from mixcr's output (for more details see 'alignments.txt' file)
    dna_sequence = 0
    quality_col = 1
    # description_col = 2
    AA_seq_CDR3_col = 13
    best_v_family_col = 23
    best_d_family_col = 24
    best_j_family_col = 25
    # dna_mutation_col = [27, 29, 31]
    # aa_mutation_col = [28, 30, 32]
    best_v_alignment_col = 33

    # num_of_CDR3_regions=3
    # num_of_CDR3_regions = len(dna_mutation_col)

    # dictionary to convert ASCII code to quality values
    ascii_to_quality_dict = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, '+': 10,
                             ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                             '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                             '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
                             'J': 41, 'K': 42}

    t1 = time.time()

    total_lines = 1 # Shifted by one because of the header. More convenient when looking in notepad++...
    sequences_frequency_counter = {}

    # don't use dict.fromkeys here. Causes a BUG!!!
    chain_to_aa_read_to_meta_data_dict = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))
    chain_to_num_of_mutations_to_counts_dict = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))

    chain_to_count_dict = dict.fromkeys(allowed_chain_types + ['unknown'], 0)
    isotypes_count_dict = dict.fromkeys(['A', 'A1', 'A2', 'D', 'E', 'G', 'M', 'unknown'], 0)
    errors_count_dict = dict.fromkeys(['cdr', 'len', 'quality', 'not_overlapped', 'inappropriate_end_j_seq', 'stop_codon'], 0)

    alignments_txt_path = os.path.join(mixcr_output_path, 'alignments.txt')

    logger.info('Start parsing {}'.format(alignments_txt_path))
    with open(alignments_txt_path) as f:

        logger.info('File was opened succssefully.')
        # skip header-related variables- use to extract specified fields from alignments file
        logger.info('First line of file is:\n{}'.format(f.readline()))

        #iterate over alignments file line by line            
        for line in f:

            #logger.debug('Next line of file is:\n{}'.format(line))
            line_tokens = line.split('\t')
            #count total number of entries provided by mixcr alignment
            total_lines += 1

            # If the first token contains two sequences (separated by a comma) it means that
            # MiXCR was unable to find an overlap between the two paired-end reads.
            if ',' in line_tokens[0]:
                errors_count_dict['not_overlapped'] += 1
                continue

            # verify CDR3 is present and in an ORF
            cdr3 = line_tokens[AA_seq_CDR3_col]
            if cdr3 == '' or '*' in cdr3 :
                errors_count_dict['cdr'] += 1
                continue
                
            # discard too short reads
            dna_read = line_tokens[dna_sequence]
            read_len = len(dna_read)
            if read_len < len_threshold:
                errors_count_dict['len'] += 1
                continue

            # discard low quality reads
            sequencing_quality = line_tokens[quality_col]

            #calculate average quality of read
            average_quality = sum([ascii_to_quality_dict[k] for k in sequencing_quality]) / read_len
            if average_quality < qlty_threshold:
                errors_count_dict['quality'] += 1
                continue

            chain = line_tokens[best_v_family_col][:3]

            # sanity check
            if line_tokens[20][:3] != line_tokens[best_v_family_col][:3]:
                logger.debug(line)
                logger.debug(line[20][:3])
                logger.debug(line[best_v_family_col][:3])
                logger.debug('line[20][:3] != line[best_v_family_col][:3]')

            if chain not in allowed_chain_types:
                logger.error("chain_type {} not in {}".format(chain, str(allowed_chain_types)))
                logger.error(line_tokens)
                chain = 'unknown'

            #update chain counter
            chain_to_count_dict[chain] += 1

            # a combination that should generate the relevant part of the antibody dna
            # (from the end of the 5' primer until the end of the end_j_seq)
            core_dna = line_tokens[6] + line_tokens[4] + line_tokens[10] + line_tokens[8] + line_tokens[14] + line_tokens[12] + line_tokens[16]

            # this should be the translation of the core_dna
            core_aa = line_tokens[7] + line_tokens[5] + line_tokens[11] + line_tokens[9] + line_tokens[15] + line_tokens[13] + line_tokens[17]

            if len(line_tokens[16]) % 3 != 0: # fill partial ORF. Probably a bug in MiXCR
                core_dna += get_end_of_ORF_after_end_j_seq(dna_read, line_tokens[16])
                core_aa = Bio.Seq.translate(core_dna)

            #sanity check
            if core_aa != Bio.Seq.translate(core_dna):
                logger.debug('core_aa is NOT identical to the translated core_dna')
                logger.debug('core_aa:\n{}'.format(core_aa))
                logger.debug('translated core_dna:\n{}'.format(Bio.Seq.translate(core_dna)))

            #sanity check
            if (core_dna not in dna_read) and (not core_aa.endswith('VTVS_')):
                logger.debug('dna_read:\n{}'.format(dna_read))
                logger.debug('core dna:\n{}'.format(core_dna))

            #sanity check
            if core_aa[-len(aa.end_j_seq): ] != aa.end_j_seq:
                logger.debug('end_j_seq after fixation is: {}'.format(core_aa[-len(aa.end_j_seq): ]))

            if not match_with_up_to_k_mismatches(core_aa[-len(aa.end_j_seq): ], aa.end_j_seq):
                logger.error('No end_j_seq in core_aa:\n{}'.format(core_aa))
                errors_count_dict['inappropriate_end_j_seq'] += 1
                continue

            if '*' in core_aa:
                logger.debug('line {} in alignment.txt file: STOP codon in core_aa!!!\n{}'.format(total_lines, core_aa))
                errors_count_dict['stop_codon'] += 1
                continue

            isotype = get_isotype(dna_read, core_dna, aa.end_j_seq)
            isotypes_count_dict[isotype] += 1

            sequences_frequency_counter[core_aa] = sequences_frequency_counter.get(core_aa, 0) + 1

            if core_aa not in chain_to_aa_read_to_meta_data_dict[chain]:
                chain_to_aa_read_to_meta_data_dict[chain][core_aa] = get_meta_data(line_tokens, core_aa, core_dna, cdr3, isotype, best_v_family_col, best_d_family_col, best_j_family_col)

            update_mutation_count(line_tokens[best_v_alignment_col], chain_to_num_of_mutations_to_counts_dict[chain])

    #create nucleotide mutation frequency chart
    for chain in allowed_chain_types:
        write_dict_to_file(parsed_mixcr_output_path + '/' + chain + mutation_count_file_suffix,
                           chain_to_num_of_mutations_to_counts_dict[chain])

    for chain in chain_to_aa_read_to_meta_data_dict:
        aa_read_to_meta_data_dict = chain_to_aa_read_to_meta_data_dict[chain]
        if aa_read_to_meta_data_dict != {}:
            with open(parsed_mixcr_output_path + '/' + chain + sequence_annotation_file_suffix, 'w') as f:
                for core_aa in aa_read_to_meta_data_dict:
                    f.write('\t'.join(aa_read_to_meta_data_dict[core_aa] + [str(sequences_frequency_counter[core_aa])]) + '\n')

    t2 = time.time()

    logger.debug('sum(isotypes_count_dict.values():' + str(sum(isotypes_count_dict.values())))
    #assert total_iso == sum(isotypes_count_dict.values())
    #assert total_chains == sum(chain_to_count_dict.values())

    outfile_report = parsed_mixcr_output_path + '/alignment_report.log'
    write_reports(outfile_report, t2 - t1, errors_count_dict, total_lines, chain_to_count_dict, isotypes_count_dict)


# '''identify relevant file according to chain and translate read to amino acid sequence'''
# #input: chain type, 3 chain-specific files, amino acid sequence
# #output: file of input chain and chain isotype (if exists)
# def translate_dna_to_amino_acid_and_get_isotype_if_any(isotypes_count_dict, read, chain_type, allowed_chain_types):
#
#     global core_aa
#     isotype = ''
#
#     if chain_type not in allowed_chain_types:
#         #sanity check. Shouldn't crash here!
#         logger.error('chain type: ' + chain_type)
#         raise Exception('Unrecognized chain in alignments file!')
#
#     #recognize chain
#     if (chain_type == "IGH"):
#         isotype, read_AA = get_isotype_and_translation(read, isotypes_count_dict)
#
#     elif (chain_type == "IGK"):
#
#         read_AA = Bio.Seq.translate(read[2:])
#         # check_reads_length_divides_by_3(read[2:], read_AA)
#
#         #edit end of IGK sequence (keep E/DIK before "RTVAA")
#         start_of_match = get_start_index_of_substring_in_string(aa.IGK_id, read_AA, 1)
#         if start_of_match != None:
#             read_AA = read_AA[:start_of_match-1] + "K" #ensure that IGK ends with K
#
#     elif (chain_type == "IGL"):
#         read_AA = Bio.Seq.translate(read[2:])
#         #check_reads_length_divides_by_3(read[2:], read_AA)
#
#         #edit end of IGL sequence (add "GQPK" after "VTVLL")
#         start_of_match = get_start_index_of_substring_in_string(aa.IGL_id, read_AA, 1)
#         if start_of_match != None:
#             read_AA = read_AA[:start_of_match] + aa.IGL_id + "GQPK"
#
#     return read_AA, isotype
            
def get_isotype(dna_read, core_dna, end_j_seq):

    '''if isinstance(end_j_seq, tuple):

        len_end_j_seq
    else:
        len_end_j_seq = len(end_j_seq)
    '''
    dna_of_end_j_seq = core_dna[-3*len(end_j_seq):] # dna is 3 times longer than aa

    dna_after_end_j_seq = dna_read[dna_read.find(dna_of_end_j_seq) + len(dna_of_end_j_seq):]

    aa_after_end_j_seq = Bio.Seq.translate(dna_after_end_j_seq)

    if match_with_up_to_k_mismatches(aa.isotype_A_id, aa_after_end_j_seq) != None:
        if match_with_up_to_k_mismatches(aa.isotype_A1_id, aa_after_end_j_seq) != None:
            return 'A1'
        elif match_with_up_to_k_mismatches(aa.isotype_A2_id, aa_after_end_j_seq) != None:
            return 'A2'
        else:
            return 'A'
    elif match_with_up_to_k_mismatches(aa.isotype_D_id, aa_after_end_j_seq) != None:
        return 'D'
    elif match_with_up_to_k_mismatches(aa.isotype_E_id, aa_after_end_j_seq) != None:
        return 'E'
    elif match_with_up_to_k_mismatches(aa.isotype_G_id, aa_after_end_j_seq) != None:
        return 'G'
    elif match_with_up_to_k_mismatches(aa.isotype_M_id, aa_after_end_j_seq) != None:
        return 'M'

    return 'unknown'


def get_end_of_ORF_after_end_j_seq(dna_read, core_dna):
    '''
    used when the length of core_dna is not a multiplication of 3.
    :param dna_read:
    :param core_dna: a substring of the dna_read
    :return: next one (or two) nucleotide(s) in dna_read after core_dna (in order to fix the length of core_dna)
    '''
    end_of_core_dna_in_dna_read = dna_read.find(core_dna) + len(core_dna)
    dna_read_after_core_dna = dna_read[end_of_core_dna_in_dna_read:]
    return dna_read_after_core_dna[:3-(len(core_dna)%3)]

# '''recognize isotype by markers'''
##input: amino acid sequence
##output: isotype of sequence
# def get_isotype_and_translation(read, isotypes_count_dict):
#
#     #global total_iso
#     #total_iso += 1
#
#     isotype = ''
#
#     for i in range(0,3): #count isotype in ALL reading frames
#
#         aa_read = Bio.Seq.translate(read[i:])
#
#         #check_reads_length_divides_by_3(read[i:], read_AA)
#         start = max(0, aa_read.find(aa.end_j_seq)) # -1 if no end_j_seq
#         if get_start_index_of_substring_in_string(aa.isotype_A_id, aa_read[start:], 2) != None:
#             if get_start_index_of_substring_in_string(aa.isotype_A1_id, aa_read[start:], 2, 2) != None:
#                 isotype = 'A1'
#                 aa_read = get_final_aa_seq(aa_read, aa.isotype_A1_id, aa.mass_spec_seq, 2, 1, aa.end_j_seq)
#             elif get_start_index_of_substring_in_string(aa.isotype_A2_id, aa_read[start:], 2, 2) != None:
#                 isotype = 'A2'
#                 aa_read = get_final_aa_seq(aa_read, aa.isotype_A2_id, aa.mass_spec_seq, 2, 1, aa.end_j_seq)
#             else:
#                 isotype = 'A'
#                 aa_read = get_final_aa_seq(aa_read, aa.isotype_A_id, aa.mass_spec_seq, -1, 1, aa.end_j_seq)
#
#         elif get_start_index_of_substring_in_string(aa.isotype_D_id, aa_read[start:], 2) != None:
#             isotype = 'D'
#             aa_read = get_final_aa_seq(aa_read, aa.isotype_D_id, aa.mass_spec_seq, -1, 1, aa.end_j_seq)
#
#         elif get_start_index_of_substring_in_string(aa.isotype_E_id, aa_read[start:], 2) != None:
#             isotype = 'E'
#             aa_read = get_final_aa_seq(aa_read, aa.isotype_E_id, aa.mass_spec_seq, -1, 1, aa.end_j_seq)
#
#         elif get_start_index_of_substring_in_string(aa.isotype_G_id, aa_read[start:], 2) != None:
#             isotype = 'G'
#             aa_read = get_final_aa_seq(aa_read, aa.isotype_G_id, aa.mass_spec_seq, -1, 1, aa.end_j_seq)
#
#         elif get_start_index_of_substring_in_string(aa.isotype_M_id, aa_read[start:], 2) != None:
#             isotype = 'M'
#             aa_read = get_final_aa_seq(aa_read, aa.isotype_M_id, aa.mass_spec_seq, -1, 1, aa.end_j_seq)
#
#         if isotype:
#             break #an isotype was detected!
#
#     has_end_j_seq = False
#     #edit IGH sequences with unknown isotype
#     if not isotype:
#         for i in range(0, 3):  # count isotype in ALL reading frames
#             aa_read = Bio.Seq.translate(read[i:])
#             start_of_match = get_start_index_of_substring_in_string(aa.end_j_seq, aa_read, 1)
#             if start_of_match != None:
#                 aa_read = aa_read[:start_of_match] + aa.end_j_seq + aa.mass_spec_seq
#                 isotype = 'unknown'
#                 has_end_j_seq=True
#                 break
#
#         if not has_end_j_seq:
#             aa_read = ''
#             logger.info('The following read was not recognized as a subisotype and has no end_j_seq:')
#             logger.info(read)
#             logger.info('No aa_read to return...')
#             logger.info('Frame 0:\n' + Bio.Seq.translate(read[0:]))
#             logger.info('Frame 1:\n' + Bio.Seq.translate(read[1:]))
#             logger.info('Frame 2:\n' + Bio.Seq.translate(read[2:]))
#
#     '''
#     if total_iso != sum(isotypes_count_dict.values()):
#         logger.info(total_iso)
#         logger.info(sum(isotypes_count_dict.values()))
#         logger.info(isotypes_count_dict)
#         raise AssertionError
#     '''
#
#     if isotype:
#         isotypes_count_dict[isotype] += 1
#
#     return isotype, aa_read


# def check_reads_length_divides_by_3(read, read_AA):
#     if len(read) % 3 != 0:
#         logger.info('read: ' + read)
#         logger.info('read length: ' + str(len(read)))
#         logger.info('read_AA: ' + read_AA)
#         logger.info('read_AA length: ' + str(len(read_AA)))


def match_with_up_to_k_mismatches(sub_string, string, max_mismatches_allowed = 1, mandatory = -1):
    match = regex.match('(' + sub_string + '){s<=' + str(max_mismatches_allowed) + '}', string)
    return match and (mandatory<0 or string[mandatory] == sub_string[mandatory])



'''find (partial )sub-string in string'''
#find substring(iso id) in another string(read_AA) 
#allowing "allow_mismatch" miss matches as long as their not in index "must_match" index 
def get_start_index_of_substring_in_string(sub_string, string, max_mismatches_allowed, mandatory = -1):

    len_sub_string = len(sub_string)
    for i in range(len(string) - len_sub_string + 1):
        potential_match = string[i: i+len_sub_string]
        match = regex.match('(' + sub_string + '){s<=' + str(max_mismatches_allowed) + '}', potential_match)
        if not match: # in case there's no match
            continue
        if mandatory<0 or potential_match[mandatory] == sub_string[mandatory]:
            return i

    return None



# '''edit amino acid sequence to mass-spec'''
# #find IGH sequence identifier and edit the end of the sequence
# def get_final_aa_seq(read_AA, seq_id, mass_spec_seq, must_match, mismatch, end_J_seq):
#
#     start_of_match = get_start_index_of_substring_in_string(end_J_seq, read_AA, 1)
#
#     if start_of_match != None: #IGH identifier found
#         read_AA = read_AA[:start_of_match + len(end_J_seq)] + mass_spec_seq # discards seq_id (downstream to end_j_seq)
#     else: #IGH identifier not found => edit chain according to isotype identifier
#         start_of_match = get_start_index_of_substring_in_string(seq_id, read_AA, mismatch, must_match)
#         read_AA = read_AA[:start_of_match] + mass_spec_seq
#
#     return read_AA
    

# '''remove adapter sequence PPLIP from start of the sequence'''
# def trim_adapter_sequence_if_any(dna_read, aa_read, adapter_seq):
#
#     start_of_match = get_start_index_of_substring_in_string(adapter_seq, aa_read, 1)
#     if start_of_match != None:
#         length_to__trim = start_of_match+len(adapter_seq)
#         aa_read = aa_read[length_to__trim:]
#         dna_read = dna_read[length_to__trim * 3:]
#
#     return dna_read, aa_read

    
'''trim prefix to stop codon'''
def edit_stop_codon(read_AA):
    
    stop_codon_index = read_AA.rfind('*')       #last * in sequence

    #there is a stop codon in the sequence
    if stop_codon_index > -1:
        read_AA = read_AA[stop_codon_index+1:]      #trim sequence start (up until last '*')
        
        R_index = read_AA.find('R')
        K_index = read_AA.find('K')
        
        #find first occurance of R/K and start sequence from it
        if R_index > -1 and K_index > -1: 
            new_start_index = min (R_index, K_index)
        elif R_index > -1 or K_index > -1:
            new_start_index = max(R_index, K_index)
            
        #safety check- no 'R' and no 'K' left in sequence
        else:
            new_start_index = 0
        
        read_AA = read_AA[new_start_index:]
            
    return read_AA
      
    
'''write entry to FASTA file'''
#input: output file, read id, read Nuc sequence, headers, alignment entry, read AA sequence, isotype
#output: none. writing entry information to its proper file 
def get_meta_data(line, aa_read, dna_read, cdr3, isotype, best_v_family_col, best_d_family_col, best_j_family_col):

    # verify legal pairing- by comparing reads' ids
    # slice from start until this seqeunce ' 1:N:0:1'
    # line[description_col][0:-8]

    chain = line[best_v_family_col][:3]

    if chain != 'IGH':
        if not isotype:
            logger.info('Changing isotype from ' + isotype + 'to NONE')
        isotype = 'NONE'

    if not re.match(chain + 'V\d+', line[best_v_family_col]):
        print(line)
        print(line[best_v_family_col])
    v_type = re.match(chain + 'V\d+', line[best_v_family_col]).group()
    d_type = 'unknown'
    if line[best_d_family_col]: # d assignment is sometimes missing
        d_type = re.match(chain + 'D\d+', line[best_d_family_col]).group()
    j_type = re.match(chain + 'J\d+', line[best_j_family_col]).group()


    return [aa_read, chain, cdr3, v_type, d_type, j_type, dna_read, isotype]


'''write report for procedure (statistics)'''
#input: out file path, time measurements, counters information
#output: none. writing all gathered data to report file
def write_reports(out_file_path, total_time, errors_count_dict, total_lines, chain_to_count_dict, isotypes_count_dict):

    statistics_precent_dict = percent_calculator(errors_count_dict, total_lines)
    chain_to_percent_dict = percent_calculator(chain_to_count_dict, sum(chain_to_count_dict.values()))
    isotype_to_precent_dict = percent_calculator(isotypes_count_dict, sum(isotypes_count_dict.values()))

    logger.info('chain_to_count_dict: {}'.format(chain_to_count_dict))
    logger.info('chain_to_percent_dict: {}'.format(chain_to_percent_dict))
    logger.info('errors_count_dict: {}'.format(errors_count_dict))
    logger.info('statistics_precent_dict: {}'.format(statistics_precent_dict))
    logger.info('isotypes_count_dict: {}'.format(isotypes_count_dict))

    A_sub_isotypes = ['A', 'A1', 'A2']
    A_sub_isotypes_to_count_dict = dict(
        [(A_sub_isotype, isotypes_count_dict[A_sub_isotype]) for A_sub_isotype in A_sub_isotypes])
    A_count = isotypes_count_dict['A1'] + isotypes_count_dict['A2'] + isotypes_count_dict['A']
    A_sub_isotype_to_precent_dict = percent_calculator(A_sub_isotypes_to_count_dict, A_count)
    #A_precent = A_sub_isotype_to_precent_dict['A1'] + A_sub_isotype_to_precent_dict['A2'] + A_sub_isotype_to_precent_dict['A']

    with open(out_file_path, 'w') as f:
        #information about valid data that was extracted from mixcr alignments file


        f.write('Total time (seconds) = {:.3f}\n'.format(total_time))
        f.write('Mixcr alignment file provides {} results\n'.format(total_lines))
        f.write('% of paired reads with no CDR3 region = {:.3f}\n'.format(statistics_precent_dict['cdr']))
        f.write('% of entries shorter than specified threshold = {:.3f}\n'.format(statistics_precent_dict['len']))
        f.write('% of entries with lower quality than specified threshold = {:.3f}\n'.format(statistics_precent_dict['quality']))
        f.write('% of entries that are not overlapped = {:.3f}\n'.format(statistics_precent_dict['not_overlapped']))
        f.write('% of entries that have an inappropriate end_j_seq = {:.3f}\n'.format(statistics_precent_dict['inappropriate_end_j_seq']))
        f.write('% of entries that have a stop codon = {:.3f}\n\n'.format(statistics_precent_dict['stop_codon']))

        for chain in sorted(chain_to_count_dict):
            f.write('% of {} chains = {:.3f}\n'.format(chain, chain_to_percent_dict[chain]))

        f.write('\n')
        logger.debug('% of A isotype = {:.3f}\n'.format(isotype_to_precent_dict['A']))
        f.write('% of A isotype = {:.3f}\n'.format(isotype_to_precent_dict['A']))
        f.write('% of A1 isotype = {:.3f}\n'.format(isotype_to_precent_dict['A1']))
        f.write('% of A2 isotype = {:.3f}\n'.format(isotype_to_precent_dict['A2']))
        f.write('% of D isotype = {:.3f}\n'.format(isotype_to_precent_dict['D']))
        f.write('% of E isotype = {:.3f}\n'.format(isotype_to_precent_dict['E']))
        f.write('% of G isotype = {:.3f}\n'.format(isotype_to_precent_dict['G']))
        f.write('% of M isotype = {:.3f}\n'.format(isotype_to_precent_dict['M']))
        f.write('% of unknown isotypes = {:.3f}\n'.format(isotype_to_precent_dict['unknown']))

        f.write('\n')
        f.write('% of A isotype (out of A) = {:.3f}\n'.format(A_sub_isotype_to_precent_dict['A']))
        f.write('% of A1 isotype (out of A) = {:.3f}\n'.format(A_sub_isotype_to_precent_dict['A1']))
        f.write('% of A2 isotype (out of A) = {:.3f}\n'.format(A_sub_isotype_to_precent_dict['A2']))

    outfile_pie_chart = out_file_path.replace('log', 'png')
    if isotypes_count_dict:
        generate_alignment_report_pie_chart(outfile_pie_chart, isotypes_count_dict)

    #TODO: generate A_subisotypes pie chart

    
'''calculate perectage for report file'''
def percent_calculator(fractions, complete):

    percentage = dict.fromkeys(fractions, 0)

    if complete != 0:
        for key in fractions:
            percentage[key] = 100 * fractions[key]/complete


    if sum(fractions.values()) != complete:
        logger.debug('complete: ' + str(complete))
        logger.debug('fractions: ' + str(fractions))
        logger.debug('sum(fractions.values()): ' + str(sum(fractions.values())))
        logger.warn('Sum of fractions does not sum up to complete!', exc_info=True)

    return percentage



def update_mutation_count(best_v_alignment, mutation_count_frequency_dict):
    # count number of nucleotide mutations
    # for more info regarding the field's format see the bottom of this page:
    # http://mixcr.readthedocs.io/en/latest/appendix.html#ref-encoding
    mutations_of_read = best_v_alignment.split("|")[5]
    positions_of_mutations = re.findall('\d+', mutations_of_read)
    mutations_count = len(positions_of_mutations)
    mutation_count_frequency_dict[mutations_count] = mutation_count_frequency_dict.get(mutations_count, 0) + 1

'''
def plot_nuc_mut_cnt_dict(sample, chain, out_file_path, nuc_mut_cnt_dict, MUT_BY_CDR3):
    'plot mutation count frequency'
    max_mut_cnt = max(nuc_mut_cnt_dict)  # maximum number of mutations
    nuc_mut_cnt_values = range(max_mut_cnt + 1)  # list of numbers up until the maximum

    # pad with zeros mutations amounts that are absent
    for i in range(max_mut_cnt):
        if i not in nuc_mut_cnt_dict:
            nuc_mut_cnt_dict[i] = 0

    mut_sum = sum(nuc_mut_cnt_dict.values())

    # cretae list of frequencies of mutations amount
    nuc_mut_freq = list()
    for nuc_mut_cnt in sorted(nuc_mut_cnt_dict):
        crnt_freq = round((nuc_mut_cnt_dict[nuc_mut_cnt] / mut_sum) * 100, 3)
        nuc_mut_freq.append(crnt_freq)

    # keep raw data in xls file
    if MUT_BY_CDR3:
        xls_file_name = out_file_path + "/" + chain + "_mutation_freq_by_cdr3.xls"
    else:
        xls_file_name = out_file_path + "/" + chain + "_mutation_freq_by_seq.xls"

    with open(xls_file_name, "w") as f:
        f.write("# Mutations\t# reads\t# Mutations frequency\n")
        for i in nuc_mut_cnt_values:
            f.write('%d\t%d\t%f\n' % (nuc_mut_cnt_values[i], nuc_mut_cnt_dict[i], nuc_mut_freq[i]))

    # set plot parameters
    pos = np.arange(len(nuc_mut_cnt_values))
    width = 1.0

    ax = plt.axes()
    ax.set_xticks(pos + (width / 2))
    ax.set_xticklabels(nuc_mut_cnt_values)

    plt.figure()
    plt.title(sample + " " + chain + " Frequency of nucleotide mutations amount")
    plt.xlabel("Number of nucleotide mutations")
    plt.ylabel("Frequency (%)")
    plt.bar(pos, nuc_mut_freq, width, color='b')
    plt.tight_layout()

    if MUT_BY_CDR3:
        fig_name = out_file_path + "/" + chain + "_mutation_freq_by_cdr3.png"
    else:
        fig_name = out_file_path + "/" + chain + "_mutation_freq_by seq.png"

    plt.savefig(fig_name)
    plt.close('all')
'''