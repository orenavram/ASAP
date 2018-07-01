import os
import re

# import warnings
# warnings.filterwarnings('error', category=Bio.Seq.BiopythonWarning)
import time
import Bio.Seq
import aa_sequences as aa
import regex
from plots_generator import generate_alignment_report_pie_chart
from text_handler import write_dict_to_file, logger, write_mapping_file
from auxiliaries import measure_time

def parse_alignment_file(mixcr_output_path, parsed_mixcr_output_path, sequence_annotation_file_suffix, mutations_file_suffix, len_threshold, qlty_threshold):
    '''parse alignment procedure'''
    # input: alignments file, path for output files, length threshold, quality thresholds of total sequence and of CDR3 region
    # output: none. creates output files as specified in "notes" file

    # column indices of the relevant data from mixcr's output (for more details see 'alignments.txt' file)
    overlapped_reads= 0
    quality = 1
    accession_number = 2
    DNA_FR1 = 4
    DNA_FR4 = 10
    AA_FR1 = 11
    AA_CDR3 = 16
    AA_FR4 = 17
    best_v_family = 23
    best_d_family = 24
    best_j_family = 25
    best_v_alignment = 31


    # dictionary to convert ASCII code to quality values
    ascii_to_quality_dict = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, '+': 10,
                             ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                             '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                             '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
                             'J': 41, 'K': 42}

    t1 = time.time()

    allowed_chain_types = ['IGH', 'IGK', 'IGL', 'unknown'] #do not use += or append here!

    total_lines = 1 # Shifted by one because of the header. More convenient when looking in notepad++...
    sequences_frequency_counter: {str: int} = {}

    # don't use dict.fromkeys here. Causes a BUG!!!
    chain_to_aa_read_to_meta_data_dict = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))
    chain_to_core_dna_to_mutations_info_dict = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))
    chain_to_core_aa_to_dna_reads_and_accession_numbers = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))
    pseudo_count = 1

    chain_to_count_dict = dict.fromkeys(allowed_chain_types, 0)
    isotypes_count_dict = dict.fromkeys(['A', 'A1', 'A2', 'D', 'E', 'G', 'M', 'unknown'], 0)
    errors_count_dict = dict.fromkeys(['no_overlap', 'too_short_length', 'too_low_quality', 'missing_cdr3', 'nonsense_stop_codon', 'inappropriate_end_j_seq'], 0)

    alignments_txt_path = os.path.join(mixcr_output_path, 'alignments.txt')
    #alignments_filtered_txt_path = os.path.join(parsed_mixcr_output_path, 'alignments_filtered.txt')

    logger.info('Start parsing {}'.format(alignments_txt_path))
    with open(alignments_txt_path) as f:

        logger.info('File was opened succssefully.')
        # skip header-related variables- use to extract specified fields from alignments file
        header = f.readline()
        logger.info('First line of file is:\n{}'.format(header))

        #alignments_filtered_txt = header

        #iterate over alignments file line by line            
        for line in f:

            #logger.debug('Next line of file is:\n{}'.format(line))
            line_tokens = line.split('\t')
            #count total number of entries provided by mixcr alignment
            total_lines += 1

            if not total_lines%100000:
                logger.info('total_lines: {}'.format(total_lines))
            # If the first token contains two sequences (separated by a comma) it means that
            # MiXCR was unable to find an overlap between the two paired-end reads.
            if ',' in line_tokens[overlapped_reads]:
                errors_count_dict['no_overlap'] += 1
                continue

            chain = line_tokens[best_v_family][:3]

            # sanity check
            # if line_tokens[20][:3] != line_tokens[best_v_family_col][:3]:
            #     logger.debug(line)
            #     logger.debug(line[20][:3])
            #     logger.debug(line[best_v_family_col][:3])
            #     logger.debug('line[20][:3] != line[best_v_family_col][:3]')

            dna_read = line_tokens[overlapped_reads]
            # a combination that should generate the relevant part of the antibody dna
            # (from the end of the 5' primer until the end of the end_j_seq)
            #core_dna = line_tokens[6] + line_tokens[4] + line_tokens[10] + line_tokens[8] + line_tokens[14] + line_tokens[12] + line_tokens[16]
            core_dna = ''.join(line_tokens[DNA_FR1:DNA_FR4+1])

            # discard too short core dna's
            read_len = len(core_dna)
            if read_len < len_threshold:
                errors_count_dict['too_short_length'] += 1
                continue

            # discard low quality reads
            sequencing_quality = line_tokens[quality]
            #calculate average quality of read
            average_quality = sum([ascii_to_quality_dict[k] for k in sequencing_quality]) / read_len
            if average_quality < qlty_threshold:
                errors_count_dict['too_low_quality'] += 1
                continue

            #verify CDR3 is present
            cdr3 = line_tokens[AA_CDR3]
            if cdr3 == '': # or '*' in cdr3 :
                errors_count_dict['missing_cdr3'] += 1
                continue

            # this should be the translation of the core_dna
            #core_aa = line_tokens[7] + line_tokens[5] + line_tokens[11] + line_tokens[9] + line_tokens[15] + line_tokens[13] + line_tokens[17]
            core_aa = ''.join(line_tokens[AA_FR1:AA_FR4+1])

            # fixed_core_dna = ''
            # partial_ORF = False
            # for i in range(DNA_FR1,DNA_FR4+1):
            #     fixed_core_dna += line_tokens[i]
            #     if len(line_tokens[i]) % 3 != 0: # fill partial ORFs (Probably a bug in MiXCR)
            #         partial_ORF = True
            #         fixed_core_dna += get_end_of_ORF_after_core_dna_fragment(dna_read, line_tokens[i])
            # if partial_ORF:
            #     # fix dna only if ORF is partial and translate to aa
            #     core_dna = fixed_core_dna
            #     core_aa = Bio.Seq.translate(core_dna)

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
            if not core_aa.endswith(aa.end_j_seq):
                logger.debug('end_j_seq after fixation is: {}'.format(core_aa[-len(aa.end_j_seq[0]): ]))

            # verify that core_aa is not non-sense
            if '*' in core_aa:
                logger.debug('line {} in alignment.txt file: STOP codon in core_aa!!!\n{}'.format(total_lines, core_aa))
                errors_count_dict['nonsense_stop_codon'] += 1
                continue

            # verify that there is a proper end_j_seq.
            # MUST be after making sure that '*' is NOT in core_aa (otherwise it makes problems with the regex).
            if chain == 'IGH':
                has_end_j_seq = False
                for end_j_seq in aa.end_j_seq: # aa.end_j_seq is a tuple with at least one string
                    if match_with_up_to_k_mismatches(core_aa[-len(end_j_seq): ], end_j_seq):
                        has_end_j_seq = True
                        break
                    if match_with_up_to_k_mismatches(core_aa[-len(end_j_seq)-1: -1], end_j_seq):
                        # in case of 'VTVSS_', remove last (full/partial) "codon"
                        core_aa = core_aa[:-1]
                        end_j_seq = core_aa[-len(end_j_seq):]  # update current_end_j_seq. Maybe it's with one mismatch
                        core_dna = trim_ORF(core_dna, end_j_seq)  # sometimes it's a full codon, sometimes partial
                        has_end_j_seq = True
                        break

                if not has_end_j_seq:
                    logger.debug('IGH with no end_j_seq in core_aa:\n{}'.format(core_aa))
                    errors_count_dict['inappropriate_end_j_seq'] += 1
                    # if core_aa[-6:-1] in aa.end_j_seq:
                    #     errors_count_dict['VTVSS_'] = errors_count_dict.get('VTVSS_',0) + 1
                    continue

            #no more filtrations after this point!!
            #alignments_filtered_txt += line

            if chain not in allowed_chain_types:
                logger.error('chain_type {} not in {}'.format(chain, allowed_chain_types))
                logger.error(line_tokens)
                chain = 'unknown'

            # update chain counts
            chain_to_count_dict[chain] += 1


            if chain == 'IGH':
                isotype = get_isotype(dna_read, core_dna, end_j_seq)
                # update isotype counts
                isotypes_count_dict[isotype] += 1
            else:
                #No need to count these
                isotype = 'NONE'

            # update aa_sequence counts
            sequences_frequency_counter[core_aa] = sequences_frequency_counter.get(core_aa, 0) + 1

            # set annotation for the (unique) aa_sequence (only for the first time)
            if core_aa not in chain_to_aa_read_to_meta_data_dict[chain]:
                chain_to_aa_read_to_meta_data_dict[chain][core_aa] = get_meta_data(line_tokens, chain, isotype, core_dna, core_aa, cdr3, best_v_family, best_d_family, best_j_family)

            # update mutation counts and Ka_Ks for the (unique) dna_sequence (only for the first time)
            if core_dna not in chain_to_core_dna_to_mutations_info_dict[chain]:
                #extract mutations field from column number $best_v_alignment_col that looks like this:
                #1|292|312|21|313|SG5CI8ASG15CSA36CSG90ASA91GDC95I98GSC143TSC148ASC218TSC259A|1288.0
                mutations_field = line_tokens[best_v_alignment].split("|")[5]
                update_mutation_count(core_dna, mutations_field, chain_to_core_dna_to_mutations_info_dict[chain], pseudo_count) # chain_to_core_dna_to_num_of_non_synonymous_mutations[chain], pseudo_count)

            # track mapping between each aa sequence and the reads behind it
            chain_to_core_aa_to_dna_reads_and_accession_numbers[chain][core_aa] = chain_to_core_aa_to_dna_reads_and_accession_numbers[chain].get(core_aa, []) + [(core_dna, line_tokens[accession_number])]

    # for chain in chain_to_core_dna_to_num_of_mutations:
    #     core_dna_to_num_of_mutations = chain_to_core_dna_to_num_of_mutations[chain]
    #     if core_dna_to_num_of_mutations != {}:
    #         mutation_counts_file = parsed_mixcr_output_path + '/' + chain + mutations_file_suffix
    #         write_dict_to_file(mutation_counts_file, core_dna_to_num_of_mutations)

    for chain in allowed_chain_types:
        core_dna_to_mutations_info_dict = chain_to_core_dna_to_mutations_info_dict[chain]
        if core_dna_to_mutations_info_dict != {}:
            mutations_info_file = parsed_mixcr_output_path + '/' + chain + mutations_file_suffix
            write_dict_to_file(mutations_info_file, core_dna_to_mutations_info_dict, value_type=list, header='dna' + '\t' + ';'.join(['Ka_per_codon', 'Ks_per_codon', 'number_of_baspair_mutations']))

    #for chain in chain_to_aa_read_to_meta_data_dict:
        aa_read_to_meta_data_dict = chain_to_aa_read_to_meta_data_dict[chain]
        if aa_read_to_meta_data_dict != {}:
            with open(parsed_mixcr_output_path + '/' + chain + sequence_annotation_file_suffix, 'w') as f:
                f.write('\t'.join(['chain', 'isotype', 'dna', 'aa', 'missing_cdr3', 'v_type', 'd_type', 'j_type', 'counts']) + '\n')
                for core_aa in aa_read_to_meta_data_dict:
                    f.write('\t'.join(aa_read_to_meta_data_dict[core_aa] + [str(sequences_frequency_counter[core_aa])]) + '\n')

        core_aa_to_dna_reads_and_accession_numbers = chain_to_core_aa_to_dna_reads_and_accession_numbers[chain]
        if core_aa_to_dna_reads_and_accession_numbers != {}:
            aa_to_read_and_accession_path = os.path.join(parsed_mixcr_output_path, chain + '_AA_to_DNA_reads.fasta')
            write_mapping_file(core_aa_to_dna_reads_and_accession_numbers, aa_to_read_and_accession_path)

    t2 = time.time()

    logger.debug('sum(isotypes_count_dict.values():' + str(sum(isotypes_count_dict.values())))

    outfile_report = parsed_mixcr_output_path + '/alignment_report.log'
    write_reports(outfile_report, t1, t2, errors_count_dict, total_lines, chain_to_count_dict, isotypes_count_dict)

    outfile_pie_chart = outfile_report.replace('log', 'png')
    if isotypes_count_dict:
        generate_alignment_report_pie_chart(outfile_pie_chart, isotypes_count_dict)

    #save alignments filtered file for debugging and fast different statistics....
    #with open(alignments_filtered_txt_path, 'w') as f:
    #    f.write(alignments_filtered_txt)


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


def trim_ORF(core_dna, end_j_seq):
    for i in range(1,4): #1, 2, 3
        if end_j_seq == Bio.Seq.translate(core_dna[len(core_dna)-3*len(end_j_seq)-i: len(core_dna)-i]):
            return core_dna[:len(core_dna)-i]
    logger.error(f'\n{"#"*50}\nSomething went wrong with the ORF!!\n{core_dna}\n{Bio.Seq.transcribe(core_dna)}\n{end_j_seq}\n{"#"*50}')
    #raise ValueError


def get_end_of_ORF_after_core_dna_fragment(dna_read, core_dna):
    '''
    used when the length of core_dna is not a multiplication of 3.
    :param dna_read:
    :param core_dna: a substring of the dna_read
    :return: next one (or two) nucleotide(s) in dna_read after core_dna (in order to fix the length of core_dna)
    '''
    end_of_core_dna_in_dna_read = dna_read.find(core_dna) + len(core_dna)
    dna_read_after_core_dna = dna_read[end_of_core_dna_in_dna_read:]
    return dna_read_after_core_dna[:3-(len(core_dna)%3)]


def match_with_up_to_k_mismatches(sub_string, string, max_mismatches_allowed = 1):#, mandatory = None):
    '''
    :param sub_string:
    :param string:
    :param max_mismatches_allowed:
    :param mandatory:
    :return: return a match of substring in string with up to k mismatches where position $mandatory must be a match
    '''
    if '*' in sub_string:
        # prevent problems with '*' (stop codons) in substring
        sub_string = sub_string.replace('*','X')
    match = regex.match('(' + sub_string + '){s<=' + str(max_mismatches_allowed) + '}', string)
    return match #and (mandatory == None or string[mandatory] == sub_string[mandatory])


    
# '''trim prefix to stop codon'''
# def edit_stop_codon(read_AA):
#
#     stop_codon_index = read_AA.rfind('*')       #last * in sequence
#
#     #there is a stop codon in the sequence
#     if stop_codon_index > -1:
#         read_AA = read_AA[stop_codon_index+1:]      #trim sequence start (up until last '*')
#
#         R_index = read_AA.find('R')
#         K_index = read_AA.find('K')
#
#         #find first occurance of R/K and start sequence from it
#         if R_index > -1 and K_index > -1:
#             new_start_index = min (R_index, K_index)
#         elif R_index > -1 or K_index > -1:
#             new_start_index = max(R_index, K_index)
#
#         #safety check- no 'R' and no 'K' left in sequence
#         else:
#             new_start_index = 0
#
#         read_AA = read_AA[new_start_index:]
#
#     return read_AA
      
    
def get_meta_data(line_tokens, chain, isotype, core_dna, core_aa, cdr3, best_v_family_col, best_d_family_col, best_j_family_col):

    if not re.match(chain + 'V\d+', line_tokens[best_v_family_col]):
        logger.error(line_tokens)
        logger.error(line_tokens[best_v_family_col])
    v_type = re.match(chain + 'V\d+', line_tokens[best_v_family_col]).group()
    d_type = 'unknown'
    if line_tokens[best_d_family_col]: # d assignment is sometimes missing
        d_type = re.match(chain + 'D\d+', line_tokens[best_d_family_col]).group()
    j_type = re.match(chain + 'J\d+', line_tokens[best_j_family_col]).group()

    return [chain, isotype, core_dna, core_aa, cdr3, v_type, d_type, j_type]


'''write report for procedure (statistics)'''
#input: out file path, time measurements, counters information
#output: none. writing all gathered data to report file
def write_reports(out_file_path, start_time, end_time, errors_count_dict, total_lines, chain_to_count_dict, isotypes_count_dict):
    '''
    :return: write information about the valid data that were extracted from mixcr's alignments.txt file
    '''
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
    A_sub_isotypes_to_precent_dict = percent_calculator(A_sub_isotypes_to_count_dict, A_count)

    with open(out_file_path, 'w') as f:
        f.write(f'Parsing alignment.txt took {measure_time(int(end_time-start_time))}\n')
        f.write('Mixcr alignment file provided {} results\n\n'.format(total_lines))

        # to keep the filtrations order...
        errors = ['no_overlap', 'too_short_length', 'too_low_quality', 'missing_cdr3', 'nonsense_stop_codon',
                  'inappropriate_end_j_seq']#, 'VTVSS_']
        f.write('Detailed filtrations:\n')
        for error in errors:
            f.write('{} entries with {} ({:.3f}%)\n'.format(errors_count_dict[error], error, statistics_precent_dict[error]))
        f.write('\n')

        f.write('Chains dispersion:\n')
        for chain in sorted(chain_to_count_dict):
            f.write('{} of {} chains ({:.3f}%)\n'.format(chain_to_count_dict[chain], chain, chain_to_percent_dict[chain]))
        f.write('\n')

        f.write('Isotypes dispersion:\n')
        for isotype in isotypes_count_dict:
            f.write('{} of isotype {} ({:.3f}%)\n'.format(isotypes_count_dict[isotype], isotype, isotype_to_precent_dict[isotype]))
        f.write('\n')

        f.write('"A" sub-isotypes dispersion:\n')
        for isotype in A_sub_isotypes_to_precent_dict:
            f.write('{} of isotype {} ({:.3f}%)\n'.format(isotypes_count_dict[isotype], isotype, isotype_to_precent_dict[isotype]))
        f.write('\n')

    
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


def update_mutation_count(dna, mutations_of_read, dna_to_Ka_Ks_per_codon, pseudo_count):
    # count number of nucleotide mutations
    # for more info regarding the field's format see the bottom of this page:
    # http://mixcr.readthedocs.io/en/latest/appendix.html#ref-encoding
    synonymous_mutation = non_synonymous_mutation = pseudo_count # pseudo counts to avoid zero division
    for match in re.finditer(r'S[ACGT](\d+)[ACGT]', mutations_of_read):
        if silent_mutation(match, dna):
            synonymous_mutation += 1
        else:
            non_synonymous_mutation += 1

    Ka_per_codon = non_synonymous_mutation/(len(dna)/3)
    Ks_per_codon = synonymous_mutation/(len(dna)/3)
    num_of_substitutions = synonymous_mutation + non_synonymous_mutation
    dna_to_Ka_Ks_per_codon[dna] = [Ka_per_codon, Ks_per_codon, num_of_substitutions]


def silent_mutation(match, dna):
    mutation_position_in_read = int(match.group(1)) #e.g., catches the 207 in SA207G

    begining_of_derived_codon = 3*(mutation_position_in_read//3)
    end_of_derived_codon = 3*(mutation_position_in_read//3 + 1)
    derived_codon = dna[begining_of_derived_codon:end_of_derived_codon]

    mutation_position_in_codon = mutation_position_in_read%3
    ancestral_allele = match.group()[1] #e.g., catches the A in SA207G
    ancestral_codon = derived_codon[:mutation_position_in_codon] + ancestral_allele + derived_codon[mutation_position_in_codon+1:]

    return Bio.Seq.translate(ancestral_codon) == Bio.Seq.translate(derived_codon)