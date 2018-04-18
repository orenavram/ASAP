import re,regex, os
#import warnings
#warnings.filterwarnings('error', category=Bio.Seq.BiopythonWarning)

import matplotlib     # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay https://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay
import time
import Bio.Seq
import aa_sequences as aa
from plots_generator import generate_alignment_report_pie_chart, generate_mutations_boxplots
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
    chain_to_core_dna_to_Ka_Ks_dict = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))
    #chain_to_core_dna_to_num_of_non_synonymous_mutations = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))
    pseudo_count = 1
    # old code for SHM frequency distribution TODO: delete?
    #chain_to_num_of_mutations_to_counts_dict = dict(zip(allowed_chain_types, [{} for chain in allowed_chain_types]))

    #TODO: should the 'unknown' category be in all places where allowed_chain_types is used??
    chain_to_count_dict = dict.fromkeys(allowed_chain_types + ['unknown'], 0)
    isotypes_count_dict = dict.fromkeys(['A', 'A1', 'A2', 'D', 'E', 'G', 'M', 'unknown'], 0)
    errors_count_dict = dict.fromkeys(['cdr', 'len', 'quality', 'not_overlapped', 'inappropriate_end_j_seq', 'stop_codon'], 0)

    alignments_txt_path = os.path.join(mixcr_output_path, 'alignments.txt')
    alignments_filtered_txt_path = os.path.join(parsed_mixcr_output_path, 'alignments_filtered.txt')

    logger.info('Start parsing {}'.format(alignments_txt_path))
    with open(alignments_txt_path) as f:

        logger.info('File was opened succssefully.')
        # skip header-related variables- use to extract specified fields from alignments file
        alignments_filtered_txt = f.readline()

        logger.info('First line of file is:\n{}'.format(alignments_filtered_txt))

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

            #no more filtrations after this point!!
            alignments_filtered_txt += line

            if chain not in allowed_chain_types:
                logger.error("chain_type {} not in {}".format(chain, str(allowed_chain_types)))
                logger.error(line_tokens)
                chain = 'unknown'

            # update chain counts
            chain_to_count_dict[chain] += 1

            isotype = get_isotype(dna_read, core_dna, aa.end_j_seq)
            if chain != 'IGH':
                logger.info('Changing isotype from ' + isotype + 'to NONE')
                isotype = 'NONE'
                #TODO: should these be counted?
            else:
                # update isotype counts
                isotypes_count_dict[isotype] += 1

            # update aa_sequence counts
            sequences_frequency_counter[core_aa] = sequences_frequency_counter.get(core_aa, 0) + 1

            # set annotation for the (unique) aa_sequence (only for the first time)
            if core_aa not in chain_to_aa_read_to_meta_data_dict[chain]:
                chain_to_aa_read_to_meta_data_dict[chain][core_aa] = get_meta_data(line_tokens, chain, isotype, core_dna, core_aa, cdr3, best_v_family_col, best_d_family_col, best_j_family_col)

            # set annotation for unique aa_sequence (only for the first time)
            if core_dna not in chain_to_core_dna_to_Ka_Ks_dict[chain]:
                #extract mutations field from column number $best_v_alignment_col that looks like this:
                #1|292|312|21|313|SG5CI8ASG15CSA36CSG90ASA91GDC95I98GSC143TSC148ASC218TSC259A|1288.0
                mutations_field = line_tokens[best_v_alignment_col].split("|")[5]
                update_mutation_count(core_dna, mutations_field, chain_to_core_dna_to_Ka_Ks_dict[chain], pseudo_count) # chain_to_core_dna_to_num_of_non_synonymous_mutations[chain], pseudo_count)
                # old code for SHM frequency distribution TODO: delete?
                #update_mutation_count(core_dna, line_tokens[best_v_alignment_col], chain_to_num_of_mutations_to_counts_dict[chain], chain_to_core_dna_to_num_of_mutations[chain], chain_to_core_dna_to_num_of_non_synonymous_mutations[chain])

    #create nucleotide mutation frequency chart
    for chain in chain_to_core_dna_to_Ka_Ks_dict:
        core_dna_to_Ka_Ks_dict = chain_to_core_dna_to_Ka_Ks_dict[chain]
        if core_dna_to_Ka_Ks_dict != {}:
            mutations_file = parsed_mixcr_output_path + '/' + chain + mutation_count_file_suffix
            write_dict_to_file(mutations_file, core_dna_to_Ka_Ks_dict, value_type=list)
            #write_dict_to_file(parsed_mixcr_output_path + '/' + chain + '_non_synonymous' + mutation_count_file_suffix,
            #                   chain_to_core_dna_to_num_of_non_synonymous_mutations[chain])

            # old code for SHM frequency distribution TODO: delete?
            #write_dict_to_file(parsed_mixcr_output_path + '/' + chain + mutation_count_file_suffix, chain_to_num_of_mutations_to_counts_dict[chain])

    for chain in chain_to_aa_read_to_meta_data_dict:
        aa_read_to_meta_data_dict = chain_to_aa_read_to_meta_data_dict[chain]
        if aa_read_to_meta_data_dict != {}:
            with open(parsed_mixcr_output_path + '/' + chain + sequence_annotation_file_suffix, 'w') as f:
                for core_aa in aa_read_to_meta_data_dict:
                    f.write('\t'.join(aa_read_to_meta_data_dict[core_aa] + [str(sequences_frequency_counter[core_aa])]) + '\n')

    t2 = time.time()

    logger.debug('sum(isotypes_count_dict.values():' + str(sum(isotypes_count_dict.values())))

    outfile_report = parsed_mixcr_output_path + '/alignment_report.log'
    write_reports(outfile_report, t2 - t1, errors_count_dict, total_lines, chain_to_count_dict, isotypes_count_dict)

    outfile_pie_chart = outfile_report.replace('log', 'png')
    if isotypes_count_dict:
        generate_alignment_report_pie_chart(outfile_pie_chart, isotypes_count_dict, minimal_portion=0.0)

    #TODO: generate A_subisotypes pie chart

    #save alignments filtered file for debugging and fast different statistics....
    with open(alignments_filtered_txt_path, 'w') as f:
        f.write(alignments_filtered_txt)


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
def write_reports(out_file_path, total_time, errors_count_dict, total_lines, chain_to_count_dict, isotypes_count_dict):
    '''
    :param out_file_path:
    :param total_time:
    :param errors_count_dict:
    :param total_lines:
    :param chain_to_count_dict:
    :param isotypes_count_dict:
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
    A_sub_isotype_to_precent_dict = percent_calculator(A_sub_isotypes_to_count_dict, A_count)

    with open(out_file_path, 'w') as f:
        f.write('Total time (seconds) = {:.3f}\n'.format(total_time))
        f.write('Mixcr alignment file provides {} results\n'.format(total_lines))
        f.write('% of paired reads with non-sense CDR3 region = {:.3f}\n'.format(statistics_precent_dict['cdr']))
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


# old code for SHM frequency distribution TODO: delete?
#def update_mutation_count(dna, mutations_field, mutation_count_frequency_dict, dna_to_num_of_synonymous_mutations, dna_to_num_of_non_synonymous_mutations):
def update_mutation_count(dna, mutations_of_read, dna_to_num_of_mutations, pseudo_count):#, dna_to_num_of_non_synonymous_mutations, pseudo_count):
    # count number of nucleotide mutations
    # for more info regarding the field's format see the bottom of this page:
    # http://mixcr.readthedocs.io/en/latest/appendix.html#ref-encoding
    synonymous_mutation = non_synonymous_mutation = pseudo_count # pseudo counts to avoid zero division
    for match in re.finditer(r'S[ACGT](\d+)[ACGT]', mutations_of_read):
        if silent_mutation(match, dna):
            synonymous_mutation += 1
        else:
            non_synonymous_mutation += 1

    Ks_per_codon = synonymous_mutation/(len(dna)/3)
    Ka_per_codon = non_synonymous_mutation/(len(dna)/3)
    dna_to_num_of_mutations[dna] = [Ka_per_codon, Ks_per_codon]

    #old code for SHM frequency distribution TODO: delete?
    #positions_of_mutations = re.findall('\d+', mutations_of_read)
    #mutations_count = len(positions_of_mutations)
    #mutation_count_frequency_dict[mutations_count] = mutation_count_frequency_dict.get(mutations_count, 0) + 1


def silent_mutation(match, dna):
    mutation_position_in_read = int(match.group(1)) #e.g., catches the 207 in SA207G

    begining_of_derived_codon = 3*(mutation_position_in_read//3)
    end_of_derived_codon = 3*(mutation_position_in_read//3 + 1)
    derived_codon = dna[begining_of_derived_codon:end_of_derived_codon]

    mutation_position_in_codon = mutation_position_in_read%3
    ancestral_allele = match.group()[1] #e.g., catches the A in SA207G
    ancestral_codon = derived_codon[:mutation_position_in_codon] + ancestral_allele + derived_codon[mutation_position_in_codon+1:]

    return Bio.Seq.translate(ancestral_codon) == Bio.Seq.translate(derived_codon)