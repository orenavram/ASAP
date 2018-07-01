from sys import argv
from Bio.Seq import IUPAC
from fastq_parser import FastqParser
import os
import argparse
import re
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

"""
Sai FID:
HHHHACHHHHACHHHNGCAG
Sai RID:
GHHHHHACAHHHHHACAHHHHN
"""

def handle_UMI(fastq_path, filtered_fastq_path, umi_iupac, fastq2_path, filtered_fastq2_path, log_file_path=None, relevant_records_path=None):

    pair_mate = '1' if 'R1' in fastq_path else '2'
    if not log_file_path:
        log_file_path = f'{os.path.split(filtered_fastq_path)[0]}/filter_by_R{pair_mate}.log'

    if not relevant_records_path:
        relevant_records_path = f'{os.path.split(filtered_fastq_path)[0]}/relevant_records_based_on_R{pair_mate}.log'

    filtrations_scenarios = ['no_umi', 'same_umi_different_read_no_major', 'same_read', 'same_umi_different_read']
    filtrations_count_dict = dict.fromkeys(filtrations_scenarios, 0)

    logger.info(f'Generating UMI to INFO dict for fastq file: {fastq_path}')
    umi_to_info, total_number_of_reads = process_fastq_file(fastq_path, umi_iupac, filtrations_count_dict)

    logger.info(f"Filtering reads by UMI.\nTotal number of different UMIs is: {len(umi_to_info)}\nTotal number of reads that didn't contain {umi_iupac} UMI is: {filtrations_count_dict['no_umi']}")
    filter_info_by_UMIs(umi_to_info, filtrations_count_dict)

    logger.info('Writing filtered fastq file...')
    sorted_umi_by_rec_id = sorted(umi_to_info, key=lambda umi: umi_to_info[umi][2])
    sorted_rec_nums = [umi_to_info[umi][2] for umi in sorted_umi_by_rec_id]
    with open(relevant_records_path, 'w') as f:
        f.write('\n'.join(str(x) for x in sorted_rec_nums))
    result = ''
    for umi in sorted_umi_by_rec_id:
        result += '\n'.join(umi_to_info[umi][0]) + '\n'
    with open(filtered_fastq_path, 'w') as f:
        f.write(result)
    if fastq2_path:
        logger.info(f'Applying filtrations over paired-end mate file: {fastq2_path}')
        handle_fastq2(fastq2_path, filtered_fastq2_path, sorted_rec_nums)

    logger.info(f'Writing filtrations statistics to: {log_file_path}')
    write_log_file(log_file_path, filtrations_scenarios, filtrations_count_dict, total_number_of_reads)


def write_log_file(log_file_path, filtrations_scenarios, filtrations_count_dict, total_number_of_reads):
    result = f'Total number of reads BEFORE filtration: {total_number_of_reads} (100.000%)\n\n'
    i = 0
    for scenario in filtrations_scenarios:
        result += f'Number of reads {"filtered" if i<2 else "collapsed"} due to {scenario}: {filtrations_count_dict[scenario]} ({100*filtrations_count_dict[scenario]/total_number_of_reads:.3f}%)\n'
        i += 1
    number_of_records_left_after_filtration = total_number_of_reads - sum(filtrations_count_dict.values())
    result += f'\nTotal number of reads AFTER filtration: {number_of_records_left_after_filtration} ({100*number_of_records_left_after_filtration/total_number_of_reads:.3f}%)\n'
    with open(log_file_path, 'w') as f:
        f.write(result)


def filter_info_by_UMIs(umi_to_info, filtrations_count_dict):
    umis = list(umi_to_info)
    for umi in umis:
        info = umi_to_info[umi]
        supporting_reads_to_counts = info[1]
        if len(supporting_reads_to_counts) == 1:  # count collapsed reads with the same umi and same sequence
            filtrations_count_dict['same_read'] += list(supporting_reads_to_counts.values())[0] - 1  # we still keep 1
        else:
            collapse_reads_with_same_umi_by_major_sequence(supporting_reads_to_counts, filtrations_count_dict)


def collapse_reads_with_same_umi_by_major_sequence(supporting_reads_to_counts, filtrations_count_dict):
    major_count = max(supporting_reads_to_counts.values())
    major_sequence = [read for read in supporting_reads_to_counts if supporting_reads_to_counts[read]==major_count]
    if len(major_sequence) == 1:
        # collapse to major sequence
        filtrations_count_dict['same_umi_different_read'] += sum(supporting_reads_to_counts.values()) - 1  # we still keep 1
    else:
        # tie between the majors
        filtrations_count_dict['same_umi_different_read_no_major'] += sum(supporting_reads_to_counts.values()) # discard all


def process_fastq_file(fastq_path, umi_iupac, filtrations_count_dict, ambiguous_dna_values = IUPAC.IUPACData.ambiguous_dna_values):
    umi_to_info = {}  # each umi will be mapped to a list of 3 elements -> [rec, {rec[1]:0}, rec_num]
    rec_num = 0
    umi_pattern = ''.join([f'[{ambiguous_dna_values[c]}]' for c in umi_iupac])
    for rec in FastqParser(fastq_path):  # rec is a tuple of 4 rows of a fastq file -> (header, sequence, +, quality)
        rec_num += 1
        umi = extract_umi(umi_pattern, rec[1])
        if not umi:
            filtrations_count_dict['no_umi'] += not umi # adds 1 when umi is missing (umi == '')
        else:
            if umi not in umi_to_info:
                umi_to_info[umi] = [rec, {rec[1]: 0}, rec_num]
            info = umi_to_info[umi]
            info[1][rec[1]] = info[1].get(rec[1], 0) + 1
        if rec_num % 100000 == 0:
            logger.info(f'{rec_num} reads were parsed')
    return umi_to_info, rec_num


def extract_umi(umi_pattern, read):
    m = re.search(umi_pattern, read)
    if not m:
        return ''
    return m.group()


def handle_fastq2(fastq2_path, filtered_fastq2_path, sorted_rec_nums):
    next_rec_num_in_fastq2 = 1
    result2 = ''
    fastq2_iterator = FastqParser(fastq2_path)
    for rec_num_to_save in sorted_rec_nums:
        while next_rec_num_in_fastq2 < rec_num_to_save:  # run until next relevant record
            fastq2_iterator.__next__()
            next_rec_num_in_fastq2 += 1
        result2 += '\n'.join(fastq2_iterator.__next__()) + '\n'
        next_rec_num_in_fastq2 += 1
    with open(filtered_fastq2_path, 'w') as f:
        f.write(result2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_path', help='A path to a fastq file',
                        type=lambda file_path: str(file_path) if os.path.exists(file_path) else parser.error(
                            f'{file_path} does not exist!'))
    parser.add_argument('UMIs_pattern', help='The pattern of the Unique Molecular Identifier according to IUPAC system')
    parser.add_argument('filtered_fastq_path', help='A path to a filtered fastq file',
                        type=lambda file_path: str(file_path) if os.path.exists(os.path.split(file_path)[0]) else parser.error(
                            f'{os.path.split(file_path)[0]} does not exist!'))
    parser.add_argument('--fastq2_path', help='A path to a R2 file (for paired-end reads)', required='--filtered_fastq2_path' in argv,
                        type=lambda file_path: str(file_path) if os.path.exists(file_path) else parser.error(
                            f'{file_path} does not exist!'))
    parser.add_argument('--filtered_fastq2_path', help='A path to a filtered fastq file', required='--fastq2_path' in argv,
                        type=lambda file_path: str(file_path) if os.path.exists(os.path.split(file_path)[0]) else parser.error(
                            f'{os.path.split(file_path)[0]} does not exist!'))
    # parser.add_argument('--override', help='Increase output verbosity', action='store_true')
    parser.add_argument('--log_file_path', help='A path to a file to which the filtration statistics will be written')

    args = parser.parse_args()

    handle_UMI(args.fastq_path, args.filtered_fastq_path, args.UMIs_pattern, args.fastq2_path, args.filtered_fastq2_path, args.log_file_path)