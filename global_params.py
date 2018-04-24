from text_handler import get_next_relevant_rstriped_line
from sys import argv
import os
import logging
logger = logging.getLogger('main')

if len(argv)<2:
    params_file_path = '/Users/Oren/Dropbox/Projects/wine/parameters.txt'
else:
    params_file_path = argv[1]

logger.info('Parsing parameter file: {}'.format(params_file_path))

with open(params_file_path) as f:

    line = get_next_relevant_rstriped_line(f)
    #String that represents the path to a folder where the output dir will be generated
    working_dir = line.rstrip('/')

    line = get_next_relevant_rstriped_line(f)
    #String that represents the path to MiXCR executable file
    path_to_mixcr = line

    #Integer that represents number of runs (the output will be +1 because of the joint)
    line = get_next_relevant_rstriped_line(f)
    number_of_runs = int(line)

    #List of strings that represent the chains
    line = get_next_relevant_rstriped_line(f)
    chains = line.split(',')

    #Integer that represents minimal threshold of reads' length (nucleotides)
    line = get_next_relevant_rstriped_line(f)
    len_threshold = int(line)

    #Integer that represents minimal threshold for reads' average quality
    #(reads with average quality lower than that are filtered out)
    line = get_next_relevant_rstriped_line(f)
    sequencing_quality_threshold = int(line)

    #Integer that represents the k-top clones to be further analyzed
    line = get_next_relevant_rstriped_line(f)
    top_cdr3_clones_to_further_analyze = int(line)

    #String that indicates whether the samples originated in mice (and not human)
    line = get_next_relevant_rstriped_line(f)
    MMU = True if line == 'Mouse' else False

    #String that represent the raw data files suffix. txt / xls / etc...
    line = get_next_relevant_rstriped_line(f)
    raw_data_file_suffix = line

    #String (yes/no) that represent whether to add mass_spec_seq to each aa sequence in the fasta output
    line = get_next_relevant_rstriped_line(f)
    add_mass_spec_seq = True if line == 'yes' else False

output_path = os.path.join(working_dir, 'outputs')

run_output_paths = []

mixcr_output_paths = []
parsed_mixcr_output_paths = []
assignments_paths = []
cdr3_analysis_paths = []
top_cdr3_clones_to_clonal_expansion_graph = 100
top_cdr3_clones_to_further_analyze = min(top_cdr3_clones_to_further_analyze, top_cdr3_clones_to_clonal_expansion_graph)
sequence_annotation_file_suffix = '_aa_sequence_annotations.' + raw_data_file_suffix
top_cdr3_annotation_file_suffix = '_top_{}_cdr3_extended_annotations.{}'.format(top_cdr3_clones_to_further_analyze, raw_data_file_suffix)
cdr3_annotation_file_suffix = '_cdr3_annotations.' + raw_data_file_suffix
mutations_file_suffix = '_mutations.' + raw_data_file_suffix
#Ka_Ks_file_suffix = '_Ka_Ks_analysis_2.' + raw_data_file_suffix
