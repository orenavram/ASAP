from text_handler import get_next_relevant_rstriped_line
from sys import argv
import os
import logging
logger = logging.getLogger('main')

if len(argv)<3:
    params_file_path = '/Users/Oren/Dropbox/Projects/wine/parameters_242.txt'
    debug_run = ''
    # debug_run = False
else:
    params_file_path = argv[1]
    debug_run = argv[2]

logger.info('DEBUG MODE?', debug_run if debug_run else False)

with open(params_file_path) as f:

    line = get_next_relevant_rstriped_line(f)
    #String that represents the path to MiXCR executable file
    working_dir = line.rstrip('/')

    line = get_next_relevant_rstriped_line(f)
    #String that represents the path to MiXCR executable file
    path_to_mixcr = line

    #String that represent a sample number. Should match the FASTQ file names
    line = get_next_relevant_rstriped_line(f)
    sample = line

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

    #Boolean that indicates whether the samples originated in mice (and not human)
    line = get_next_relevant_rstriped_line(f)
    MMU = True if line == 'True' else False

        #String that represent the raw data files suffix. txt / xls / etc...
    line = get_next_relevant_rstriped_line(f)
    raw_data_file_suffix = line

output_path = os.path.join(working_dir, 'outputs')

run_output_paths = []

