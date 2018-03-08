#!/shared/python/anaconda3.5/bin/python

import GENERAL_CONSTANTS as GC
import os

ASAP_RESULTS_DIR = os.path.join(GC.SERVERS_RESULTS_DIR, 'ASAP')
ASAP_LOGS_DIR = os.path.join(GC.SERVERS_LOGS_DIR, 'ASAP')
ASAP_RESULTS_URL = os.path.join(GC.ASAP_URL, 'results')
ASAP_HTML_DIR = '/data/www/html/asap'
ASAP_EXEC = '/bioseq/asap'

MAIN_SCRIPT = '/bioseq/asap/ASAP/NGS_analyzer.py'

#path to example runs
EXAMPLE_FILE_RUN1_R1 = os.path.join(ASAP_HTML_DIR, 'example/run1/242_R1.fastq')
EXAMPLE_FILE_RUN1_R2 = os.path.join(ASAP_HTML_DIR, 'example/run1/242_R2.fastq')
EXAMPLE_FILE_RUN2_R1 = os.path.join(ASAP_HTML_DIR, 'example/run2/242_R1.fastq')
EXAMPLE_FILE_RUN2_R2 = os.path.join(ASAP_HTML_DIR, 'example/run2/242_R2.fastq')