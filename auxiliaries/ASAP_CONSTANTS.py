#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

# general paths
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'

# external programs
# IMPORTANT: one must run the command: setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}" ahead of this mafft command so all components will be found...
MAFFT_v7_222 = '/bioseq/Programs/MAFFT_7.222/installation/bin/mafft' # v7.222

ASAP_URL = 'https://asap.tau.ac.il'
#ASAP_LOG = '/bioseq/asap/ASAP_runs.log'

RELOAD_INTERVAL = 30
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'


ASAP_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'ASAP')
ASAP_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'ASAP')
ASAP_RESULTS_URL = os.path.join(ASAP_URL, 'results')
ASAP_HTML_DIR = '/data/www/html/asap'
ASAP_EXEC = '/bioseq/asap/ASAP'

MAIN_SCRIPT = os.path.join(ASAP_EXEC, 'pipeline/NGS_analyzer.py')

#path to example runs
EXAMPLE_FILE_RUN1_R1 = os.path.join(ASAP_HTML_DIR, 'example', 'run1', '242_R1.fastq')
EXAMPLE_FILE_RUN1_R2 = os.path.join(ASAP_HTML_DIR, 'example', 'run1', '242_R2.fastq')
EXAMPLE_FILE_RUN2_R1 = os.path.join(ASAP_HTML_DIR, 'example', 'run2', '242_R1.fastq')
EXAMPLE_FILE_RUN2_R2 = os.path.join(ASAP_HTML_DIR, 'example', 'run2', '242_R2.fastq')

#path to mass_spec_db
MASS_SPEC_DB_MOUSE = os.path.join(ASAP_HTML_DIR, 'MassSpecDB', 'murine.fasta')
MASS_SPEC_DB_HUMAN = os.path.join(ASAP_HTML_DIR, 'MassSpecDB', 'human.fasta')

#path to IMGT reference library
IMGT_LIB = os.path.join(ASAP_EXEC, 'imgt.201822-5.sv4.json')