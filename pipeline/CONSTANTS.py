#!/groups/pupko/modules/python-anaconda3.6.5/bin/python

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

#ASAP_LOG = '/bioseq/asap/ASAP_runs.log'

RELOAD_INTERVAL = 10
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

WEBSERVER_NAME = 'asap'

ASAP_URL = f'https://{WEBSERVER_NAME}.tau.ac.il'

ASAP_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, WEBSERVER_NAME)
ASAP_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, WEBSERVER_NAME)
ASAP_RESULTS_URL = os.path.join(ASAP_URL, 'results')
ASAP_HTML_DIR = f'/data/www/html/{WEBSERVER_NAME}'
ASAP_EXEC = f'/bioseq/{WEBSERVER_NAME}'

SUBMISSIONS_LOG = f'/bioseq/{WEBSERVER_NAME}/submissions_log.txt'

MAIN_SCRIPT = os.path.join(ASAP_EXEC, 'pipeline/NGS_analyzer.py')

RESULT_WEBPAGE_NAME = 'result.html'

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