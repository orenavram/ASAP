#!/shared/python/anaconda3.5/bin/python

import os, sys
# if os.path.exists('/Users/Oren/'): #local run
#     sys.path.append('./auxiliaries') # this is where both ASAP_CONSTANTS and GENERAL_CONSTANTS are located in my comp
# else: #run on host-ibis
#     sys.path.append('/bioseq/bioSequence_scripts_and_constants') # this is where GENERAL_CONSTANTS is located in host-ibis3
#     sys.path.append('/bioseq/asap/ASAP/auxiliaries') # this is where ASAP_CONSTANTS is located in host-ibis3

#import GENERAL_CONSTANTS as GC


#!/shared/python/anaconda3.5/bin/python

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

# general paths
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'

# external programs
# IMPORTANT: one must run the command: setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}" ahead of this mafft command so all components will be found...
MAFFT_v7_222 = '/bioseq/Programs/MAFFT_7.222/installation/bin/mafft' # v7.222
MiXCR_dir = '/bioseq/asap/MiXCR_2.0.4'

ASAP_URL = 'http://asap.tau.ac.il'
#ASAP_LOG = '/bioseq/asap/ASAP_runs.log'

RELOAD_INTERVAL = 30
RELOAD_TAGS = """
    <META HTTP-EQUIV="REFRESH" CONTENT="{}"/>
    <META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE"/>
""".format(str(RELOAD_INTERVAL))


ASAP_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'ASAP')
ASAP_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'ASAP')
ASAP_RESULTS_URL = os.path.join(ASAP_URL, 'results')
ASAP_HTML_DIR = '/data/www/html/asap'
#ASAP_EXEC = '/bioseq/asap'

MAIN_SCRIPT = '/bioseq/asap/ASAP/pipeline/NGS_analyzer.py'

#path to example runs
EXAMPLE_FILE_RUN1_R1 = os.path.join(ASAP_HTML_DIR, 'example/run1/242_R1.fastq')
EXAMPLE_FILE_RUN1_R2 = os.path.join(ASAP_HTML_DIR, 'example/run1/242_R2.fastq')
EXAMPLE_FILE_RUN2_R1 = os.path.join(ASAP_HTML_DIR, 'example/run2/242_R1.fastq')
EXAMPLE_FILE_RUN2_R2 = os.path.join(ASAP_HTML_DIR, 'example/run2/242_R2.fastq')