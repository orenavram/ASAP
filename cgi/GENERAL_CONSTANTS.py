#!/shared/python/anaconda3.5/bin/python

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
ADMIN_USER_NAME = 'bioSequence'
ADMIN_PASSWORD = 'elana'
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
RELOAD_TAG = '<META HTTP-EQUIV="REFRESH" CONTENT="{}"/>'.format(str(RELOAD_INTERVAL))

# # copy-pasted from email_sender.py
# def send_email(sender, receiver, subject='', content=''):
#     from email.mime.text import MIMEText
#     from smtplib import SMTP
#     msg = MIMEText(content)
#     msg['Subject'] = subject
#     msg['From'] = sender
#     msg['To'] = receiver
#     s = SMTP(SMTP_SERVER)
#     s.send_message(msg)
#     s.quit()
#
#
# import logging, os
# # copy-pasted from email_sender.py
# def create_dir(path, logger = logging.getLogger('main')):
#     if not os.path.exists(path):
#         logger.info("Creating directory:\n" + path)
#         os.makedirs(path)