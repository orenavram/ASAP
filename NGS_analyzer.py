from time import time
start = time()
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

import sys, os
argv = sys.argv
import shutil
from email_sender import send_email
from sample_analyzer import analyze_samples
import global_params as gp
from directory_creator import create_dir
from html_editor import edit_html

logger.info('Starting '+argv[0]+'!')
logger.debug('argv = ' + str(argv))
logger.info('Usage: python3 ' + argv[0] + ' <?parameters_file_name [parameters.txt]>')
if len(argv) < 2:
    parameters_file_name = os.path.join(gp.working_dir, 'parameters.txt')
else:
    parameters_file_name = os.path.join(gp.working_dir, argv[1])
logger.info('parameter file: ' + parameters_file_name)

create_dir(gp.output_path)
gp.mixcr_output_paths = []
gp.parsed_mixcr_output_paths = []
gp.assignments_paths = []
gp.top_cdr3_clones_to_polarization_graph = 100
gp.top_cdr3_clones_to_further_analyze = 10
gp.sequence_annotation_file_suffix = '_sequence_annotations.' + gp.raw_data_file_suffix
gp.top_cdr3_annotation_file_suffix = '_top_{}_cdr3_extended_annotations.{}'.format(gp.top_cdr3_clones_to_further_analyze, gp.raw_data_file_suffix)
gp.cdr3_annotation_file_suffix = '_cdr3_annotations.' + gp.raw_data_file_suffix
gp.mutation_count_file_suffix = '_mutation_counts_distribution.' + gp.raw_data_file_suffix


#main code!
succeeded = True
server_main_url = 'http://asap.tau.ac.il/'
if os.path.exists('/Users/Oren/'):
    html_path = os.path.join(gp.working_dir, 'output.html')
    html_mode = 'w'
    smtp_server = 'mxout.tau.ac.il'
else:
    sys.path.append('/bioseq/bioSequence_scripts_and_constants')
    sys.path.append('/bioseq/asap')
    import ASAP_CONSTANTS as CONSTS
    smtp_server = CONSTS.GC.SMTP_SERVER
    html_path = os.path.join(gp.working_dir, 'output.php')
    html_mode = 'a'


#try:
analyze_samples(gp)
logger.info('Zipping results...')
out_path = os.path.join(gp.working_dir, 'outputs')
if not os.path.exists('/Users/Oren/Dropbox/Projects/wine/outputs.zip'):
    shutil.make_archive(out_path, 'zip', out_path)
else:
    logger.info('Skipping (zip already exists..)')

run_number = gp.working_dir.split('/')[-1]
#wd='/Users/Oren/Dropbox/Projects/wine/output'

logger.info('Editing html file...')
edit_html(gp, html_path, html_mode, server_main_url, run_number)

with open(html_path) as f:
    html_content = f.read()
html_content = html_content.replace(CONSTS.GC.RELOAD_TAG, '')
html_content = html_content.replace('<META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE">', '')

if succeeded:
    html_content = html_content.replace('RUNNING', 'FINISHED')
else:
    html_content = html_content.replace('RUNNING', 'FAILED')

with open(html_path, 'w') as f:
    html_content = f.write(html_content)

logger.info('Sending email...')
user_email_file = os.path.join(gp.working_dir, 'user_email.txt')
if os.path.exists(user_email_file):
    with open(user_email_file) as f:
        addressee = f.read().rstrip()

    server_url = 'http://asap.tau.ac.il/results'
    email_content = '''Hello,

The results for your ASAP run are ready at:
{}


Please note: the results will be kept on the server for three months.

Thanks
ASAP Team
    '''.format(os.path.join(server_url, run_number, 'output.php'))

    send_email(smtp_server=smtp_server, sender='TAU BioSequence <bioSequence@tauex.tau.ac.il>',
               receiver=addressee, subject='Your ASAP run is ready!', content=email_content)

logger.info(argv[0] + ' is DONE!!')


end = time()
hours = int((end - start) / 3600)
minutes = int(((end - start) % 3600) / 60)
seconds = int((end - start) % 60)
logger.info('Finished joining samples. Took {}:{}:{} hours.'.format(hours, minutes, seconds))


print('Bye.')