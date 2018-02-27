from time import time
start = time()
from sys import argv
import os
import shutil
from email_sender import send_email
from sample_analyzer import analyze_samples
import global_params as gp
from directory_creator import create_dir
from text_handler import logger
from html_editor import edit_html
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

if len(argv) < 2:
    logger.error('Usage: python3 '+argv[0]+' <working_directory>; <?parameters_file_name [parameters.txt]>')
    exit(-1)

logger.info('Starting '+argv[0]+'!')
logger.debug('argv = ' + str(argv))

create_dir(gp.output_path)
gp.mixcr_output_paths = []
gp.parsed_mixcr_output_paths = []
gp.assignments_paths = []
gp.sequence_annotation_file_suffix = '_sequence_annotations.' + gp.raw_data_file_suffix
gp.cdr3_annotation_file_suffix = '_cdr3_annotations.' + gp.raw_data_file_suffix
gp.mutation_count_file_suffix = '_mutation_counts_distribution.' + gp.raw_data_file_suffix
gp.top_cdr3_clones = 10

#main code!
succeeded = True
server_main_url = 'http://asap.tau.ac.il/'
if os.path.exists('/Users/Oren/'):
    html_path = os.path.join(gp.working_dir, 'output.html')
    html_mode = 'w'
else:
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
'''
except Exception as e:
    succeeded = False
    print('$'*100 + '\nanalyze_samples crashed!\n' + '$'*100)
    logger.error(e.args)
'''

with open(html_path) as f:
    html_content = f.read()
html_content = html_content.replace('<META HTTP-EQUIV="REFRESH" CONTENT=30> </HEAD>', '')
html_content = html_content.replace('<META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE"> </HEAD>', '')

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
    ASAP Team'''.format(os.path.join(server_url, run_number, 'output.php'))

    send_email(sender='TAU BioSequence <bioSequence@tauex.tau.ac.il>', receiver=addressee,
               subject='Your ASAP run is ready!', content=email_content)

logger.info(argv[0] + ' is DONE!!')

'''
f.write('<li><a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/2_consensus.fasta" target="_blank">2nd clone cluster</a></li>\n')
f.write('<li><a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/3_consensus.fasta" target="_blank">3rd clone cluster</a></li>\n')
f.write('<li><a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/4_consensus.fasta" target="_blank">4th clone cluster</a></li>\n')
f.write('<li><a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/5_consensus.fasta" target="_blank">5th clone cluster</a></li>\n')
'''

end = time()
hours = int((end - start) / 3600)
minutes = int(((end - start) % 3600) / 60)
seconds = int((end - start) % 60)
logger.info('Finished joining samples. Took {}:{}:{} hours.'.format(hours, minutes, seconds))


print('Bye.')