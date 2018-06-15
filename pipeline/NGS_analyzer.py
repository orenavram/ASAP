try:
    import matplotlib  # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay https://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
    matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay

    import sys, os, traceback, shutil

    if os.path.exists('/bioseq/'): #remote run
        #sys.path.append('/bioseq/bioSequence_scripts_and_constants') # this is where GENERAL_CONSTANTS is located in host-ibis3
        sys.path.append('/bioseq/asap/ASAP/auxiliaries') # this is where ASAP_CONSTANTS is located in host-ibis3
    else: #run on host-ibis
        sys.path.append('./auxiliaries')  # this is where ASAP_CONSTANTS (and GENERAL_CONSTANTS, currently unused) is located in my comp

    from time import time, ctime, sleep
    import logging

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    logger.info(f'sys.path: {sys.path}')
    import global_params as gp

    start = time()

    argv = sys.argv
    from email_sender import send_email
    from directory_creator import create_dir
    from html_editor import edit_success_html, edit_failure_html

    logger.info('Starting ' + argv[0] + '!')
    logger.debug('argv = ' + str(argv))
    if len(argv) < 2:
        logger.info('Usage: python ' + argv[0] + ' <parameters_file_path>')
        exit()
    else:
        parameters_file_path = argv[1]
    logger.info('parameters file: ' + parameters_file_path)

    create_dir(gp.output_path)

    # main code!
    output_path = os.path.join(gp.working_dir, 'output.html')

    import ASAP_CONSTANTS as CONSTS

    gp.initial_db_path = CONSTS.MASS_SPEC_DB_MOUSE if gp.MMU else CONSTS.MASS_SPEC_DB_HUMAN

    from sample_analyzer import analyze_samples
    analyze_samples(gp)
    succeeded = True
except Exception as e:
    error_msg = 'ASAP calculation crashed :('
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logger.error('\n\n' + '$'*60 + '\n\n')
    logger.error(ctime() + ': ' + error_msg + '\n\n')
    logger.error(str(fname) +': ' + str(exc_type) + ', at line: ' + str(exc_tb.tb_lineno) + '\n\n')
    logger.error('$'*60)
    succeeded = False

run_number = gp.working_dir.split('/')[-1]
if succeeded:
    logger.info('Zipping results...')
    if not os.path.exists(os.path.join(gp.output_path + '.zip')):
        shutil.make_archive(gp.output_path, 'zip', gp.output_path)
    else:
        logger.info('Skipping (zip already exists..)')
    logger.info('Editing html file...')
    edit_success_html(gp, output_path, CONSTS.ASAP_URL, run_number)
else:
    edit_failure_html(output_path, error_msg)

#Change running status
# with open(output_path) as f:
#     html_content = f.read()
# if succeeded:
#     html_content = html_content.replace('RUNNING', 'FINISHED')
# else:
#     html_content = html_content.replace('RUNNING', 'FAILED')
# with open(output_path, 'w') as f:
#     html_content = f.write(html_content)

output_url = os.path.join(CONSTS.ASAP_RESULTS_URL, run_number, 'output.html')

logger.info('Sending email...')
user_email_file = os.path.join(gp.working_dir, 'user_email.txt')
if os.path.exists(user_email_file):
    with open(user_email_file) as f:
        addressee = f.read().rstrip()

    results_page = output_url
    if succeeded:
        email_content = '''Hello,
    
The results for your ASAP run are ready at:
{}

Please note: the results will be kept on the server for three months.
        '''.format(results_page)
    else:
        email_content = 'ASAP calculation failed. For further information please visit: {}'.format(results_page)

    email_content += '\n\nThanks, ASAP Team'

    send_email(smtp_server=CONSTS.SMTP_SERVER, sender=CONSTS.ADMIN_EMAIL,
               receiver=addressee, subject='Your ASAP run is ready!', content=email_content)

logger.info(argv[0] + ' is DONE!!')


end = time()
hours = int((end - start) / 3600)
minutes = int(((end - start) % 3600) / 60)
seconds = int((end - start) % 60)
logger.info('Finished joining samples. Took {}:{}:{} hours.'.format(hours, minutes, seconds))
with open(os.path.join(gp.output_path, 'time.txt'), 'w') as f:
    f.write('{}:{}:{}'.format(hours, minutes, seconds))

logger.info(f'Waiting {2*CONSTS.RELOAD_INTERVAL} seconds to remove html refreshing headers...')
# Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
sleep(2*CONSTS.RELOAD_INTERVAL)
with open(output_path) as f:
    html_content = f.read()
html_content = html_content.replace(CONSTS.RELOAD_TAGS, '')
with open(output_path, 'w') as f:
    html_content = f.write(html_content)

print('Done. Bye.')
