try:
    from time import time

    start = time()
    import logging

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    import sys, os, traceback

    argv = sys.argv
    import shutil
    from email_sender import send_email
    from directory_creator import create_dir
    from html_editor import edit_success_html, edit_failure_html

    import global_params as gp

    logger.info('Starting ' + argv[0] + '!')
    logger.debug('argv = ' + str(argv))
    logger.info('Usage: python3 ' + argv[0] + ' <?parameters_file_name [parameters.txt]>')
    if len(argv) < 2:
        parameters_file_name = os.path.join(gp.working_dir, 'parameters.txt')
    else:
        parameters_file_name = os.path.join(gp.working_dir, argv[1])
    logger.info('parameter file: ' + parameters_file_name)

    create_dir(gp.output_path)

    # main code!
    output_path = os.path.join(gp.working_dir, 'output.html')
    if os.path.exists('/Users/Oren/'): #local run
        sys.path.append('./cgi')  # this is where GENERAL_CONSTANTS.py is located in my comp
        html_mode = 'w'
    else: #run on host-ibis
        sys.path.append('/bioseq/bioSequence_scripts_and_constants')  # this is where GENERAL_CONSTANTS.py is located in host-ibis3
        sys.path.append('/bioseq/asap')  # this is where ASAP_CONSTANTS is located in host-ibis3
        html_mode = 'a'

    import ASAP_CONSTANTS as CONSTS

    from sample_analyzer import analyze_samples
    analyze_samples(gp)
    succeeded = True
except Exception as e:
    error_msg = 'ASAP calculation crashed :('
    logger.error(error_msg)
    succeeded = False

run_number = gp.working_dir.split('/')[-1]
if succeeded:
    logger.info('Zipping results...')
    if not os.path.exists('/Users/Oren/Dropbox/Projects/wine/outputs.zip'):
        shutil.make_archive(gp.output_path, 'zip', gp.output_path)
    else:
        logger.info('Skipping (zip already exists..)')
    logger.info('Editing html file...')
    edit_success_html(gp, output_path, html_mode, CONSTS.GC.ASAP_URL, run_number)
else:
    edit_failure_html(output_path, html_mode, error_msg)

#Change running status
with open(output_path) as f:
    html_content = f.read()
html_content = html_content.replace(CONSTS.GC.RELOAD_TAGS, '')
if succeeded:
    html_content = html_content.replace('RUNNING', 'FINISHED')
else:
    html_content = html_content.replace('RUNNING', 'FAILED')
with open(output_path, 'w') as f:
    html_content = f.write(html_content)

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

    send_email(smtp_server=CONSTS.GC.SMTP_SERVER, sender='TAU BioSequence <bioSequence@tauex.tau.ac.il>',
               receiver=addressee, subject='Your ASAP run is ready!', content=email_content)

logger.info(argv[0] + ' is DONE!!')


end = time()
hours = int((end - start) / 3600)
minutes = int(((end - start) % 3600) / 60)
seconds = int((end - start) % 60)
logger.info('Finished joining samples. Took {}:{}:{} hours.'.format(hours, minutes, seconds))
with open(os.path.join(gp.output_path, 'time.txt'), 'w') as f:
    f.write('{}:{}:{}'.format(hours, minutes, seconds))
print('Bye.')