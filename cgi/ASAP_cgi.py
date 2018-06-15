#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

# old shebang: #!/shared/python/anaconda3.5/bin/python

import os
import sys
import sh
import cgi
import cgitb
import subprocess
from time import time, ctime
from random import randint

#sys.path.append('/bioseq/bioSequence_scripts_and_constants')
sys.path.append('/bioseq/asap/ASAP/auxiliaries')
import ASAP_CONSTANTS as CONSTS
from directory_creator import create_dir
from email_sender import send_email

def print_hello_world(output_path = '', run_number = 'NO_RUN_NUMBER'):

    hello_world_html = """
<html>
    <body>
        <h2>Hello World! """ + run_number + """</h2>
    </body>
</html>
    """
    if not output_path:
        print(hello_world_html)
    else:
        with open(output_path, 'w') as f:
            f.write(hello_world_html)


def write_html_prefix(output_path, run_number):
    # ESCAPE CURLY BRACES BY DOUBLING THEM!! {{ or }}
    with open(output_path, 'w') as f:
        f.write(f'''<html><head>
        
    <meta http-equiv="cache-control" content="no-cache, must-revalidate, post-check=0, pre-check=0" />
    <meta http-equiv="cache-control" content="max-age=0" />
    <meta http-equiv="expires" content="0" />
    <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
    <meta http-equiv="pragma" content="no-cache" />
    {CONSTS.RELOAD_TAGS}
    
    <title>ASAP Job {run_number}</title>
    <link rel="shortcut icon" type="image/x-icon" href="{CONSTS.ASAP_URL}/pics/ASAP_logo.gif" />
    
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
    <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">
    
    <link rel="stylesheet" href="{CONSTS.ASAP_URL}/css/general.css">
    
    </head><body>
    <nav role="navigation" class="navbar navbar-fixed-top">
        <div class="jumbotron" id="jumbo">
            <div class="container">
                <div class="row" id="title-row">
                    <div class="col-md-1">
                    </div>
                    <div class="col-md-1">
                        <img src="{CONSTS.ASAP_URL}/pics/ASAP_logo.gif" id="antibody_image" class="img-rounded">
                    </div>
                    <div class="col-md-9">
                        &nbsp;&nbsp;&nbsp;&nbsp;<span id="asap-title">ASAP</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title"><b>A</b> web server for Ig-<b>S</b>eq <b>A</b>nalysis <b>P</b>ipeline</span><br>
                    </div>
                </div>
            </div>       
        </div>
    </nav>
    <div id="behind-nav-bar-results">
    </div>
''')


def write_info_paragraph_to_html(output_path):
    with open(output_path, 'a') as f:
        f.write(f"""<br><div class="container" style="font-size: 20px;" align="justify"> 
<H1 align=center>Job Status - <FONT color='red'>RUNNING</FONT></h1>
<br>ASAP is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for three months.
</div>
<br>""")


def write_running_parameters_to_html(output_path, job_title, number_of_duplicates, urls_to_reads_files, files_names, MMU, len_threshold, qlty_threshold, number_of_clones_to_analyze, raw_data_suffix, mass_spec_seq):

    with open(output_path, 'a') as f:

        #regular params row
        f.write("""<div class="container"><u><h3>Running Parameters:</h3></u><br>""")

        f.write('<div class="row" style="font-size: 20px;">')
        if job_title != '':
            f.write('<div class="col-md-6">')
            f.write(f'<b>Job title: </b>{job_title}<br><br>')
            f.write('</div>')
            f.write('</div><div class="row" style="font-size: 20px;">')

        f.write('<div class="col-md-3">')
        f.write(f'<b>Organism: </b>{MMU}')
        f.write('</div>')

        # f.write('<div class="col-md-3">')
        # f.write('<b>Chains: </b>{}<br>'.format(chains))
        # f.write('</div>')

        f.write('</div><br>')

        f.write('<div class="row" style="font-size: 20px;">')
        for i in range(number_of_duplicates):
            # show only path suffix... http://asap.tau.ac.il/results/1520442230/reads/run1/242_R1.fastq
            f.write('<div class="col-md-4">')
            f.write('<b>Run {} sequence data:</b><br>'.format(i + 1))
            f.write(f'R1 = <a href="{urls_to_reads_files[i * 2]}" target=_blank>{files_names[i*2]}</A><br>')
            f.write(f'R2 = <a href="{urls_to_reads_files[i * 2 + 1]}" target=_blank>{files_names[i*2+1]}</A><br>')
            f.write('</div>')

        f.write('</div>')

        #Advanced params row
        f.write('<br><u><h3>Advanced Parameters:</h3></u><br>')

        f.write('<div class="row">')
        f.write('<div class="col-md-2">')
        f.write(f'<b>Min read length: </b>{len_threshold}<br>')
        f.write('</div>')

        f.write('<div class="col-md-2">')
        f.write(f'<b>Min read quality: </b>{qlty_threshold}<br>')
        f.write('</div>')

        f.write('<div class="col-md-3">')
        f.write(f'<b>MassSpec sequence: </b>{mass_spec_seq}<br>')
        f.write('</div>')

        f.write('<div class="col-md-3">')
        f.write(f'<b>#clones to analyze: </b>{number_of_clones_to_analyze}<br>')
        f.write('</div>')

        f.write('<div class="col-md-2">')
        f.write(f'<b>Raw data format: </b>{raw_data_suffix}<br>')
        f.write('</div>')

        f.write('</div></div>')


def write_pair_file(debug_path, pair, run_content, run_filename, run_dir):
    if run_filename.endswith('gz'):
        #TODO: This option is not working properly. Files are being uploaded but they're non sense :(
        # open_operator = gzip.open
        local_file_name = f'R{pair}.fastq.gz'
    else:
        local_file_name = f'R{pair}.fastq'
    open_operator = open
    with open(debug_path, 'a') as f:
        f.write(f'{local_file_name} is being handled with {open_operator}\n')
    local_file_path = os.path.join(run_dir, local_file_name)
    with open_operator(local_file_path, 'wb') as f:
        f.write(run_content)

    #avoid double zipping:
    if local_file_path.endswith('gz'):
        sh.gunzip(local_file_path)

    with open(debug_path, 'a') as f:
        f.write(f'R{pair} was handled successfully\n')


def process_uploaded_files(run, debug_path):
    #for debugging
    with open(debug_path, 'a') as f:
        f.write(f'uploading file of {run}\n')
    run_R1_filename = form[run + '_R1'].filename
    run_R2_filename = form[run + '_R2'].filename
    # for debugging
    with open(debug_path, 'a') as f:
        f.write(f'file names are:\n{run_R1_filename} (of type {type(form[run + "_R1"].value)}) and {run_R2_filename} (of type {type(form[run + "_R2"].value)})\n')
    run_R1_content = form[run + '_R1'].value
    run_R2_content = form[run + '_R2'].value
    # for debugging
    with open(debug_path, 'a') as f:
        f.write(f'{run_R1_filename} first 100 chars are: {run_R1_content[:100]}\n{run_R2_filename} first 100 chars are: {run_R2_content[:100]}\n')
    run_dir = os.path.join(wd, 'reads', run)
    create_dir(run_dir)

    write_pair_file(debug_path, '1', run_R1_content, run_R1_filename, run_dir)

    write_pair_file(debug_path, '2', run_R2_content, run_R2_filename, run_dir)

    return run_R1_filename, run_R2_filename


def prepare_parameters_file(parameters_file, parameters):
        with open(parameters_file, 'w') as f:
            f.write("""
    #String that represents the path to a folder where the output dir will be generated
    {}

    #String that represents the path to MiXCR executable file
    {}

    #Integer that represents number of runs (the output will be +1 because of the joint)
    {}

    #List of strings that represent the chains
    # comma delimited (if more than one) without spaces!
    {}

    #Integer that represents minimal threshold of reads' length (nucleotides)
    {}

    #Integer that represents minimal threshold for reads' average quality
    #(reads with average quality lower than that are filtered out)
    {}

    #Integer that represents the k-top clones to be further analyzed
    {}

    #String that indicates whether the samples originated in mice (and not human)
    {}

    #String that represent the raw data files suffix. txt / xls / etc...
    {}

    #String that represent the mass_spec_seq to add for the proteomics DB
    {}

    """.format(*parameters))


def write_cmds_file(cmds_file, run_number, parameters_file):
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
    new_line_delimiter = ';!@#'
    with open(cmds_file, 'w') as f:
        f.write('setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}"')
        f.write(new_line_delimiter)
        f.write('module load python/anaconda_python-3.6.4')
        f.write(new_line_delimiter)
        f.write(' '.join(['python', CONSTS.MAIN_SCRIPT, parameters_file, '""']))
        f.write('\t' + 'ASAP_' + run_number)


# prints detailed error report on BROWSER when cgi crashes
# This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
cgitb.enable()

#print_hello_world() # for debugging
form = cgi.FieldStorage() # extract POSTed object

#random_chars = "".join(choice(string.ascii_letters + string.digits) for x in range(20))
run_number = str(round(time())) + str(randint(10**19,10**20-1)) # adding 20 random digits to prevent users see data that are not their's
if False:
    run_number = 'debug'#str(round(time())) + str(randint(1000,9999)) # adding 4 random figures to prevent users see data that are not their's

results_url = os.path.join(CONSTS.ASAP_RESULTS_URL, run_number)
output_url = os.path.join(results_url, 'output.html')

wd = os.path.join(CONSTS.ASAP_RESULTS_DIR, run_number)
create_dir(wd)
output_path = os.path.join(wd, 'output.html')
cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')

#print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
#print_hello_world() # comment out for debugging
#print_hello_world(output_path, run_number) # comment out for debugging

write_html_prefix(output_path, run_number) # html's prefix must be written BEFORE redirecting...

print('Location: ' + output_url) #Redirects to the results url. MUST appear before any other print.
print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
sys.stdout.flush() #must be flushed immediately!!!

# Send me a notification email every time there's a new request
send_email(smtp_server=CONSTS.SMTP_SERVER, sender=CONSTS.ADMIN_EMAIL,
           receiver='orenavram@gmail.com', subject=f'ASAP - A new job has been submitted: {run_number}', content=os.path.join(CONSTS.ASAP_URL,'results',run_number))

try:
    write_info_paragraph_to_html(output_path)

    with open(cgi_debug_path, 'a') as f:
        # form debugging
        f.write(f'{"#"*50}\n{ctime()}: A new CGI request has been recieved!\n')
        f.write(f'These are the keys that the CGI received:\n{"; ".join(sorted(form.keys()))}\n\n')
        f.write('Form values are:\n')
        for key in sorted(form.keys()):
            if 'run' not in key:
                f.write(f'{key} = {form[key]}\n')
        # for key in sorted(form.keys()):
        #     if 'run' in key:
        #         f.write('100 first characters of {} = '.format(key))
        #         f.write('{}\n'.format(form[key].value[:100]))
        f.write('\n\n')

    #extract form's values:
    user_email = form['email'].value.strip()
    example_page = form['example_page'].value
    chains = ','.join(form.getlist('chains'))
    MMU = form['MMU'].value.lower()
    len_threshold = form['len_threshold'].value
    qlty_threshold = form['qlty_threshold'].value
    number_of_clones_to_analyze = form['number_of_clones_to_analyze'].value
    raw_data_suffix = form['raw_data_suffix'].value
    mass_spec_seq = form['mass_spec_seq'].value.strip()
    # #if this option is unchecked, it won't be send in the json (i.e., form['add_mass_spec_seq'].value might not work...)
    # if 'add_mass_spec_seq' in form:
    #     add_mass_spec_seq = 'yes'
    job_title = ''
    if form['job_title'].value != '':
        job_title = form['job_title'].value.strip()

    files_names = []
    if example_page == 'no':
        # handling uploaded files:
        number_of_duplicates = 1

        #at least one run should exist:
        run1_R1_filename, run1_R2_filename = process_uploaded_files('run1', cgi_debug_path)
        files_names.extend([run1_R1_filename, run1_R2_filename])
        # additional files might not exist
        run2_R1_filename, run2_R2_filename, run3_R1_filename, run3_R2_filename, = '', '', '', ''
        if form['run2_R1'].filename != '': #handle run2 if any
            run2_R1_filename, run2_R2_filename = process_uploaded_files('run2', cgi_debug_path)
            number_of_duplicates += 1
            files_names.extend([run2_R1_filename, run2_R2_filename])
        if form['run3_R1'].filename != '': #handle run3 if any
            run3_R1_filename, run3_R2_filename = process_uploaded_files('run3', cgi_debug_path)
            number_of_duplicates += 1
            files_names.extend([run3_R1_filename, run3_R2_filename])
    else:
        # no files to handle, just copy them to the example folder (if needed).
        number_of_duplicates = 2

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write(f'Number of duplicates is {number_of_duplicates}\n')

    paths_to_reads_files = []
    urls_to_reads_files = []
    for i in range(number_of_duplicates):
        run = 'run' + str(i + 1)
        path_to_reads = os.path.join(wd, 'reads', run)
        url_to_reads = os.path.join(results_url, 'reads', run)
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write(f'{ctime()}: Creating path for reads...\n')
        create_dir(path_to_reads)
        paths_to_reads_files.append(os.path.join(path_to_reads, 'R1.fastq'))
        paths_to_reads_files.append(os.path.join(path_to_reads, 'R2.fastq'))
        urls_to_reads_files.append(os.path.join(url_to_reads, 'R1.fastq'))
        urls_to_reads_files.append(os.path.join(url_to_reads, 'R2.fastq'))

    if example_page == 'yes':
        files_names = ['rep1_R1.fastq', 'rep1_R2.fastq', 'rep2_R1.fastq', 'rep2_R2.fastq']
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write(f'{ctime()}: Copying example files...\n')
        # copy example data
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write(f'Fetching: rsync -avz {CONSTS.EXAMPLE_FILE_RUN1_R1} {paths_to_reads_files[0]}\n')
        os.system(f'rsync -avz {CONSTS.EXAMPLE_FILE_RUN1_R1} {paths_to_reads_files[0]}')
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write(f'{ctime()}: ls of {os.path.join(wd, "reads", "run1")} yields:\n{os.listdir(path_to_reads)}\n')
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write(f'Fetching: rsync -avz {CONSTS.EXAMPLE_FILE_RUN1_R2} {paths_to_reads_files[1]}\n')
        os.system(f'rsync -avz {CONSTS.EXAMPLE_FILE_RUN1_R2} {paths_to_reads_files[1]}')
        os.system(f'rsync -avz {CONSTS.EXAMPLE_FILE_RUN2_R1} {paths_to_reads_files[2]}')
        os.system(f'rsync -avz {CONSTS.EXAMPLE_FILE_RUN2_R2} {paths_to_reads_files[3]}')
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('{}: Files were copied successfully.\n'.format(ctime()))

    # copy mass_spec db to the results folder
    if MMU == 'mouse':
        os.system(f'rsync -avz {CONSTS.MASS_SPEC_DB_MOUSE} {wd}')
    else:
        os.system(f'rsync -avz {CONSTS.MASS_SPEC_DB_HUMAN} {wd}')

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write(f'{ctime()}: ls of {path_to_reads} yields:\n{os.listdir(path_to_reads)}\n')

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write(f'{ctime()}: write_running_parameters_to_html...\n')

    write_running_parameters_to_html(output_path, job_title, number_of_duplicates, urls_to_reads_files, files_names, MMU, len_threshold, qlty_threshold, number_of_clones_to_analyze, raw_data_suffix, mass_spec_seq)

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write(f'{ctime()}: Running parameters were written to html successfully.\n')

    # This is hidden field that only spammer bots might fill in...
    confirm_email_add = form['confirm_email'].value  # if it is contain a value it is a spammer.

    parameters_file = os.path.join(wd, 'parameters.txt')
    parameters = [wd, CONSTS.MiXCR_dir, number_of_duplicates, chains, len_threshold, qlty_threshold, number_of_clones_to_analyze, MMU, raw_data_suffix, mass_spec_seq]
    prepare_parameters_file(parameters_file, parameters)

    cmds_file = os.path.join(wd, 'qsub.cmds')
    write_cmds_file(cmds_file, run_number, parameters_file)

    log_file = cmds_file.replace('cmds', 'log')
    # complex command with more than one operation (module load + python q_submitter.py)
    #submission_cmd = 'ssh bioseq@lecs2login "module load python/anaconda_python-3.6.4; python /bioseq/bioSequence_scripts_and_constants/q_submitter.py {} {} -q {} --verbose > {}"'.format(cmds_file, wd, queue_name, log_file)

    # simple command when using shebang header
    submission_cmd = f'ssh bioseq@lecs2login /bioseq/bioSequence_scripts_and_constants/q_submitter.py {cmds_file} {wd} -q bioseq --verbose > {log_file}'

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write(f'\nSSHing and SUBMITting the JOB to the QUEUE:\n{submission_cmd}\n\n')

    subprocess.call(submission_cmd, shell=True)

    if user_email != '':
        with open(os.path.join(wd, 'user_email.txt'), 'w') as f:
            f.write(user_email+'\n')

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write(f'{ctime()}: CGI finished running.\n\n')

except Exception as e:
    msg = 'CGI crashed before the job was submitted :('
    with open(output_path) as f:
        html_content = f.read()
    html_content = html_content.replace('RUNNING', 'FAILED')
    html_content += '<br><br><br><center><h2><font color="red">' + msg + '</font><br><br>Please try to re-run your job or <a href="mailto:bioSequence@tauex.tau.ac.il?subject=ASAP%20Run%20Number%2015249296875723">contact us</a> for further information</h2></center><br><br>\n</body>\n</html>\n'
    with open(output_path, 'w') as f:
        html_content = f.write(html_content)

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        f.write('\n\n' + '$'*60 + '\n\n')
        f.write(ctime() + ': ' + msg + '\n\n')
        f.write(f'{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\n')
        f.write('$'*60)

    #make sure the page will refresh until the end of editing
    with open(output_path) as f:
        html_content = f.read()
    html_content = html_content.replace(CONSTS.RELOAD_TAGS, '')
    with open(output_path, 'w') as f:
        html_content = f.write(html_content)
