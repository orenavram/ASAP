#!/shared/python/anaconda3.5/bin/python
import gzip
import os, sys
try:
    import sh
except:
    pass
import cgi, cgitb, traceback
from time import time, ctime
from random import randint
sys.path.append('/bioseq/bioSequence_scripts_and_constants')
sys.path.append('/bioseq/asap')

from email_sender import send_email
from directory_creator import create_dir
import ASAP_CONSTANTS as CONSTS

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
        f.write('''<html><head>{0}
<title>ASAP Job {1}</title>
<link rel="shortcut icon" type="image/x-icon" href="{2}/pics/ASAP_logo.gif" />

<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
<link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">

<link rel="stylesheet" href="{2}/css/general.css">

</head><body>
<nav role="navigation" class="navbar navbar-fixed-top">
    <div class="jumbotron" id="jumbo">
        <div class="container">
            <div class="row" id="title-row">
                <div class="col-md-1">
                </div>
                <div class="col-md-1">
                    <img src="{2}/pics/ASAP_logo.gif" id="antibody_image" class="img-rounded">
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
'''.format(CONSTS.GC.RELOAD_TAGS, run_number, CONSTS.GC.ASAP_URL))


def write_info_paragraph_to_html(output_path):
    with open(output_path, 'a') as f:
        f.write("""<br><div class="container" style="font-size: 20px;" align="justify"> 
<H1 align=center>Job Status - <FONT color='red'>RUNNING</FONT></h1>
<br>ASAP is now processing your request. This page will be automatically updated every {} seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for three months.
</div>
<br>""".format(CONSTS.GC.RELOAD_INTERVAL))


def write_running_parameters_to_html(output_path, job_title, number_of_duplicates, urls_to_reads_files, files_names, MMU, chains, len_threshold, qlty_threshold, number_of_clones_to_analyze, raw_data_suffix, add_mass_spec_seq):

    with open(output_path, 'a') as f:

        #regular params row
        f.write("""<div class="container"><u><h3>Running Parameters:</h3></u><br>""")

        f.write('<div class="row" style="font-size: 20px;">')
        if job_title != '':
            f.write('<div class="col-md-6">')
            f.write('<b>Job title: </b>{}<br><br>'.format(job_title))
            f.write('</div>')
            f.write('</div><div class="row" style="font-size: 20px;">')

        f.write('<div class="col-md-3">')
        f.write('<b>Organism: </b>{}'.format(MMU))
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
            f.write('R1 = <a href=\"{}\" target=_blank>{}</A><br>'.format(urls_to_reads_files[i * 2], files_names[i*2]))
            f.write('R2 = <a href=\"{}\" target=_blank>{}</A><br>'.format(urls_to_reads_files[i * 2 + 1], files_names[i*2+1]))
            f.write('</div>')

        f.write('</div>')

        #Advanced params row
        f.write('<br><u><h3>Advanced Parameters:</h3></u><br>')

        f.write('<div class="row">')
        f.write('<div class="col-md-2">')
        f.write('<b>Min read length: </b>{}<br>'.format(len_threshold))
        f.write('</div>')

        f.write('<div class="col-md-2">')
        f.write('<b>Min read quality: </b>{}<br>'.format(qlty_threshold))
        f.write('</div>')

        f.write('<div class="col-md-2">')
        f.write('<b>Raw data format: </b>{}<br>'.format(raw_data_suffix))
        f.write('</div>')

        f.write('<div class="col-md-2">')
        f.write('<b>Add mass_spec: </b>{}<br>'.format(add_mass_spec_seq))
        f.write('</div>')

        f.write('<div class="col-md-3">')
        f.write('<b>#clones to analyze: </b>{}<br>'.format(number_of_clones_to_analyze))
        f.write('</div>')

        f.write('</div></div>')


def write_pair_file(debug_path, pair, run_content, run_filename, run_dir):
    if run_filename.endswith('gz'):
        #TODO: This option is not working properly. Files are being uploaded but they're non sense :(
        # open_operator = gzip.open
        local_file_name = 'R{}.fastq.gz'.format(pair)
    else:
        local_file_name = 'R{}.fastq'.format(pair)
    open_operator = open
    with open(debug_path, 'a') as f:
        f.write('{} is being handled with {}\n'.format(local_file_name, str(open_operator)))
    local_file_path = os.path.join(run_dir, local_file_name)
    with open_operator(local_file_path, 'wb') as f:
        f.write(run_content)

    #avoid double zipping:
    if local_file_path.endswith('gz'):
        sh.gunzip(local_file_path)

    with open(debug_path, 'a') as f:
        f.write('R{} was handled successfully\n'.format(pair))


def process_uploaded_files(run, debug_path):
    #for debugging
    with open(debug_path, 'a') as f:
        f.write('uploading file of {}\n'.format(run))
    run_R1_filename = form[run + '_R1'].filename
    run_R2_filename = form[run + '_R2'].filename
    # for debugging
    with open(debug_path, 'a') as f:
        f.write('file names are:\n{} (of type {}) and {} (of type {})\n'.format(run_R1_filename, type(form[run + '_R1'].value), run_R2_filename, type(form[run + '_R2'].value),))
    run_R1_content = form[run + '_R1'].value
    run_R2_content = form[run + '_R2'].value
    # for debugging
    with open(debug_path, 'a') as f:
        f.write('{} first 100 chars are: {}\n{} first 100 chars are: {}\n'.format(run_R1_filename, run_R1_content[:100], run_R2_filename, run_R2_content[:100]))
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

    #String (yes/no) that represent whether to add mass_spec_seq to each aa sequence in the fasta output
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

run_number = str(round(time())) + str(randint(1000, 9999))  # adding 4 random figures to prevent users see data that are not their's
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
send_email(smtp_server=CONSTS.GC.SMTP_SERVER, sender='TAU BioSequence <bioSequence@tauex.tau.ac.il>',
           receiver='orenavram@gmail.com', subject='ASAP - A new job has been submitted: {}'.format(run_number), content=os.path.join(CONSTS.GC.ASAP_URL,'results',run_number))

try:
    write_info_paragraph_to_html(output_path)

    with open(cgi_debug_path, 'a') as f:
        # form debugging
        f.write('{}\n{}: A new CGI request has been recieved!\n'.format('#'*50, ctime()))
        f.write('These are the keys that the CGI received:\n{}\n\n'.format('; '.join(sorted(form.keys()))))
        f.write('Form values are:\n')
        for key in sorted(form.keys()):
            if 'run' not in key:
                f.write('{} = {}\n'.format(key, form[key]))
        # for key in sorted(form.keys()):
        #     if 'run' in key:
        #         f.write('100 first characters of {} = '.format(key))
        #         f.write('{}\n'.format(form[key].value[:100]))
        f.write('\n\n')

    #extract form's values:
    user_email = form['email'].value.strip()
    example_page = form['example_page'].value
    chains = ','.join(form.getlist('chains'))
    MMU = form['MMU'].value
    len_threshold = form['len_threshold'].value
    qlty_threshold = form['qlty_threshold'].value
    number_of_clones_to_analyze = form['number_of_clones_to_analyze'].value
    raw_data_suffix = 'xls'
    if form['raw_data_suffix'].value == 'txt':
        raw_data_suffix = 'txt'
    add_mass_spec_seq = 'no'
    #if this option is unchecked, it won't be send in the json (i.e., form['add_mass_spec_seq'].value might not work...)
    if 'add_mass_spec_seq' in form:
        add_mass_spec_seq = 'yes'
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
        f.write('Number of duplicates is {}\n'.format(number_of_duplicates))

    paths_to_reads_files = []
    urls_to_reads_files = []
    for i in range(number_of_duplicates):
        run = 'run' + str(i + 1)
        path_to_reads = os.path.join(wd, 'reads', run)
        url_to_reads = os.path.join(results_url, 'reads', run)
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('{}: Creating path for reads...\n'.format(ctime()))
        create_dir(path_to_reads)
        paths_to_reads_files.append(os.path.join(path_to_reads, 'R1.fastq'))
        paths_to_reads_files.append(os.path.join(path_to_reads, 'R2.fastq'))
        urls_to_reads_files.append(os.path.join(url_to_reads, 'R1.fastq'))
        urls_to_reads_files.append(os.path.join(url_to_reads, 'R2.fastq'))

    if example_page == 'yes':
        files_names = ['rep1_R1.fastq', 'rep1_R2.fastq', 'rep2_R1.fastq', 'rep2_R2.fastq']
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('{}: Copying example files...\n'.format(ctime()))
        # copy example data
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('Fetching: rsync -ravz {} {}\n'.format(CONSTS.EXAMPLE_FILE_RUN1_R1, paths_to_reads_files[0]))
        os.system('rsync -ravz {} {}'.format(CONSTS.EXAMPLE_FILE_RUN1_R1, paths_to_reads_files[0]))
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('{}: ls of {} yields:\n{}\n'.format(ctime(), os.path.join(wd, 'reads', 'run1'), os.listdir(path_to_reads)))
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('Fetching: rsync -ravz {} {}\n'.format(CONSTS.EXAMPLE_FILE_RUN1_R2, paths_to_reads_files[1]))
        os.system('rsync -ravz {} {}'.format(CONSTS.EXAMPLE_FILE_RUN1_R2, paths_to_reads_files[1]))
        os.system('rsync -ravz {} {}'.format(CONSTS.EXAMPLE_FILE_RUN2_R1, paths_to_reads_files[2]))
        os.system('rsync -ravz {} {}'.format(CONSTS.EXAMPLE_FILE_RUN2_R2, paths_to_reads_files[3]))
        with open(cgi_debug_path, 'a') as f: # for cgi debugging
            f.write('{}: Files were copied successfully.\n'.format(ctime()))

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write('{}: ls of {} yields:\n{}\n'.format(ctime(), path_to_reads, os.listdir(path_to_reads)))

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write('{}: write_running_parameters_to_html...\n'.format(ctime()))

    write_running_parameters_to_html(output_path, job_title, number_of_duplicates, urls_to_reads_files, files_names, MMU, chains, len_threshold, qlty_threshold, number_of_clones_to_analyze, raw_data_suffix, add_mass_spec_seq)

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write('{}: Running parameters were written to html successfully.\n'.format(ctime()))

    # This is hidden field that only spammer bots might fill in...
    confirm_email_add = form['confirm_email'].value  # if it is contain a value it is a spammer.

    parameters_file = os.path.join(wd, 'parameters.txt')
    parameters = [wd, CONSTS.GC.MiXCR_dir, number_of_duplicates, chains, len_threshold, qlty_threshold, number_of_clones_to_analyze, MMU, raw_data_suffix, add_mass_spec_seq]
    prepare_parameters_file(parameters_file, parameters)

    cmds_file = os.path.join(wd, 'qsub.cmds')
    write_cmds_file(cmds_file, run_number, parameters_file)

    queue_name = 'bioseq'
    os.system('ssh bioseq@lecs2 python /bioseq/bioSequence_scripts_and_constants/q_submitter.py {} {} {} > {}'.format(cmds_file, wd, queue_name, cmds_file.replace('cmds', 'log')))

    if user_email != '':
        with open(os.path.join(wd, 'user_email.txt'), 'w') as f:
            f.write(user_email)

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        f.write('{}: CGI finished running.\n\n'.format(ctime()))

except Exception as e:
    msg = 'CGI crashed before the job was submitted :('
    with open(output_path) as f:
        html_content = f.read()
    html_content = html_content.replace('RUNNING', 'FAILED')
    html_content += '<br><br><br><center><h2><font color="red">' + msg + '</font><br><br>Please try to re-run your job or <a href="mailto:bioSequence@tauex.tau.ac.il?subject=ASAP%20Run%20Number%2015249296875723">contact us</a> for further information</h2></center><br><br>\n</body>\n</html>\n'

    with open(cgi_debug_path, 'a') as f: # for cgi debugging
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        f.write('\n\n' + '$'*60 + '\n\n')
        f.write(ctime() + ': ' + msg + '\n\n')
        f.write(str(fname) +': ' + str(exc_type) + ', at line: ' + str(exc_tb.tb_lineno) + '\n\n')
        f.write('$'*60)

    #make sure the page will refresh until the end of editing
    with open(output_path) as f:
        html_content = f.read()
    html_content = html_content.replace(CONSTS.GC.RELOAD_TAGS, '')
    with open(output_path, 'w') as f:
        html_content = f.write(html_content)


