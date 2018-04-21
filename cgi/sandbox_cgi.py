#!/shared/python/anaconda3.5/bin/python

#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!
#THIS IS SANDBOX!!!!!!

import os, sys
import cgi, cgitb

sys.path.append('/bioseq/bioSequence_scripts_and_constants')
sys.path.append('/bioseq/asap')

from directory_creator import create_dir
import ASAP_CONSTANTS as CONSTS


def hello_world(output_path='', run_number='NO_RUN_NUMBER'):
    with open(output_path, 'w') as f:
        f.write("""<html><body><h1 align=center><FONT color='red'>""" + run_number + """ was POSTed successfully!!!</FONT></h1>""")

def write_running_parameters_to_html(output_path, form):
    with open(output_path, 'a') as f:

        f.write("""<font face=Verdana><u><h4>Running Parameters:</h4></u></font>""")
        for key in form:
            if 'run' not in key:
                f.write('<ul>{} = {}<br></ul>'.format(key, form[key]))
        f.write('<br><br>')
        f.write('</body></html>')

# def print_running_parameters_to_html(form):
#     print("""<font face=Verdana><u><h4>Running Parameters:</h4></u></font>""")
#     for key in form:
#         if 'run' not in key:
#             print('<ul>{} = {}<br></ul>'.format(key, form[key]))
#     print('<br><br>')
#     print('</body></html>')
#
# def print_hello_world(run_number='NO_RUN_NUMBER'):
#     print("""<html><body><h1 align=center><FONT color='red'>""" + run_number + """ was POSTed successfully!!!</FONT></h1>""")


# prints detailed error report on BROWSER when cgi crashes
# This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
cgitb.enable()

# print_hello_world() # for debugging
form = cgi.FieldStorage()  # extract POSTed object

run_number = 'sandbox'

results_url = os.path.join(CONSTS.ASAP_RESULTS_URL, run_number)
output_url = os.path.join(results_url, 'output.php')

wd = os.path.join(CONSTS.ASAP_RESULTS_DIR, run_number)
create_dir(wd)
output_path = os.path.join(wd, 'output.php')
cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')

hello_world(output_path, run_number)
write_running_parameters_to_html(output_path, form) # html's prefix must be written BEFORE redirecting...

print('Location: ' + output_url)  # Redirects to the results url. MUST appear before any other print.
print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
sys.stdout.flush()  # must be flushed immediately!!!

# print_hello_world(output_path, run_number)
# print_running_parameters_to_html(form) # html's prefix must be written BEFORE redirecting...





