import logging
import os
import sys

logger = logging.getLogger('main')

import ASAP_CONSTANTS as CONSTS


def edit_success_html(gp, html_path, server_main_url, run_number):
    html_text = ''
    if os.path.exists('/bioseq/'): # run on ibis. cgi generated the initial file
        with open(html_path) as f:
            html_text = f.read()
    html_text = html_text.replace('RUNNING', 'FINISHED')
    html_text += '<br><br><br><center><h2>RESULTS:<h2><a href=\'outputs.zip\' target=\'_blank\'><h3><b>Download zipped full results</b></h3></a></center><br><br>\n'
    html_text += '<div><table class="table">'
    html_text += '<thead><tr><th></th>'

    runs = ['run' + str(i + 1) for i in range(gp.number_of_runs)]
    if gp.joint_run_is_needed:
        runs += ['joint']
    for run in runs:
        html_text += '<th align="center">'
        #if os.path.exists(gp.working_dir + '/outputs/' + run):
        html_text += '<h2><b>' + run.replace('run','replicate ').title() + ' results</b></h2>\n'
        html_text += '</th>'

    html_text += '</tr><thead><tbody>'
    for chain in gp.chains:

        html_text += '<tr>'
        html_text += '<td><H3><b>Chain ' + chain + ' results</b></H3></td>\n'
        for run in runs:

            if not os.path.exists(gp.working_dir + '/outputs/' + run):
                logger.info(run + ' does not exist...')
                continue

            html_text += '<td>'
            html_text += '<ul>'
            if chain == 'IGH':
                link = '<a href="outputs/' + run + '/parsed_mixcr_output/alignment_report.png" target="_blank">Isotypes distribution</a>'
                html_text += '<li>' + link + ' ; '
                html_text += '</li>'

            link = '<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.sequence_annotation_file_suffix + '" target="_blank">Sequence annotations</a>'
            html_text += '<li>' + link + ' ;</li>'

            link = 'SHM analysis: <a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutations_file_suffix.replace(
                gp.raw_data_file_suffix,'1.png') + '" target="_blank">nucleotide substitution frequency</a> ; <a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutations_file_suffix.replace(
                gp.raw_data_file_suffix,
                '2.png') + '" target="_blank">Ka Ks analysis</a>'
            raw_link = '(<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutations_file_suffix + '" target="_blank">raw_data</a>)'
            html_text += '<li>' + link + ' ; ' + raw_link + '</li>'

            link = '<a href="outputs/' + run + '/cdr3_analysis/' + chain + '_cdr3_len_counts.png" target="_blank">CDR3 length distribution (AA level)</a>'
            raw_link = '(<a href="outputs/' + run + '/cdr3_analysis/' + chain + '_cdr3_len_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
            html_text += '<li>' + link + ' ; ' + raw_link + '</li>'

            family_distributions = ''
            for group_combination in ['V', 'D', 'J', 'VD', 'VJ', 'DJ', 'VDJ']:
                if chain == 'IGH' or 'D' not in group_combination: # no 'D' fragment in IGK/IGL
                    link = '<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + group_combination + '_counts.png" target="_blank">' + group_combination + ' family subgroup distribution</a>'
                    raw_link = '(<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + group_combination + '_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                    family_distributions += '<li>' + link + ' ; ' + raw_link + '</li>\n'
            html_text += family_distributions

            #if run == 'joint':
            link = '<a href="outputs/'+run+'/cdr3_analysis/' + chain + gp.cdr3_annotation_file_suffix.replace(gp.raw_data_file_suffix, 'png') + '" target="_blank">Clonal expansion graph</a>'
            raw_link = '(<a href="outputs/'+run+'/cdr3_analysis/' + chain + gp.cdr3_annotation_file_suffix + '" target="_blank">raw_data</a>)'
            html_text += '<li>' + link + ' ; ' + raw_link + '</li>\n'

            cdr3_analysis_raw_file = os.path.join(gp.output_path, run, 'cdr3_analysis', chain + gp.top_cdr3_annotation_file_suffix)
            if os.path.exists(cdr3_analysis_raw_file):
                top_cdr3_analysis_html_url = os.path.join(server_main_url, 'results', run_number, 'outputs', run, chain + gp.top_cdr3_annotation_file_suffix).replace(gp.raw_data_file_suffix, 'html')
                link = '<a href="' + top_cdr3_analysis_html_url + '" target="_blank">Top {} clones annotations</a>'.format(gp.top_cdr3_clones_to_further_analyze)
                raw_link = '(<a href="outputs/'+run+'/cdr3_analysis/' + chain + gp.top_cdr3_annotation_file_suffix + '" target="_blank">raw_data</a>)'
                html_text += '<li>' + link + ' ; ' + raw_link + '</li>\n'
                top_cdr3_analysis_html_path = os.path.join(gp.output_path, run, chain + gp.top_cdr3_annotation_file_suffix).replace(gp.raw_data_file_suffix, 'html')
                edit_top_cdr3_analysis_html_page(top_cdr3_analysis_html_path, gp, server_main_url, run_number, chain, run)

            link = '<a href="outputs/' + run + '/' + chain + '_final.fasta" target="_blank">Final fasta</a>'
            html_text += '<li>' + link + ' ;</li>\n\n'

            if run == 'joint':
                link = '<a href="outputs/'+run+'/parsed_mixcr_output/' + chain + '_runs_intersections.png" target="_blank">Runs intersection</a>'
                html_text += '<li>' + link + ' ;</li>\n'

                joint_path = os.path.join(gp.output_path, 'joint')
                for correlation_file in os.listdir(joint_path):
                    if 'correlation' in correlation_file and chain in correlation_file:
                        #correlation_path = os.path.join(joint_path, correlation_file)
                        c_runs = correlation_file.split('_')[:2] #e.g., 'run1_run2_IGH_correlation.png'
                        link = '<a href="outputs/'+run+'/' + correlation_file + '" target="_blank">Correlation of {} and {}</a>'.format(*c_runs)
                        html_text += '<li>' + link + ' ;</li>\n'

            html_text += '</ul>'
            html_text += '</td>'

        html_text += '</tr></tbody>'
    html_text += '</table></div>'
    html_text += '<br><br><br>\n<hr>\n<h4 class=footer><p align=\'center\'>Questions and comments are welcome! Please <span class="admin_link"><a href="mailto:bioSequence@tauex.tau.ac.il?subject=ASAP%20Run%20Number%20' + run_number + '">contact us</a></span></p></h4>\n'
    # html_text += '<div id="bottom_links" align="center"><span class="bottom_link"><a href="' + server_main_url + '" target="_blank">Home</a>&nbsp;|&nbsp<a href="' + server_main_url + 'overview.html" target="_blank">Overview</a>\n</span>\n<br>\n</div>\n</body>\n</html>\n'
    html_text += '<br><br><br>\n</body>\n</html>\n'

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()


def edit_top_cdr3_analysis_html_page(top_cdr3_analysis_html_path, gp, server_main_url, run_number, chain, run):
    with open(top_cdr3_analysis_html_path, 'w') as f:
        annotation_file_path = os.path.join(os.path.split(top_cdr3_analysis_html_path)[0], 'cdr3_analysis', chain + gp.top_cdr3_annotation_file_suffix)
        f.write('''<html>
<head>
    <title>ASAP top clones</title>
    <link rel="shortcut icon" type="image/x-icon" href="{0}/pics/ASAP_logo.gif" />

    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
    <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">
    
    <link rel="stylesheet" href="{0}/css/general.css">

</head>
<body>
    <nav role="navigation" class="navbar navbar-fixed-top">
        <div class="jumbotron" id="jumbo">
            <div class="container">
                <div class="row" id="title-row">
                    <div class="col-md-1">
                    </div>
                    <div class="col-md-1">
                        <img src="{0}/pics/ASAP_logo.gif" id="antibody_image" class="img-rounded">
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
    <div class="container" align="center" style="width: 650px">
        <a href="{1}"><h2>Top {2} clones of {3}, chain {4}</h2></a><br>
        <table class="table">
            <thead>
                <tr>
                    <th>Clone</th>
                    <th>MSA</th>
                    <th>Sequence logo</th>
                    <th>Aligned raw data</th>
                    <th>Unaligned raw data</th>
                </tr>
            </thead>
            <tbody>
        '''.format(CONSTS.ASAP_URL, annotation_file_path, gp.top_cdr3_clones_to_further_analyze, run, chain))
        for i in range(gp.top_cdr3_clones_to_further_analyze):
            wasabi = '<a href="' + server_main_url + '/wasabi/index_general.html?url=' + server_main_url + '/results/' + run_number + '/outputs/'+run+'/cdr3_analysis/'+chain+'_top_' + str(gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(i) + '_msa.aln" target="_blank">+</a> '
            sequence_logo = '<a href="cdr3_analysis/' + chain + '_top_' + str(gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(i) + '_weblogo.pdf" target="_blank">+</a>'
            msa = '<a href="cdr3_analysis/'+chain+'_top_' + str(gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(i) + '_msa.aln" target="_blank">+</a>'
            ms = '<a href="cdr3_analysis/' + chain + '_top_' + str(gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(i) + '_ms.fasta" target="_blank">+</a>'
            f.write('''<tr>
    <td align="center">{}</td>
    <td align="center">{}</td>
    <td align="center">{}</td>
    <td align="center">{}</td>
    <td align="center">{}</td>
</tr>'''.format(i+1, wasabi, sequence_logo, msa, ms))
        f.write('</tbody></table></body></html>')
        f.flush()


def edit_failure_html(html_path, msg):
    html_text = ''
    if os.path.exists(CONSTS.SERVERS_RESULTS_DIR): # run on ibis. cgi generated the initial file
        with open(html_path) as f:
            html_text = f.read()
    html_text = html_text.replace('RUNNING', 'FAILED')
    html_text +='<br><br><br>'
    html_text +='<center><h2>'
    html_text +='<font color="red">{}</font><br><br>'.format(msg)
    html_text +='Please try to re-run your job or <a href="mailto:bioSequence@tauex.tau.ac.il?subject=ASAP%20Run%20Number%2015249296875723">contact us</a> for further information'
    html_text +='</h2></center><br><br>'
    html_text +='\n</body>\n</html>\n'
    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()