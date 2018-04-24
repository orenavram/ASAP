import os
import logging
logger = logging.getLogger('main')

def edit_success_html(gp, html_path, html_mode, server_main_url, run_number):
    with open(html_path, html_mode) as f:
        html_finish_header = '<br><br><br><center><H1><a name=finish>ASAP calculation is finished:</a></H1>\n'
        html_finish_header += '<a href=\'outputs.zip\' target=\'_blank\'><h3><b>Zipped full results</b></h3></a></center><br><br>\n'
        f.write(html_finish_header)
        f.write('<div><table align="center">')
        f.write('<tr><td></td>')

        runs = ['run' + str(i + 1) for i in range(gp.number_of_runs)] + ['joint']
        for run in runs:
            f.write('<td align="center">')
            #if os.path.exists(gp.working_dir + '/outputs/' + run):
            f.write('<H2><b>' + run.replace('run','replicate ').title() + ' results</b></H2>\n')
            f.write('</td>')

        f.write('</tr>')
        for chain in gp.chains:

            f.write('<tr>')
            f.write('<td><H3><b>Chain ' + chain + ' results</b></H3></td>\n')
            for run in runs:

                if not os.path.exists(gp.working_dir + '/outputs/' + run):
                    logger.info(run + ' does not exist...')
                    continue

                f.write('<td>')
                f.write('<ul>')
                link = '<a href="outputs/' + run + '/parsed_mixcr_output/alignment_report.png" target="_blank">Alignment report pie chart</a>'
                f.write('<li>' + link + ' ; ')
                # if run != 'joint': #TODO: maybe this block should be removed?
                #     raw_link = '(<a href="outputs/' + run + '/parsed_mixcr_output/alignment_report.log" target="_blank">raw_data</a>)'
                #     f.write(raw_link)
                f.write('</li>')

                link = '<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.sequence_annotation_file_suffix + '" target="_blank">Sequence annotations</a>'
                f.write('<li>' + link + ' ;</li>')

                link = '<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutations_file_suffix.replace(
                    gp.raw_data_file_suffix,
                    '1.png') + '" target="_blank">SHM analysis: nucleotide substitution frequency</a> ; <a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutations_file_suffix.replace(
                    gp.raw_data_file_suffix,
                    '2.png') + '" target="_blank">Ka Ks analysis</a>'
                raw_link = '(<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutations_file_suffix + '" target="_blank">raw_data</a>)'
                f.write('<li>' + link + ' ; ' + raw_link + '</li>')

                link = '<a href="outputs/' + run + '/cdr3_analysis/' + chain + '_cdr3_len_counts.png" target="_blank">CDR3 length distribution (AA level)</a>'
                raw_link = '(<a href="outputs/' + run + '/cdr3_analysis/' + chain + '_cdr3_len_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                f.write('<li>' + link + ' ; ' + raw_link + '</li>')

                family_distributions = ''
                for group_combination in ['V', 'D', 'J', 'VD', 'VJ', 'DJ', 'VDJ']:
                    link = '<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + group_combination + '_counts.png" target="_blank">' + group_combination + ' family subgroup distribution</a>'
                    raw_link = '(<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + group_combination + '_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                    family_distributions += '<li>' + link + ' ; ' + raw_link + '</li>\n'
                f.write(family_distributions)

                #if run == 'joint':
                link = '<a href="outputs/'+run+'/cdr3_analysis/' + chain + gp.cdr3_annotation_file_suffix.replace(gp.raw_data_file_suffix, 'png') + '" target="_blank">Clonal expansion graph</a>'
                raw_link = '(<a href="outputs/'+run+'/cdr3_analysis/' + chain + gp.cdr3_annotation_file_suffix + '" target="_blank">raw_data</a>)'
                f.write('<li>' + link + ' ; ' + raw_link + '</li>\n')

                cdr3_analysis_raw_file = os.path.join(gp.output_path, run, 'cdr3_analysis', chain + gp.top_cdr3_annotation_file_suffix)
                if os.path.exists(cdr3_analysis_raw_file):
                    top_cdr3_analysis_html_path = os.path.join(gp.output_path, run, chain + gp.top_cdr3_annotation_file_suffix).replace(gp.raw_data_file_suffix, 'html')
                    link = '<a href="' + top_cdr3_analysis_html_path + '" target="_blank">Top {} clones annotations</a>'.format(gp.top_cdr3_clones_to_further_analyze)
                    raw_link = '(<a href="outputs/'+run+'/cdr3_analysis/' + chain + gp.top_cdr3_annotation_file_suffix + '" target="_blank">raw_data</a>)'
                    f.write('<li>' + link + ' ; ' + raw_link + '</li>\n')
                    edit_top_cdr3_analysis_html_page(top_cdr3_analysis_html_path, gp, server_main_url, run_number, chain, run)

                if run == 'joint':
                    link = '<a href="outputs/'+run+'/parsed_mixcr_output/' + chain + '_runs_intersections.png" target="_blank">Runs intersection</a>'
                    f.write('<li>' + link + ' ;</li>\n')

                    joint_path = os.path.join(gp.output_path, 'joint')
                    for correlation_file in os.listdir(joint_path):
                        if 'correlation' in correlation_file:
                            #correlation_path = os.path.join(joint_path, correlation_file)
                            c_runs = correlation_file.split('_')[:2] #e.g., 'run1_run2_IGH_correlation.png'
                            link = '<a href="outputs/'+run+'/' + correlation_file + '" target="_blank">Correlation of {} and {}</a>'.format(*c_runs)
                            f.write('<li>' + link + ' ;</li>\n')

                    link = '<a href="outputs/'+run+'/' + chain + '_final.fasta" target="_blank">Final '+run+' fasta</a>'
                    f.write('<li>' + link + ' ;</li>\n\n')

                f.write('</ul>')
                f.write('</td>')

            f.write('</tr>')
        f.write('</table></div>')
        html_footer = '<br><br><br>\n<hr>\n<h4 class=footer><p align=\'center\'>Questions and comments are welcome! Please <span class="admin_link"><a href="mailto:bioSequence@tauex.tau.ac.il?subject=ASAP%20Run%20Number%20' + run_number + '">contact us</a></span></p></h4>\n'
        html_footer += '<div id="bottom_links" align="center"><span class="bottom_link"><a href="' + server_main_url + '" target="_blank">Home</a>&nbsp;|&nbsp<a href="' + server_main_url + 'overview.php" target="_blank">Overview</a>\n</span>\n<br>\n</div>\n</body>\n</html>\n'
        f.write(html_footer)


def edit_top_cdr3_analysis_html_page(top_cdr3_analysis_html_path, gp, server_main_url, run_number, chain, run):
    with open(top_cdr3_analysis_html_path, 'w') as f:
        for i in range(gp.top_cdr3_clones_to_further_analyze):
            msa_link = '<a href="' + server_main_url + 'wasabi/index_general.html?url=' + server_main_url + 'results/' + run_number + '/cdr3_analysis/'+chain+'_top_' + str(
                gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(
                i) + '_msa.aln" target="_blank">MSA</a> '
            msa_link += '(raw_data: <a href="cdr3_analysis/'+chain+'_top_' + str(
                gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(
                i) + '_ms.fasta" target="_blank">unaligned</a>, '
            msa_link += '<a href="cdr3_analysis/'+chain+'_top_' + str(
                gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(
                i) + '_msa.aln" target="_blank">aligned</a>)'
            # consensus_link = '<a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/'+str(i)+'_consensus.fasta" target="_blank">Consensus</a>'
            weblogo_link = '<a href="cdr3_analysis/'+chain+'_top_' + str(
                gp.top_cdr3_clones_to_further_analyze) + '_clones/cluster_' + str(
                i) + '_weblogo.pdf" target="_blank">Logo</a>'
            f.write(
                '<li>Clone #' + str(i+1) + ' cluster: ' + ' ; '.join([msa_link, weblogo_link]) + '</li>\n')


def edit_failure_html(html_path, html_mode, msg):
    with open(html_path, html_mode) as f:
        f.write('<br><br><br>')
        f.write('<center><H1><a name=finish>')
        f.write(msg)
        f.write('Please try to re-run your job or contact us for further information')
        f.write('</a></H1></center>')