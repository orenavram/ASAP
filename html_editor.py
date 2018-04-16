import os
import logging
logger = logging.getLogger('main')

def edit_html(gp, html_path, html_mode, server_main_url, run_number):
    with open(html_path, html_mode) as f:
        html_finish_header = '<br><br><br><center><H1><a name=finish>ASAP calculation is finished:</a></H1>\n'
        html_finish_header += '<a href=\'outputs.zip\' target=\'_blank\'><h3><b>Zipped full results</b></h3></a></center><br><br>\n'
        f.write(html_finish_header)
        f.write('<div><table align="center">')
        f.write('<tr><td></td>')

        runs = ['run' + str(i + 1) for i in range(gp.number_of_runs)] + ['joint']
        for run in runs:
            f.write('<td align="center">')
            if os.path.exists(gp.working_dir + '/outputs/' + run):
                f.write('<H2><b>' + run.title() + ' results</b></H2>\n')
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
                if run != 'joint':
                    link = '<a href="outputs/' + run + '/parsed_mixcr_output/alignment_report.png" target="_blank">Alignment report pie chart</a>'
                    raw_link = '(<a href="outputs/' + run + '/parsed_mixcr_output/alignment_report.log" target="_blank">raw_data</a>)'
                    f.write('<li>' + link + ' ; ' + raw_link + '</li>')

                link = '<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.sequence_annotation_file_suffix + '" target="_blank">Sequence annotations</a>'
                f.write('<li>' + link + ' ;</li>')

                link = '<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutation_count_file_suffix.replace(
                    gp.raw_data_file_suffix,
                    'png') + '" target="_blank">Somatic hypermutations - nucleotide substitution frequency</a>'
                raw_link = '(<a href="outputs/' + run + '/parsed_mixcr_output/' + chain + gp.mutation_count_file_suffix + '" target="_blank">raw_data</a>)'
                f.write('<li>' + link + ' ; ' + raw_link + '</li>')

                link = '<a href="outputs/' + run + '/vdj_assignments/' + chain + '_cdr3_len_counts.png" target="_blank">CDR3 length distribution (AA level)</a>'
                raw_link = '(<a href="outputs/' + run + '/vdj_assignments/' + chain + '_cdr3_len_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                f.write('<li>' + link + ' ; ' + raw_link + '</li>')

                singles = ''
                for char in 'VDJ':
                    link = '<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + char + '_counts.png" target="_blank">' + char + ' families distribution</a>'
                    raw_link = '(<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + char + '_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                    singles += '<li>' + link + ' ; ' + raw_link + '</li>\n'
                f.write(singles)

                pairs = ''
                for chars_pair in ['VD', 'VJ', 'DJ']:
                    link = '<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + chars_pair + '_counts.png" target="_blank">' + chars_pair + ' families distribution</a>'
                    raw_link = '(<a href="outputs/' + run + '/vdj_assignments/' + chain + '_' + chars_pair + '_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                    pairs += '<li>' + link + ' ; ' + raw_link + '</li>\n'
                f.write(pairs)

                link = '<a href="outputs/' + run + '/vdj_assignments/' + chain + '_VDJ_counts.png" target="_blank">VDJ families distribution</a>'
                raw_link = '(<a href="outputs/' + run + '/vdj_assignments/' + chain + '_VDJ_counts.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                f.write('<li>' + link + ' ; ' + raw_link + '</li>\n')

                if run == 'joint':
                    link = '<a href="outputs/joint/parsed_mixcr_output/' + chain + '_runs_intersections.png" target="_blank">Runs intersection</a>'
                    f.write('<li>' + link + ' ;</li>\n')

                    joint_path = os.path.join(gp.output_path, 'joint')
                    for correlation_file in os.listdir(joint_path):
                        if 'correlation' in correlation_file:
                            #correlation_path = os.path.join(joint_path, correlation_file)
                            runs = correlation_file.split('_')[:2] #e.g., 'run1_run2_IGH_correlation.png'
                            link = '<a href="outputs/joint/' + correlation_file + '" target="_blank">Correlation of {} and {}</a>'.format(*runs)
                            f.write('<li>' + link + ' ;</li>\n')

                    link = '<a href="outputs/joint/parsed_mixcr_output/' + chain + '_runs_intersections.png" target="_blank">Runs intersection</a>'
                    f.write('<li>' + link + ' ;</li>\n')

                    link = '<a href="outputs/joint/' + chain + '_final.fasta" target="_blank">Final joint fasta</a>'
                    f.write('<li>' + link + ' ;</li>\n\n')

                    link = '<a href="outputs/joint/cdr3_analysis/' + chain + gp.cdr3_annotation_file_suffix.replace(gp.raw_data_file_suffix, 'png') + '" target="_blank">Clonal expansion graph</a>'
                    raw_link = '(<a href="outputs/joint/cdr3_analysis/' + chain + gp.cdr3_annotation_file_suffix + '" target="_blank">raw_data</a>)'
                    f.write('<li>' + link + ' ; ' + raw_link + '</li>\n')

                    link = '<a href="outputs/joint/cdr3_analysis/' + chain + gp.top_cdr3_annotation_file_suffix + '" target="_blank">Top {} clones annotations</a>'.format(gp.top_cdr3_clones_to_further_analyze)
                    f.write('<li>' + link + ' ;</li>\n')

                    #todo: move this to a blank window
                    for i in range(gp.top_cdr3_clones_to_further_analyze):
                        msa_link = '<a href="' + server_main_url + 'wasabi/index_general.html?url=' + server_main_url + 'results/' + run_number + '/outputs/joint/cdr3_analysis/cluster_' + str(
                            i) + '_msa.aln" target="_blank">MSA</a> '
                        msa_link += '(raw_data: <a href="outputs/joint/cdr3_analysis/cluster_' + str(
                            i) + '_ms.fasta" target="_blank">unaligned</a>, '
                        msa_link += '<a href="outputs/joint/cdr3_analysis/cluster_' + str(
                            i) + '_msa.aln" target="_blank">aligned</a>)'
                        # consensus_link = '<a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/'+str(i)+'_consensus.fasta" target="_blank">Consensus</a>'
                        weblogo_link = '<a href="outputs/joint/cdr3_analysis/cluster_' + str(
                            i) + '_weblogo.pdf" target="_blank">Logo</a>'
                        f.write(
                            '<li>Clone #' + str(i) + ' cluster: ' + ' ; '.join([msa_link, weblogo_link]) + '</li>\n')

                f.write('</ul>')
                '''
                #TODO: generate these files for each run....
                if run == 'final_analysis':
                    f.write('</td></tr><tr><td></td><td></td><td></td><td>')
                    f.write('<center><H3><b>Clones analysis - '+chain+'</b></H3></center>\n')
                    f.write('<ul>\n')
                    reads_frequency_distribution = '<a href="outputs/clones_analysis/' + sample + '/' + chain + '/' + chain + '_clones.png" target="_blank">Clones statistics distribution</a>'
                    reads_frequency_excel = '(<a href="outputs/clones_analysis/' + sample + '/' + chain + '/' + chain + '_clones.' + gp.raw_data_file_suffix + '" target="_blank">raw_data</a>)'
                    f.write('<li>' + reads_frequency_distribution + ' ' + reads_frequency_excel + '</li>')
                    for i in range(1, gp.top_cdr_cnt+1):
                        msa_link = '<a href="'+server_main_url+'/wasabi/index_general.html?url='+server_main_url+'/results/'+run_number+'/outputs/clones_analysis/'+sample+'/'+chain+'/consensus/'+str(i)+'_cluster_msa.aln" target="_blank">MSA</a> '
                        msa_link += '(raw_data: <a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/'+str(i)+'_cluster.fasta" target="_blank">unaligned</a>, '
                        msa_link += '<a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/'+str(i)+'_cluster_msa.aln" target="_blank">aligned</a>)'
                        consensus_link ='<a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/'+str(i)+'_consensus.fasta" target="_blank">Consensus</a>'
                        weblogo_link = '<a href="outputs/clones_analysis/' + sample + '/'+chain+'/consensus/'+str(i)+'_cluster_msa.pdf" target="_blank">Logo</a>'
                        f.write('<li>Clone #' + str(i) + ' cluster:' + ' ' + ' ; '.join([msa_link, consensus_link, weblogo_link])+'</li>\n')
                    f.write('</ul>\n')
                '''
                f.write('</td>')

            f.write('</tr>')
        f.write('</table></div>')
        html_footer = '<br><br><br>\n<hr>\n<h4 class=footer><p align=\'center\'>Questions and comments are welcome! Please <span class="admin_link"><a href="mailto:bioSequence@tauex.tau.ac.il?subject=ASAP%20Run%20Number%20' + run_number + '">contact us</a></span></p></h4>\n'
        html_footer += '<div id="bottom_links" align="center"><span class="bottom_link"><a href="' + server_main_url + '" target="_blank">Home</a>&nbsp;|&nbsp<a href="' + server_main_url + 'overview.php" target="_blank">Overview</a>\n</span>\n<br>\n</div>\n</body>\n</html>\n'
        f.write(html_footer)
