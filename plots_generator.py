import matplotlib_venn
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import logging
logger = logging.getLogger('main')

from text_handler import read_table_to_dict

# def plot_nuc_mut_cnt_dict(chain, in_file_path, raw_data_file_suffix):
#     '''plot mutation count frequency'''
#
#     nuc_mut_cnt_dict = read_table_to_dict(in_file_path, key_type=int, value_type=int)
#
#     max_mut_cnt = max(nuc_mut_cnt_dict)  # maximum number of mutations
#     nuc_mut_cnt_values = range(max_mut_cnt + 1)  # list of numbers up until the maximum
#
#     # pad with zeros mutations amounts that are absent
#     for i in range(max_mut_cnt):
#         if i not in nuc_mut_cnt_dict:
#             nuc_mut_cnt_dict[i] = 0
#
#     mut_sum = sum(nuc_mut_cnt_dict.values())
#
#     # create list of frequencies of mutations amount
#     nuc_mut_freq = list()
#     for nuc_mut_cnt in sorted(nuc_mut_cnt_dict):
#         crnt_freq = round((nuc_mut_cnt_dict[nuc_mut_cnt] / mut_sum) * 100, 3)
#         nuc_mut_freq.append(crnt_freq)
#
#     # set plot parameters
#     pos = np.arange(len(nuc_mut_cnt_values))
#     width = 1.0
#
#     ax = plt.axes()
#     ax.set_xticks(pos + (width / 2))
#     ax.set_xticklabels(nuc_mut_cnt_values)
#
#     plt.figure()
#     plt.title(chain + " Frequency of nucleotide mutations amount")
#     plt.xlabel("Number of nucleotide mutations")
#     plt.ylabel("Frequency (%)")
#     plt.bar(pos, nuc_mut_freq, width, color='b')
#     plt.tight_layout()
#
#     fig_name = in_file_path.replace(raw_data_file_suffix, 'png')
#
#     plt.savefig(fig_name)
#     plt.close('all')


def plot_barplot(assignment_file ='/Users/Oren/Dropbox/Projects/wine/outputs/run1/vdj_assignments/IGH_V_counts.txt', raw_data_file_suffix ='txt', key_type = str, value_type=int, as_proportions = True, rotation = None, fontsize = None, ylim=[0, 50]):
    '''plot a simple bar chart of CDR3 lengths and VDJ assignments'''
    d = read_table_to_dict(assignment_file, key_type=key_type, value_type=value_type)

    x_values = sorted(d)
    y_values = [d[x] for x in x_values]

    if as_proportions:
        sum_y_values = sum(y_values)/100 # to get y's in percents
        y_values = [round(y/sum_y_values, 2) for y in y_values]

    logger.debug(x_values)
    logger.debug(y_values)

    # More colors can be found at https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    barplot = sns.barplot(x_values, y_values, color="mediumspringgreen")
    ylim[1] = ylim[1] if ylim[1] > max(y_values) else 1.2*max(y_values)
    barplot.set_ylim(ylim)
    barplot.set_ylabel('Frequency (%)')
    barplot.set_xticklabels(x_values, rotation=rotation, fontsize=fontsize)

    plt.tight_layout()
    save_path = assignment_file.replace(raw_data_file_suffix, 'png')
    plt.savefig(save_path, dpi=500)
    plt.close()


# def plot_intersection_histogram(runs, y_values, out_path):
#     '''plot a simple bar chart of CDR3 lengths and VDJ assignments'''
#
#     logger.debug(runs)
#     logger.debug(y_values)
#     barplot = sns.barplot(np.arange(len(runs)), y_values, color="mediumpurple")
#     barplot.set_xticklabels(runs)
#     barplot.set_title('Percent of intersection between joint and the different runs')
#     barplot.set_ylabel('Frequency (%)')
#     plt.ylim([0, 100])
#
#     plt.savefig(out_path)
#     plt.close()

def plot_venn(out_path, runs_annotations_sets, runs):
    if len(runs) == 2:
        venn = matplotlib_venn.venn2
    elif len(runs) == 3:
        venn = matplotlib_venn.venn3
    else:
        logger.error('Can\'t plot venn diagram for {} sets!! (only for 2 or 3 sets)'.format(len(runs)))
        return
    plt.figure()
    venn(runs_annotations_sets, set_labels=runs)
    plt.title('Venn Diagram')
    plt.savefig(out_path)
    plt.close()


def generate_alignment_report_pie_chart(out_path, isotype_to_precent_dict, minimal_portion=0.05):

    isotypes = sorted(isotype_to_precent_dict, key = isotype_to_precent_dict.get)
    portions = [isotype_to_precent_dict[isotype] for isotype in isotypes]
    others_portion = 0
    others_index = 0
    total = sum(portions)
    for portion in portions:
        if (others_portion + portion)/total > minimal_portion:
            break
        others_portion += portion
        others_index += 1

    isotypes = isotypes[others_index:]
    portions = portions[others_index:]

    if others_portion > 0:
        isotypes.append('others')
        portions.append(others_portion)
    #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice

    fig1, ax1 = plt.subplots()
    ax1.pie(portions, labels=isotypes, autopct='%1.1f%%', startangle=0)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    title = ('run1' if 'run1' in out_path else 'run2') + ' isotype distribution\n'
    plt.title(title)
    plt.savefig(out_path)
    plt.close()    #plt.legend('counts', 'unique')


def plot_correlation(x, y, i, j, out_path):
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y)[0]
    correlation = pearsonr(x, y)[0]
    plt.plot(x, y, 'ro', markersize=2, label='V-gene sequence')
    plt.plot(x, a*x + b, 'b', label='Fitted line (correlation: {})'.format(correlation))
    plt.xlabel('counts in run' + str(i))
    plt.ylabel('counts in run' + str(j))
    plt.title('Pearson correlation of V-gene sequence counts of run{} and run{}'.format(i,j))
    plt.legend()
    plt.savefig(out_path, dpi=500)
    plt.close()    #plt.legend('counts', 'unique')


def generate_clonal_expansion_histogram(cdr3_annotations_path, out_path, cutoff):
    #TODO: Add legend where the higherhisto is named "number of reads" and the lower "clonal expansion"
    #TODO: number of reads
    #TODO: clonal expansion
    #old: cdr3_annotations_path='/Users/Oren/Dropbox/Projects/wine/123_IGH_clones.txt'
    title = 'Top {} common reads distribution VS unique sequences distribution'.format(cutoff)
    cols = np.loadtxt(cdr3_annotations_path, usecols=(1, 2), dtype='int')
    rank = range(1, cutoff+1)
    #for style in plt.style.available:
    style = 'seaborn-paper'
    plt.figure()
    plt.style.use(style)
    plt.figure(figsize=[25,5])
    plt.bar(np.arange(1,len(rank)+1), cols[:cutoff,0], alpha=0.5)
    plt.bar(np.arange(1,len(rank)+1), cols[:cutoff,1], alpha=0.75)
    plt.title(title)
    plt.savefig(out_path)
    plt.close()    #plt.legend('counts', 'unique')


'''
if __name__ == '__main__':
    if len(argv) < 3:
        logger.info('Usage: python ' + argv[0] + ' <report_path> <out_path>; <?title>')
        exit()
    else:
        generate_alignment_report_pie_chart(*argv[1:])
'''
# plot_singles_barplot()
# plot_pairs_barplot()
# plot_triples_barplot()

