import matplotlib_venn
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import logging
logger = logging.getLogger('main')

from text_handler import read_table_to_dict


def plot_barplot(d, out_path, title='', x_label='', y_label='Frequency (%)\n', as_proportions = True, rotation = 90, ylim=[0, 50]):
    #key_type = str, value_type=int):
    '''plot a simple bar chart of CDR3 lengths and VDJ assignments'''
    #d = read_table_to_dict(assignment_file, key_type=key_type, value_type=value_type)

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
    #print(x_label, len(x_values))
    if len(x_values) < 20:
        fontsize = 10
    elif len(x_values) < 60:
        fontsize = 5
    else:
        fontsize = 2
    barplot.set_ylim(ylim)
    barplot.set_xticklabels(x_values, rotation=rotation, fontsize=fontsize)
    barplot.set_xlabel(x_label)
    barplot.set_ylabel(y_label)
    barplot.set_title(title)

    plt.tight_layout()
    plt.savefig(out_path, dpi=500, bbox_inches='tight')
    plt.close()


def generate_clonal_expansion_histogram(cdr3_annotations_path, out_path, cutoff, fontsize=15):
    title = f'{cutoff} top ranking clones prevalence\n'
    cols = np.loadtxt(cdr3_annotations_path, usecols=(1, 2), dtype='int')
    rank = range(1, cutoff+1)
    for style in plt.style.available:
        plt.figure(figsize=[25,5])
        plt.style.use(style)
        plt.bar(np.arange(1,len(rank)+1), cols[:cutoff,0], label='Reads counts')
        plt.bar(np.arange(1,len(rank)+1), cols[:cutoff,1], label='Clonal expansion')
        plt.title(title, fontsize=fontsize)
        plt.xlabel(f'\nClone number', fontsize=fontsize)
        plt.ylabel(f'Counts\n', fontsize=fontsize)
        plt.xticks(range(0,cutoff+1,10), fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(loc='best', fontsize=fontsize)
        if style=='seaborn-paper':
            plt.savefig(out_path, dpi=500, bbox_inches='tight')
        plt.close()
        if style=='seaborn-paper':
            break

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
        logger.error(f'Can\'t plot venn diagram for {len(runs)} sets!! (only for 2 or 3 sets)')
        return
    plt.figure()
    venn(runs_annotations_sets, set_labels=runs)
    plt.title('Venn Diagram')
    plt.savefig(out_path, dpi=500, bbox_inches='tight')
    plt.close()


def generate_alignment_report_pie_chart(out_path, isotype_to_precent_dict, run=None):

    isotypes = sorted(isotype_to_precent_dict, key=isotype_to_precent_dict.get, reverse=True)
    portions = [isotype_to_precent_dict[isotype] for isotype in isotypes]

    portions_percents = [100*portions[i]/sum(portions) for i in range(len(portions)) if portions[i]!=0]
    isotypes = [isotypes[i] for i in range(len(isotypes)) if portions[i]!=0]
    labels = [f'{isotypes[i]} ({portions_percents[i]:.3f}%)' for i in range(len(portions_percents))]

    patches, texts = plt.pie(portions_percents, counterclock=False)
    plt.legend(patches, labels, loc="best")
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.tight_layout()
    if not run:
        run = 'Replicate '+ out_path[out_path.find('run')+len('run')]
    title = f'{run} isotype distribution'
    plt.title(title)
    plt.savefig(out_path, dpi=500, bbox_inches='tight')
    plt.close()


def plot_correlation(x, y, i, j, out_path):
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y)[0]
    correlation = pearsonr(x, y)[0]
    plt.plot(x, y, 'ro', markersize=2, label='V-gene sequence')
    plt.plot(x, a*x + b, 'b', label=f'Fitted line (r={correlation:.3f})')
    plt.xlabel(f'\nCounts in replicate {i}')
    plt.ylabel(f'Counts in replicate {j}\n')
    plt.title(f'Pearson correlation of V-gene sequence counts of replicate {i} and replicate {j}')
    plt.legend()
    plt.savefig(out_path, dpi=500, bbox_inches='tight')
    plt.close()


def generate_mutations_boxplots(core_dna_to_num_of_mutations_dict, out_path):
    fig = plt.figure(figsize=(10,10))

    # Create an axes instance
    ax = fig.add_subplot(221)
    data = np.array(list(core_dna_to_num_of_mutations_dict.values()))
    plot_sub_figure(ax, data, ['Ka', 'Ks'])

    # Create an axes instance
    ax = fig.add_subplot(222)
    Ka_Ks = data[:,0]/data[:,1]
    plot_sub_figure(ax, Ka_Ks, ['Ka/Ks'])

    fig.savefig(out_path, dpi=500, bbox_inches='tight')
    plt.close()


def plot_sub_figure(ax, data, xticklabels):
    ## to get fill color
    bp = ax.boxplot(data, patch_artist=True, showfliers=False)
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set(color='#7570b3', linewidth=2)
        # change fill color
        box.set(facecolor='#1b9e77')

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    ax.set_xticklabels(xticklabels)

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

