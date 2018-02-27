import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import logging
logger = logging.getLogger('main')

from text_handler import read_table_to_dict

def plot_nuc_mut_cnt_dict(chain, in_file_path, raw_data_file_suffix):
    '''plot mutation count frequency'''

    nuc_mut_cnt_dict = read_table_to_dict(in_file_path, key_type=int, value_type=int)

    max_mut_cnt = max(nuc_mut_cnt_dict)  # maximum number of mutations
    nuc_mut_cnt_values = range(max_mut_cnt + 1)  # list of numbers up until the maximum

    # pad with zeros mutations amounts that are absent
    for i in range(max_mut_cnt):
        if i not in nuc_mut_cnt_dict:
            nuc_mut_cnt_dict[i] = 0

    mut_sum = sum(nuc_mut_cnt_dict.values())

    # create list of frequencies of mutations amount
    nuc_mut_freq = list()
    for nuc_mut_cnt in sorted(nuc_mut_cnt_dict):
        crnt_freq = round((nuc_mut_cnt_dict[nuc_mut_cnt] / mut_sum) * 100, 3)
        nuc_mut_freq.append(crnt_freq)

    # set plot parameters
    pos = np.arange(len(nuc_mut_cnt_values))
    width = 1.0

    ax = plt.axes()
    ax.set_xticks(pos + (width / 2))
    ax.set_xticklabels(nuc_mut_cnt_values)

    plt.figure()
    plt.title(chain + " Frequency of nucleotide mutations amount")
    plt.xlabel("Number of nucleotide mutations")
    plt.ylabel("Frequency (%)")
    plt.bar(pos, nuc_mut_freq, width, color='b')
    plt.tight_layout()

    fig_name = in_file_path.replace(raw_data_file_suffix, 'png')

    plt.savefig(fig_name)
    plt.close('all')


def plot_chart(assignment_file, raw_data_file_suffix, key_type = str, as_proportions = True):
    '''plot a simple bar chart of CDR3 lengths and VDJ assignments'''
    # initialize chart variables

    d = read_table_to_dict(assignment_file, key_type=key_type, value_type=int)

    x_values = sorted(d)
    y_values = [d[x] for x in x_values]

    if as_proportions:
        sum_y_values = sum(y_values)/100 # to get y's in percents
        y_values = [round(y/sum_y_values, 2) for y in y_values]


    # create chart
    fig, ax = plt.subplots(figsize=(len(x_values),len(x_values)/1.5))

    bars_loc = np.arange(len(x_values))
    bars = ax.bar(bars_loc, y_values, color='blue')
    ax.set_xticks(bars_loc)  # align xticks under bars

    # set chart parameters
    # ax.set_title(sample_num + ' ' + chain + ' Frequency of ' + x_label)
    ax.set_ylabel('Frequency')
    ax.set_xticklabels(x_values, rotation=70)
    ax.set_ylim([0, 100])

    # label each bar in chart with its value
    for bar, label in zip(bars, y_values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 1, label,
                ha='center', va='bottom', rotation=60)

    save_path = assignment_file.replace(raw_data_file_suffix, 'png')

    fig.savefig(save_path)

    plt.close('all')
    plt.ioff()
    #4/0


def plot_intersection_histogram(runs, y_values, out_path):
    '''plot a simple bar chart of CDR3 lengths and VDJ assignments'''

    logger.info(runs)
    logger.info(y_values)
    plt.title('Percent of intersection between joint and the different runs')
    plt.ylabel('Frequency (%)')
    plt.bar(range(len(runs)), y_values, color='b', tick_label=runs, align='center')
    plt.ylim([0, 100])

    plt.savefig(out_path)
    plt.close()


def generate_alignment_report_pie_chart(out_path, isotype_to_precent_dict, minimal_portion=0.05):

    '''
        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        isotypes = []
        portions = []
        portion_to_isotype = {}
        #subtypes = {} # isotype -> list of subtypes
        with open(in_path) as f:
            for line in f:
                match = search('% of ([a-zA-Z]+) isotypes? = (\d+\.\d+)', line)
                if not match:
                    continue
                isotype = match.group(1)
                #subtype = match.group(2)
                #subtypes[isotype] = subtypes.get(isotype, []) + [subtype]
                portion = float(match.group(2))
                isotypes.append(isotype)
                portions.append(portion)
                #portion_to_isotype[portion] = isotype
    '''

    #remove zero counts

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

    if others_portion:
        isotypes.append('others')
        portions.append(others_portion)
    #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice

    fig1, ax1 = plt.subplots()
    ax1.pie(portions, labels=isotypes, autopct='%1.1f%%', startangle=0)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    title = ('run1' if 'run1' in out_path else 'run2') + ' isotype pie chart\n'
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