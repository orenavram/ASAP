import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
from sys import argv
import logging
logger = logging.getLogger('main')

def generate_reads_histogram(statistics_file_path, out_path, cutoff = 100, title = 'Top 100 common reads distribution VS unique sequences distribution'):
    #statistics_file_path='/Users/Oren/Dropbox/Projects/wine/123_IGH_clones.txt'
    cols = np.loadtxt(statistics_file_path, skiprows=1, usecols=(0, 2, 4), dtype='int')
    rank = cols[:cutoff,0]
    #for style in plt.style.available:
    style = 'seaborn-paper'
    plt.figure()
    plt.style.use(style)
    plt.figure(figsize=[25,5])
    plt.bar(np.arange(1,len(rank)+1), cols[:cutoff,1], alpha=0.5)
    plt.bar(np.arange(1,len(rank)+1), cols[:cutoff,2], alpha=0.75)
    plt.title(title.replace('100', str(cutoff)))
    plt.savefig(out_path)
    plt.close()    #plt.legend('counts', 'unique')

if __name__ == '__main__':
    if len(argv) < 3:
        logger.error('Usage: python ' + argv[0] + ' <statistics_file_path> <out_path>; <?cutoff[100]>')
        exit()
    else:
        generate_reads_histogram(*argv[1:])

