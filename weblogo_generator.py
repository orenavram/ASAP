from weblogolib import *
from sys import argv
import logging
logger = logging.getLogger('main')

def generate_weblogo(in_file, out_file, title = "A Logo Title"):
    # fin = open('/Users/Oren/Dropbox/Projects/wine/outputs/clones_analysis/242/IGH/consensus/1_consensus.mafft.fasta')
    # fin = open('/bioseq/data/results/ASAP/15177710318949/outputs/clones_analysis/123/IGH/consensus/1_consensus.mafft.fasta')
    with open(in_file) as f:
        seqs = read_seq_data(f)
    data = LogoData.from_seqs(seqs)
    options = LogoOptions()
    options.title = title
    options.unit_name = 'probability'
    format = LogoFormat(data, options)
    # fout = open('/Users/Oren/Dropbox/Projects/wine/outputs/clones_analysis/242/IGH/consensus/1_consensus.png', 'wb')
    # fout = open('/bioseq/data/results/ASAP/15177710318949/outputs/clones_analysis/123/IGH/consensus/1_consensus.pdf', 'wb')
    with open(out_file, 'wb') as f:
        f.write(pdf_formatter(data, format))

if __name__ == '__main__':
    if len(argv) < 3:
        logger.error('Usage: python ' + argv[0] + ' <in_file> <out_file> <title>')
        exit()
    else:
        generate_weblogo(*argv[1:])


