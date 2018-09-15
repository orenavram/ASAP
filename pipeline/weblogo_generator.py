from weblogolib import *
from sys import argv
import logging
logger = logging.getLogger('main')

def generate_weblogo(in_file, out_file, title = "A Logo Title"):
    logger.info(f'Generating logo for {out_file}')
    with open(in_file) as f:
        seqs = read_seq_data(f)
    data = LogoData.from_seqs(seqs)
    options = LogoOptions()
    options.title = title
    options.unit_name = 'probability'
    format = LogoFormat(data, options)
    with open(out_file, 'wb') as f:
        f.write(pdf_formatter(data, format))

if __name__ == '__main__':
    if len(argv) < 3:
        logger.error('Usage: python ' + argv[0] + ' <in_file> <out_file> <title>')
        exit()
    else:
        generate_weblogo(*argv[1:])


