from sys import argv
import logging
logger = logging.getLogger('main')

def parse_msa(in_path):
    headers = []
    sequences = []
    with open(in_path) as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:].rstrip()) # add new header
                sequences.append('') # initialize new sequence
            else:
                sequences[-1] += line.rstrip() #update current sequence
    return headers, sequences


def remove_sparse_columns(in_path, out_path, maximal_gap_frequency_allowed_per_column=0.2):
    maximal_gap_frequency_allowed_per_column = float(maximal_gap_frequency_allowed_per_column)
    logger.info('Removing sparse column (gaps frequency > {}) from {}'.format(maximal_gap_frequency_allowed_per_column, in_path))
    headers, sequences = parse_msa(in_path)
    logger.debug('SEQUENCES: '+str(sequences))
    cleaned_sequences = [''] * len(sequences)
    for j in range(len(sequences[0])):
        column_j = [sequences[i][j] for i in range(len(sequences))]
        gap_frequency = column_j.count('-')/len(sequences)
        if gap_frequency < maximal_gap_frequency_allowed_per_column: # not a sparse column
            for i in range(len(cleaned_sequences)): # add j'th column
                cleaned_sequences[i] += sequences[i][j]
        else:
            logger.debug('Removing column #{}: {}'.format(j, column_j))

    result = ''
    for i in range(len(cleaned_sequences)):
        result += '>' + headers[i] + '\n' + cleaned_sequences[i] + '\n'

    with open(out_path, 'w') as f:
        f.write(result)

    logger.info('Shortened from %d to %d columns' % (len(sequences[0]), len(cleaned_sequences[0])))


if __name__ == '__main__':
    if len(argv) < 3:
        logger.error('Usage: python ' + argv[0] + ' <alignment_path> <out_path>; <?maximal_gap_frequency_allowed_per_column[0.2]>')
        exit()
    else:
        remove_sparse_columns(*argv[1:])