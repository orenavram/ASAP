import logging
from difflib import SequenceMatcher

logger = logging.getLogger('main')

def read_table_to_dict(file_to_read, key_type = str, value_type = str, d = None, delimiter ='\t', reverse = False, open_operator=open):
    '''
    parse a delimited file with two columns to a dictionary
    if reverse == 0 keys are taken from the first column o.w., from the second
    open_operator can also be set to gzip.open if the file is zipped
    '''
    logger.debug('Reading dict from: ' + file_to_read)
    if d == None:
        d = {}
    with open_operator(file_to_read) as f:
        for line in f:
            if line != '' and not line.isspace():
                item1, item2 = line.rstrip().split(delimiter)
                d[key_type(item1)] = value_type(item2)
    if reverse:
        d = dict((d[key],key) for key in d)
    return d #, dict([(barcode_to_name[barcode], barcode) for barcode in barcode_to_name])


def write_dict_to_file(out_path, d, delimiter='\t', sort_by=None, reverse=False):
    logger.info('Writing dict to: ' + out_path)
    if d:
        with open(out_path, 'w') as f:
            for key in sorted(d, key=sort_by, reverse=reverse):
                f.write(delimiter.join([str(key), str(d[key])]) + '\n')
    else:
        logger.info('Skipping. Dictionary is empty...')


def get_next_relevant_rstriped_line(file_handler):
    while True:
        line = file_handler.readline().strip()
        if line!='' and line[0]!='#':
            logger.info('Parsing next relevant line:' + line)
            return line
        else:
            logger.info('Skipping irrelevant (i.e., comment/white spaces) line:' + line)

def string_similarity(a, b):
    # measure similarity rate of two strings
    return SequenceMatcher(None, a, b).ratio()


