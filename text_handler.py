import logging
from difflib import SequenceMatcher

logger = logging.getLogger('main')

def read_table_to_dict(file_to_read, skip_rows=0, delimiter ='\t', key_type=str, value_type=str, secondary_delimiter=';', secondary_value_type=float, values_to_keys=False, open_operator=open):
    '''
    parse a delimited file with two columns to a dictionary
    if values_to_keys == 0 keys are taken from the first column o.w., from the second
    open_operator can also be set to gzip.open if the file is zipped
    '''
    logger.info('Reading dict from: ' + file_to_read)
    d = {}
    with open_operator(file_to_read) as f:
        for i in range(skip_rows):
            f.readline()
        for line in f:
            if line != '' and not line.isspace():
                key, value = line.rstrip().split(delimiter)
                if value_type == list:
                    value = [secondary_value_type(item) for item in value.split(secondary_delimiter)]
                d[key_type(key)] = value_type(value)
    if values_to_keys:
        d = dict((d[key],key) for key in d)
    return d


def write_dict_to_file(out_path, d, delimiter='\t', value_type=str, secondary_delimiter=';', reverse=False, open_operator=open, sort_by=None, header=''):
    logger.info('Writing dict to: ' + out_path)
    if d:
        with open_operator(out_path, 'w') as f:
            if header:
                f.write(header + '\n')
            for key in sorted(d, key=sort_by, reverse=reverse):
                if value_type == list:
                    value = secondary_delimiter.join(str(item) for item in d[key])
                else:
                    value = str(d[key])
                f.write(delimiter.join([str(key), value]) + '\n')
    else:
        logger.info('Skipping. Dictionary is empty...')


def get_next_relevant_rstriped_line(file_handler):
    while True:
        line = file_handler.readline()
        if line == '': #EOF
            return None
        line = line.strip()
        if line == '' or line[0] == '#':
            logger.debug('Skipping irrelevant (i.e., comment/white spaces) line:\n' + line)
        else:
            logger.info('Parsing next relevant line:\n' + line)
            return line


def string_similarity(a, b):
    # measure similarity rate of two strings
    return SequenceMatcher(None, a, b).ratio()


