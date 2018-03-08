
import os, logging
from sys import argv
logger = logging.getLogger('main')

def create_dir(path):
    if not os.path.exists(path):
        logger.info("Creating directory:\n" + path)
        os.makedirs(path)

if __name__ == '__main__':
    if len(argv) < 2:
        logger.error('Usage: python ' + argv[0] + ' <path>')
        exit()
    else:
        create_dir(*argv[1:])