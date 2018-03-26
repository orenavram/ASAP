import logging
logger = logging.getLogger('main')
from sys import argv

def main_func(path): #with a name that is related to the file's name
    pass # do something

def auxiliary_func1(path):
    pass # do something

def auxiliary_func2(path):
    pass # do something

#etc...

if __name__ == '__main__':
    if len(argv) < 4: # change the number of arguments according to your script
        logger.error('Usage: python ' + argv[0] + ' <arg1> <arg2> <arg3>...')
        exit()
    else:
        main_func(*argv[1:])