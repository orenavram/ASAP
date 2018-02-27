import os, logging
def create_dir(path, logger = logging.getLogger('main')):
    if not os.path.exists(path):
        logger.info("Creating directory:\n" + path)
        os.makedirs(path)