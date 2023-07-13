import sys
import logging

logging.basicConfig(level=logging.DEBUG,
                    filename='output.log',
                    datefmt='%Y/%m/%d %H:%M:%S',
                    #format='%(asctime)s - %(name)s - %(levelname)s - %(lineno)d - %(module)s - %(message)s')
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def lexit(message, code=1):
    """Log an error message and exit."""
    logger.error(message)
    sys.exit(code)

def write2log(message):
    """Log every action and record it."""
    logger.info(message)