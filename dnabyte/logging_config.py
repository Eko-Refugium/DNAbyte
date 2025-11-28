# logging_config.py
import logging
import os
from datetime import datetime

# Set up logging configuration
current_time = datetime.now().strftime('%Y%m%d_%H%M%S')
log_filename = os.path.join('app','uploads', f'{current_time}_job.log')
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
    logging.FileHandler(log_filename),
    logging.StreamHandler()
])
logger = logging.getLogger(__name__)