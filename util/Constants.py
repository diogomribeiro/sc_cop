
import os

PATH_LOG = os.path.expanduser('~') + '/.scCOP.log'

#===============================================================================
# Constants for Logging 
#===============================================================================
MODE_DEBUG = "debug"
MODE_INFO = "info"
MODE_WARNING = "warning"
MODE_ERROR = "error"
MODE_CRITICAL = "critical"
MODE_FATAL = "fatal"
VERBOSITY_LEVELS = [ MODE_DEBUG, MODE_INFO, MODE_WARNING, MODE_ERROR, MODE_CRITICAL, MODE_FATAL]

LOG_APPEND = 'a'
LOG_NO_APPEND = 'w'

