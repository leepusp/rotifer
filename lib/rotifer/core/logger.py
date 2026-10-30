import os
import sys
import logging
from tqdm import tqdm

# Rotifer's root logger handler class
# https://stackoverflow.com/questions/38543506/change-logging-print-function-to-tqdm-write-so-logging-doesnt-interfere-wit
class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            # stderr: stdout carries the data of Rotifer's command line tools
            tqdm.write(msg, file=sys.stderr)
            self.flush()
        except Exception:
            self.handleError(record)

# Create root logger and default handlers
def get_logger(name):
    '''
    Return the logger for name, making sure messages reach one handler.

    The handler is attached once, to Rotifer's root logger ("rotifer"),
    and loggers under it (e.g. "rotifer.core.functions") reach it by
    propagation. Loggers outside Rotifer's namespace, such as a script's
    "__main__", get a handler of their own, also only once.
    '''
    owner = 'rotifer' if name == 'rotifer' or name.startswith('rotifer.') else name
    target = logging.getLogger(owner)
    if not any(isinstance(h, TqdmLoggingHandler) for h in target.handlers):
        RotiferDefaultFormatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
        RotiferDefaultHandler = TqdmLoggingHandler()
        RotiferDefaultHandler.setFormatter(RotiferDefaultFormatter)
        target.addHandler(RotiferDefaultHandler)
    return logging.getLogger(name)
