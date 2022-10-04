import logging


class TaskLogger(object):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL
    DEFAULT_FORMAT = "%(levelname)s  %(asctime)s  %(message)s"

    def __init__(self, log_params):
        self.log_params = log_params
        self.is_logging = log_params['logging']
        if not self.is_logging:
            return
        self.logger = logging.getLogger()
        self.mode = self.log_params['out']
        level = self.log_params['level']
        try:
            level = self.__getattribute__(level)
        except AttributeError:
            level = self.DEBUG

        log_format = self.DEFAULT_FORMAT

        if self.mode == 'file':
            file = self.log_params['log_file']
            logging.basicConfig(filename=file,
                                format=log_format,
                                encoding='utf-8',
                                level=level)
        elif self.mode == 'console':
            logging.basicConfig(encoding='utf-8',
                                format=log_format,
                                level=level)
        else:
            self.is_logging = False

    def log_event(self, msg_level, message):
        if not self.is_logging:
            return
        if msg_level == self.DEBUG:
            logging.debug(message)
        elif msg_level == self.INFO:
            logging.info(message)
        elif msg_level == self.ERROR:
            logging.error(message)
        elif msg_level == self.CRITICAL:
            logging.critical(message)
