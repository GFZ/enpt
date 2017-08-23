# -*- coding: utf-8 -*-

import logging
import os
import warnings
import sys


class EnPT_logger(logging.Logger):
    def __init__(self, name_logfile, fmt_suffix=None, path_logfile=None, log_level='INFO', append=True):
        # type: (str, any, str, any, bool) -> None
        """Returns a logging.logger instance pointing to the given logfile path.
        :param name_logfile:
        :param fmt_suffix:      if given, it will be included into log formatter
        :param path_logfile:    if no path is given, only a StreamHandler is created
        :param log_level:       the logging level to be used (choices: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL';
                                default: 'INFO')
        :param append:          <bool> whether to append the log message to an existing logfile (1)
                                or to create a new logfile (0); default=1
        """

        # private attributes
        self._captured_stream = ''

        super(EnPT_logger, self).__init__(name_logfile)

        self.path_logfile = path_logfile
        self.formatter_fileH = logging.Formatter('%(asctime)s' + (' [%s]' % fmt_suffix if fmt_suffix else '') +
                                                 ' %(levelname)s:   %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
        self.formatter_ConsoleH = logging.Formatter('%(asctime)s' + (' [%s]' % fmt_suffix if fmt_suffix else '') +
                                                    ':   %(message)s', datefmt='%Y/%m/%d %H:%M:%S')

        if path_logfile:
            # create output directory
            while not os.path.isdir(os.path.dirname(path_logfile)):
                try:
                    os.makedirs(os.path.dirname(path_logfile))
                except OSError as e:
                    if e.errno != 17:
                        raise
                    else:
                        pass

            # create FileHandler
            fileHandler = logging.FileHandler(path_logfile, mode='a' if append else 'w')
            fileHandler.setFormatter(self.formatter_fileH)
            fileHandler.setLevel(log_level)
        else:
            fileHandler = None

        # create StreamHandler # TODO add a StringIO handler
        # self.streamObj     = StringIO()
        # self.streamHandler = logging.StreamHandler(stream=self.streamObj)
        # self.streamHandler.setFormatter(formatter)
        # self.streamHandler.set_name('StringIO handler')

        # create ConsoleHandler for logging levels DEGUG and INFO -> logging to sys.stdout
        consoleHandler_out = logging.StreamHandler(stream=sys.stdout) # by default it would go to sys.stderr
        consoleHandler_out.setFormatter(self.formatter_ConsoleH)
        consoleHandler_out.set_name('console handler stdout')
        consoleHandler_out.setLevel(log_level)
        consoleHandler_out.addFilter(LessThanFilter(logging.WARNING))

        # create ConsoleHandler for logging levels WARNING, ERROR, CRITICAL -> logging to sys.stderr
        consoleHandler_err = logging.StreamHandler(stream=sys.stderr)
        consoleHandler_err.setFormatter(self.formatter_ConsoleH)
        consoleHandler_err.setLevel(logging.WARNING)
        consoleHandler_err.set_name('console handler stderr')

        self.setLevel(log_level)

        if not self.handlers:
            if fileHandler:
                self.addHandler(fileHandler)
            # self.addHandler(self.streamHandler)
            self.addHandler(consoleHandler_out)
            self.addHandler(consoleHandler_err)

        #     if append:
        #         logfileHandler = logging.FileHandler(path_logfile, mode='a')
        #     else:
        #         logfileHandler = logging.FileHandler(path_logfile, mode='w')
        #     logfileHandler.setFormatter(formatter)
        #     logfileHandler.setLevel(logging.CRITICAL)
        #     consoleHandler_out = logging.StreamHandler()
        #     consoleHandler_out.setFormatter(formatter)
        #     consoleHandler_out.setLevel(logging.CRITICAL)
        # #    logger.setLevel(logging.DEBUG)
        #     if CPUs == 1:
        #         if not logger.handlers:
        #             logger.addHandler(logfileHandler)
        #             logger.addHandler(consoleHandler_out)
        #     else:
        #         logger.addHandler(logfileHandler)
        #         logger.addHandler(consoleHandler_out)

    @property
    def captured_stream(self):
        if not self._captured_stream:
            self._captured_stream = self.streamObj.getvalue()

        return self._captured_stream

    @captured_stream.setter
    def captured_stream(self, string):
        assert isinstance(string, str), "'captured_stream' can only be set to a string. Got %s." %type(string)
        self._captured_stream = string

    def close(self):
        # update captured_stream and flush stream
        # self.captured_stream += self.streamObj.getvalue()
        # print(self.handlers[:])

        # self.removeHandler(self.streamHandler)
        # print(dir(self.streamHandler))
        # self.streamHandler = None

        for handler in self.handlers[:]:
            try:
                # if handler.get_name() == 'StringIO handler':
                #     self.streamObj.flush()
                #     self.streamHandler.flush()

                handler.close()
                self.removeHandler(handler) # if not called with '[:]' the StreamHandlers are left open
            except PermissionError:
                warnings.warn('Could not properly close logfile due to a PermissionError: %s' % sys.exc_info()[1])

        if self.handlers[:]:
            warnings.warn('Not all logging handlers could be closed. Remaining handlers: %s' % self.handlers[:])

        # print('sh', self.streamHandler)

    def view_logfile(self):
        with open(self.path_logfile) as inF:
            print(inF.read())


def close_logger(logger):
    if logger and hasattr(logger, 'handlers'):
        for handler in logger.handlers[:]:  # if not called with '[:]' the StreamHandlers are left open
            try:
                handler.close()
                logger.removeHandler(handler)
            except PermissionError:
                warnings.warn('Could not properly close logfile due to a PermissionError: %s' %sys.exc_info()[1])

        if logger.handlers[:]:
            warnings.warn('Not all logging handlers could be closed. Remaining handlers: %s' % logger.handlers[:])


def shutdown_loggers():
    logging.shutdown()


class LessThanFilter(logging.Filter):
    # http://stackoverflow.com/questions/2302315/how-can-info-and-debug-logging-message-be-sent-to-stdout-and-higher-level-messag
    def __init__(self, exclusive_maximum, name=""):
        super(LessThanFilter, self).__init__(name)
        self.max_level = exclusive_maximum

    def filter(self, record):
        # non-zero return means we log this message
        return True if record.levelno < self.max_level else False