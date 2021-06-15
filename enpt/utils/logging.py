# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz-potsdam.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz-potsdam.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz-potsdam.de)
#
# This software was developed within the context of the EnMAP project supported
# by the DLR Space Administration with funds of the German Federal Ministry of
# Economic Affairs and Energy (on the basis of a decision by the German Bundestag:
# 50 EE 1529) and contributions from DLR, GFZ and OHB System AG.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. Please note the following exception: `EnPT` depends on tqdm, which
# is distributed under the Mozilla Public Licence (MPL) v2.0 except for the files
# "tqdm/_tqdm.py", "setup.py", "README.rst", "MANIFEST.in" and ".gitignore".
# Details can be found here: https://github.com/tqdm/tqdm/blob/master/LICENCE.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""EnPT logging module containing logging related classes and functions."""

import logging
import os
import warnings
import sys
from io import StringIO  # Python 3 only

__author__ = 'Daniel Scheffler'


class EnPT_Logger(logging.Logger):
    """Class for the EnPT logger."""

    def __init__(self,
                 name_logfile: str,
                 fmt_suffix: any = None,
                 path_logfile: str = None,
                 log_level: any = 'INFO',
                 append: bool = True) -> None:
        """Return a logging.logger instance pointing to the given logfile path.

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

        # attributes that need to be present in order to unpickle the logger via __setstate_
        self.name_logfile = name_logfile
        self.fmt_suffix = fmt_suffix
        self.path_logfile = path_logfile
        self.log_level = log_level

        super(EnPT_Logger, self).__init__(name_logfile)

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

        # create StreamHandler
        self.streamObj = StringIO()
        streamHandler = logging.StreamHandler(stream=self.streamObj)
        streamHandler.setFormatter(self.formatter_fileH)
        streamHandler.set_name('StringIO handler')

        # create ConsoleHandler for logging levels DEGUG and INFO -> logging to sys.stdout
        consoleHandler_out = logging.StreamHandler(stream=sys.stdout)  # by default it would go to sys.stderr
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
            self.addHandler(streamHandler)
            self.addHandler(consoleHandler_out)
            self.addHandler(consoleHandler_err)

    def __getstate__(self):
        self.close()
        return self.__dict__

    def __setstate__(self, ObjDict):
        """Define how the attributes of EnPT_Logger are unpickled."""
        self.__init__(ObjDict['name_logfile'], fmt_suffix=ObjDict['fmt_suffix'], path_logfile=ObjDict['path_logfile'],
                      log_level=ObjDict['log_level'], append=True)
        ObjDict = self.__dict__
        return ObjDict

    @property
    def captured_stream(self) -> str:
        """Return the already captured logging stream.

        NOTE:
            - set self.captured_stream:
                self.captured_stream = 'any string'
        """
        if not self._captured_stream:
            self._captured_stream = self.streamObj.getvalue()

        return self._captured_stream

    @captured_stream.setter
    def captured_stream(self, string: str):
        assert isinstance(string, str), "'captured_stream' can only be set to a string. Got %s." % type(string)
        self._captured_stream = string

    def close(self):
        """Close all logging handlers."""
        # update captured_stream and flush stream
        self.captured_stream += self.streamObj.getvalue()

        for handler in self.handlers[:]:
            try:
                if handler.get_name() == 'StringIO handler':
                    self.streamObj.flush()
                self.removeHandler(handler)  # if not called with '[:]' the StreamHandlers are left open
                handler.flush()
                handler.close()
            except PermissionError:
                warnings.warn('Could not properly close logfile due to a PermissionError: %s' % sys.exc_info()[1])
            except ValueError as e:
                if str(e) != 'I/O operation on closed file':
                    raise

        if self.handlers[:]:
            warnings.warn('Not all logging handlers could be closed. Remaining handlers: %s' % self.handlers[:])

    def view_logfile(self):
        """View the log file written to disk."""
        with open(self.path_logfile) as inF:
            print(inF.read())

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return True if exc_type is None else False


def close_logger(logger):
    """Close the handlers of the given logging.Logger instance.

    :param logger:  logging.Logger instance or subclass instance
    """
    if logger and hasattr(logger, 'handlers'):
        for handler in logger.handlers[:]:  # if not called with '[:]' the StreamHandlers are left open
            try:
                logger.removeHandler(handler)
                handler.flush()
                handler.close()
            except PermissionError:
                warnings.warn('Could not properly close logfile due to a PermissionError: %s' % sys.exc_info()[1])

        if logger.handlers[:]:
            warnings.warn('Not all logging handlers could be closed. Remaining handlers: %s' % logger.handlers[:])


def shutdown_loggers():
    """Shutdown any currently active loggers."""
    logging.shutdown()


class LessThanFilter(logging.Filter):
    """Filter class to filter log messages by a maximum log level.

    Based on http://stackoverflow.com/questions/2302315/
        how-can-info-and-debug-logging-message-be-sent-to-stdout-and-higher-level-messag
    """

    def __init__(self, exclusive_maximum, name=""):
        """Get an instance of LessThanFilter.

        :param exclusive_maximum:  maximum log level, e.g., logger.WARNING
        :param name:
        """
        super(LessThanFilter, self).__init__(name)
        self.max_level = exclusive_maximum

    def filter(self, record):
        """Filter funtion.

        NOTE: Returns True if logging level of the given record is below the maximum log level.

        :param record:
        :return: bool
        """
        # non-zero return means we log this message
        return True if record.levelno < self.max_level else False
