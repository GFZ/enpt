# -*- coding: utf-8 -*-

from geomultisens.misc.logging import GMS_logger


# modify the GeoMultiSens logger fit the needs of EnPT
class EnPT_logger(GMS_logger):

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

        super(EnPT_logger, self).__init__(name_logfile,
                                          fmt_suffix=fmt_suffix,
                                          path_logfile=path_logfile,
                                          log_level=log_level,
                                          append=append)

    # TODO add origin of errors module.module.submodule to logging format
