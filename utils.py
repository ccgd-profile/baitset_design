#! /usr/bin/python

import os
import sys


def setup_logger(logFnPath, name):
    """Creates the logger object and associated text file to use throughout
    the analysis.
    It first creates a log.txt file in the specified analyis directory as the
    FileHandler. The console handler is then formatted to report the time,
    name of the source, level of the message and the message.
    Args:
        log_fn_path: Absolute path for the directory that will contain the
        log file.
        name: The name of the package to initially setup the logger object.
    Returns:
        Nothing is returned.
    """

    outputPath = os.path.abspath(logFnPath)
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # FileHandler
    fileHandle = logging.FileHandler(os.path.join(outputPath, 'log.txt'), mode='w')
    fileHandle.setLevel(logging.DEBUG)
    # ConsoleHandler
    consoleHandle = logging.StreamHandler()
    consoleHandle.setLevel(logging.ERROR)

    formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fileHandle.setFormatter(formatter)
    consoleHandle.setFormatter(formatter)

    logger.addHandler(fileHandle)
    logger.addHandler(consoleHandle)