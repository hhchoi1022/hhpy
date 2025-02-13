#%%
import glob
import inspect
import json
import os
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats, sigma_clip
import numpy as np
from astropy.table import Table
from astropy.time import Time
from Research.photometry import Catalog
from Research.helper import Helper
import logging
import logging.handlers
import datetime
import os
from astropy.time import Time
import  astropy.units as u
#%%
class ProcessStatus:
    """ Manages the image processing steps status """
    def __init__(self, biascor=False, darkcor=False, flatcor=False,
                 astrometrycalc=False, scampcalc=False, zpcalc=False, bkgsub=False, 
                 combine=False, align=False, subtract= False, photometry=False):
        self.biascor = dict(status=biascor, update_time=None)
        self.darkcor = dict(status=darkcor, update_time=None)
        self.flatcor = dict(status=flatcor, update_time=None)
        self.astrometrycalc = dict(status=astrometrycalc, update_time=None)
        self.scampcalc = dict(status=scampcalc, update_time=None)
        self.zpcalc = dict(status=zpcalc, update_time=None)
        self.bkgsub = dict(status=bkgsub, update_time=None)
        self.combine = dict(status=combine, update_time=None)
        self.align = dict(status=align, update_time=None)
        self.subtract = dict(status=subtract, update_time=None)
        self.photometry = dict(status=photometry, update_time=None)

    def __repr__(self):
        """ Represent process status as a readable string """
        status_list = [f"{key}: {value['status']} (Updated: {value['update_time']})"
                       for key, value in self.__dict__.items()]
        return "ProcessStatus =====================================\n  " + "\n  ".join(status_list) + "\n==================================================="

    def to_dict(self):
        """ Convert class instance to dictionary """
        return self.__dict__

    @classmethod
    def from_dict(cls, data):
        """ Load ProcessStatus from dictionary """
        instance = cls()
        for key, value in data.items():
            setattr(instance, key, value)
        return instance
    
    
class Logger:
    """
    A class for creating and managing loggers.

    Parameters
    ----------
    unitnum : int
        The unit number.
    logger_name : str
        The name of the logger.
    **kwargs : dict, optional
        Additional keyword arguments.

    Methods
    -------
    log()
        Get the logger instance.
    createlogger(logger_name)
        Create a logger instance.
    """
    def __init__(self,
                 logger_name):
        self._log = self.createlogger(logger_name)
        self.path = logger_name
    
    def log(self):
        """
        Get the logger instance.

        Returns
        -------
        logging.Logger
            The logger instance.
        """
        return self._log
    
    def createlogger(self,
                     logger_name,
                     logger_level = 'INFO'):
        """
        Create a logger instance.

        Parameters
        ----------
        logger_name : str
            The name of the logger.

        Returns
        -------
        logging.Logger
            The created logger instance.
        """
        # Create Logger
        logger = logging.getLogger(logger_name)
        # Check handler exists
        if len(logger.handlers) > 0:
            return logger # Logger already exists
        logger.setLevel(logger_level)
        formatter = logging.Formatter(datefmt = '%Y-%m-%d %H:%M:%S',fmt = f'[%(levelname)s] %(asctime)-15s] | %(message)s')
        
        # Create Handlers
        streamHandler = logging.StreamHandler()
        streamHandler.setLevel(logger_level)
        streamHandler.setFormatter(formatter)
        logger.addHandler(streamHandler)
        fileHandler = logging.FileHandler(filename = logger_name)
        fileHandler.setLevel(logger_level)
        fileHandler.setFormatter(formatter)
        logger.addHandler(fileHandler)
        return logger
#%%
class CalibrationImage:
    """ Handles FITS image processing and tracks its status """

    def __init__(self, path: str, telinfo):
        self.helper = Helper()
        self.path = path
        dirname = os.path.dirname(path)
        filename = os.path.basename(path)
        self.statuspath = os.path.join(dirname, filename.split('.fits')[0] + '.status')
        self.telinfo = telinfo
        self._data = None
        self._header = None
        logger_name = os.path.join(dirname, filename.split('.fits')[0] + '.log')
        self.logger = Logger(logger_name = logger_name).log()
        self.loggerpath = logger_name
        
        # Initialize or load status
        if os.path.exists(self.statuspath):
            self.status = self.load_from_status()
        else:
            self.status = ProcessStatus()
            self.save_to_status()
    
    def load_from_status(self):
        """ Load processing status from a JSON file """
        with open(self.statuspath, 'r') as f:
            status_data = json.load(f)
        return ProcessStatus.from_dict(status_data)

    def save_to_status(self):
        """ Save processing status to a JSON file """
        with open(self.statuspath, 'w') as f:
            json.dump(self.status.to_dict(), f, indent=4)

    def update_status(self, process_name):
        """ Mark a process as completed and update time """
        if hasattr(self.status, process_name):
            self.status.__dict__[process_name]['status'] = True
            self.status.__dict__[process_name]['update_time'] = Time.now().iso
            self.save_to_status()
        else:
            raise ValueError(f"Invalid process name: {process_name}")

    def __repr__(self):
        return f"CalibrationImage(path = {os.path.basename(self.path)})"

    @property
    def data(self):
        """ Lazy loading of FITS data """
        if self._data is None:
            self._data = fits.getdata(self.path)
        return self._data
    
    @property
    def header(self):
        """ Lazy loading of FITS header """
        if self._header is None:
            self._header = fits.getheader(self.path)
        return self._header
