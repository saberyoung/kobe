#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/PstPrioritization.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
import numpy as np
import healpy as hp
import astropy.time
import astropy.coordinates        
import astropy.units as u
        
__all__ = ('PstPrioritization')

class PstPrioritization():
    """PstPrioritization: optimize searching strategy considering probability and visibility

    * implement the following functions:
    -> __init__()
    -> run()

    Parameters
    ----------    

    Examples
    --------
    see also https://github.com/saberyoung/pstools/blob/master/notebook/test_prio.ipynb

    >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
    """    

    # Static version info
    version = 1.0
        
    def __init__(self, hpmap=None, fithdr=None,
                 xmlinf=None, defconf=None, logger=None):
        """
        initialize
        """

        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # ----- set trigger properties ----- #
        self.data               =   {}
        self.data['hpmap']      =   hpmap
        self.data['xmlinf']     =   xmlinf
        self.data['fithdr']     =   fithdr

        # ----- define keys to check xml ----- #
        '''
        skymap_key  `string`   key for skymap url, in XML   
        radec_keys  `list`     key for coordinate and uncertainty, in XML
        obst_key    `string`   key for timeobs, in XML
        '''
        self.obst_key           =   'ISOTime'
        self.skymap_key         =   'skymap_fits'
        self.radec_keys         =   ['C1','C2','Error2Radius']
        
        # ----- define informations for report ----- #
        self.keys_to_checkxml = ('GraceID', 'AlertType', 'Group', 'FAR',
                                 'Terrestrial', 'HasNS', 'HasRemnant',
                                 'BNS', 'BBH', 'NSBH', 'Instruments',
                                 'EventPage')
        self.keys_to_checkfits = ('DISTMEAN', 'DISTSTD', 'DATE-OBS',
                                  'MJD-OBS', 'OBJECT', 'INSTRUME',
                                  'CREATOR')

        # ----- define default parameters ----- #
        self.run_config(defconf)

    def run_config(self, defconf):
        """read default value for parameters

        Parameters
        ----------
        defconf :      `dict`
          input parameter settings, if None, use default

        Examples
        --------    
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> dict = {'cls': [.5,.9,.99]}
        >>> a.run_config(dict)
        """
        
        self.conf = {
            'wdir':      '/tmp/',   # `str`   working directory
            'savefits':  None,      # `str`   save file or not
            'nside':     512,       # `int`   healpix map default resolution
            'coord':     'C',       # `str`   healpix map defualt coordinate system: G, C, E
            'nest':      False,     # `bool`  healpix map defualt ordering: nest or ring
            'cls':       [.5,.9],   # `list`  confidence levels for contours
            'style':     'sms',     # `str`   report type           
            'lat':       24.625,    # `float`   la palma chile
            'lon':       70.403,    # `float`   la palma chile
            'alt':       2635,      # `float`   la palma chile
            'obstime':   None
        }
        
        if defconf is None: return
        for _k in self.conf.keys():
            if _k in defconf.keys():
                self.conf[_k] = defconf[_k]
            else:
                self.logger.info ('### Warning: use default value for %s'%_k)
    
    def run(self, url, wdir=None, savefits=None, nest=None):
        """parser skymap url, and build PstParseTriggers.data  

        Parameters
        ----------
        url :        `string`         
          url of healpix fits
        savefits :   `string` or None    
          save healpix fits or not,
          if None, not save; 
          if a string, savefits would be used to construct the downloaded file.
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.
        nest :       `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        Examples
        --------       
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        
        # parse fits_url
        if not url is None: self.fits(url,nest=nest)

        # savefits
        if not self.conf['savefits'] is None:
            flag = self.download_skymap(self.conf['wdir'], self.conf['savefits'], url)
            if flag == 1:
                self.logger.info ('download skymap via wget')
            elif flag == 2:
                self.logger.info ('download skymap via requests')
            elif flag == 3:
                self.logger.info ('### Warning: failed to download skymap,'+\
                                  'install wget or requests first')
           
    
    @staticmethod
    def obstime(t=None):
        """ define observing time

        Parameters
        ----------   
        t :   `float` or None or astropy.time
          observing time.
          None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; `astropy.time` set specific observing time

        returns
        ----------   
        obstime :   `astropy.time`
        """                
        if t is None:
            obstime = astropy.time.Time.now()                 
        elif type(t) in [int,float]:
            obstime = astropy.time.Time.now() + \
                astropy.time.TimeDelta(t, format='sec')            
        elif type(t) == str:
            try:
                obstime = astropy.time.Time(t, scale='utc')
            except:
                obstime = astropy.time.Time.now()                
        return obstime
