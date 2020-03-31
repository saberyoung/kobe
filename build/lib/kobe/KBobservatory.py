#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : kobe/observatory.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

# Standard library
import logging

# Third-party
import astropy.coordinates
import astropy.units as u
import astropy.time
import numpy as np
from kobe import telescope

__all__ = ['observatory']

class observatory(telescope):
    """define an observaroty.    

    Parameters
    ----------
    location           : `class`
        astropy.coordinates.EarthLocation object
    strategy           : `str` 
        options: `T`, `G` 
        `T` generate a tiling network
        `G` generate a list of galaxies
    name               : `str`,
        observatory name, if None will construct with its location
    latitude           : `int`, `float`, `str`
        latitude of observatory
    longitude          : `int`, `float`, `str`
        longitude of observatory
    elevation          : `int`, `float`, `str`
        elevation of observatory   
    clobber            : `bool`
       if any products already available in self, clobber=True will redo it   
    """                
    defkwargs = {'location'           : None,                                
                 'obsname'            : None,
                 'latitude'           : None,
                 'longitude'          : None,
                 'elevation'          : 0*u.m,                                        
    }
    '''
    default parameter settings
    '''
    defkwargs = {**defkwargs, **telescope.defkwargs}
    
    location = None
    '''
    location of observaroty
    '''
    
    obsname = None
    '''
    name of observaroty
    '''
    
    telescopes = {}
    '''
    a series of telescopes located at the observaroty
    '''

    def __init__(self, logger=None, **kwargs):
        """
        initialize observatory object with input logger and parameter settings
        
        Parameters
        ----------
        logger :     `class`
            logging object

        **kwargs :     
            check kobe.tilings.defkwargs, kobe.galaxies.defkwargs
        """
        
        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger
            
        # ----- set keys ----- #        
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return                
        self.defkwargs = kwargs
        
        # ----- default runs----- #        
        self._check_iers()                

    def check_obs(self):
        assert not self.obsname is None, 'no observatory name parsed'
        if len(self.telescopes) == 0:
            self.logger.info('Warning: no telescope available in current observatory')
        
    def set_obsname(self, *args, **kwargs):
        """set observaroty name

        Parameters
        ----------
        name :      `string` 
                  if not telescope name given, construct with telescope coordinate:
                  lon-lat-alt
        """
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  

        if len(args) > 0:
            self.obsname = args[0]
            
        elif not kwargs['obsname'] is None:
            self.obsname = kwargs['obsname']

        elif kwargs['location'] is None and (kwargs['latitude'] is \
                                             not None and kwargs['longitude'] is not None):
            self.obsname = 'T%.2f-%.2f' % (kwargs['longitude'], kwargs['latitude'])
            
        elif isinstance(kwargs['location'], astropy.coordinates.EarthLocation):
            self.obsname = 'T%.2f-%.2f' % (kwargs['location'].lon.deg,\
                                           kwargs['location'].lat.deg)
            
        else:
            self.logger.info ('Warning: Observatory name set failed, '+\
                              'give either an obsname or a location')
            return
        
    def _check_iers(self, iersauto=False,clear_iers=False,
                    iers_timeout=20, **kwargs):
        '''check astropy iers settings

        Parameters
        ----------
        iersauto           : `bool`
        if False, turn off iers auto downloading
        clear_iers         : `bool`
        if True, clear iers cached files
        iers_timeout       : `int`, `float`
        maximum time for iers query timeout
        '''        
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  
        
        # check iers
        from astropy.utils import iers
        iers.remote_timeout = iers_timeout
        
        if not iersauto:
            # turn off downloading            
            iers.conf.auto_download = False
            
        if clear_iers:
            from astropy.utils.data import clear_download_cache
            clear_download_cache()

    def add_telescope(self, *args, **kwargs):
        
        if len(args) > 0:
            _telname = args[0]
            self.set_telname(_telname)
        elif not self.telname is None:
            _telname = self.telname
        else:           
            self.logger.info('Warning: telname not set')
            return
        kwargs = self.setkeys(kwargs)
        self.check_tel()       
        if _telname in self.telescopes and not kwargs['clobber']:
            self.logger.info ('telescope %s already added'%_telname)
        else:
            self.telescopes[_telname] = self.pointings            
            self.logger.info ('telescope %s added'%_telname)        
            
    def del_telescope(self, telname):
                
        assert type(self.telescopes) is dict
        del self.telescopes[telname]

    def to_telescope(self, telname):        
        assert telname in self.telescopes
        self.telname = telname        
        self.pointings = self.telescopes[telname]

    def update_telescope(self, telname):        
        assert telname in self.telescopes
        self.telescopes[telname] = self.pointings
        
    def parse_location(self, *args, **kwargs):

        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  

        if len(args) > 0: self.obsname = args[0]
                
        if not self.obsname is None:
            _name = self.obsname

        elif not kwargs['obsname'] is None:
            _name = kwargs['obsname']

        else:
            _name = None
        
        if not _name is None:
            _loc = self.EarthLocation_from_kobe(_name)            
            if _loc is None:
                self.location = astropy.coordinates.EarthLocation.of_site(_name)
            else:
                self.location = _loc
                    
        elif kwargs['location'] is None and (kwargs['latitude'] is \
                            not None and kwargs['longitude'] is not None):
            self.location = astropy.coordinates.EarthLocation.from_geodetic(\
                            kwargs['longitude'], kwargs['latitude'], kwargs['elevation'])
                
        elif isinstance(kwargs['location'], astropy.coordinates.EarthLocation):
            self.location = kwargs['location']

        else:
            self.logger.info ('Warning: Observatory location load failed')        

        self.set_obsname(*args, **kwargs)
        
    @classmethod
    def EarthLocation_from_kobe(self, string):
        
        EarthLocation = astropy.coordinates.EarthLocation
        
        subaru = EarthLocation.from_geodetic(-155.4761111111111*u.deg,
                                             19.825555555555564*u.deg,
                                             4139*u.m)

        lco = EarthLocation.from_geodetic(-70.70166666666665*u.deg,
                                          -29.003333333333327*u.deg,
                                          2282*u.m)

        aao = EarthLocation.from_geodetic(149.06608611111113*u.deg,
                                          -31.277038888888896*u.deg,
                                          1164*u.m)

        vbo = EarthLocation.from_geodetic(78.8266*u.deg,
                                          12.576659999999999*u.deg,
                                          725*u.m)

        apo = EarthLocation.from_geodetic(-105.82*u.deg,
                                          32.78*u.deg,
                                          2798*u.m)

        keck = EarthLocation.from_geodetic(-155.47833333333332*u.deg,
                                           19.828333333333326*u.deg,
                                           4160*u.m)

        kpno = EarthLocation.from_geodetic(-111.6*u.deg,
                                           31.963333333333342*u.deg,
                                           2120*u.m)

        lapalma = EarthLocation.from_geodetic(-17.879999*u.deg,
                                              28.758333*u.deg,
                                              2327*u.m)

        asiago = EarthLocation.from_geodetic(11.5096*u.deg,
                                             45.8759*u.deg,
                                             3284*u.m)
        
        observatories = dict(lco=lco, subaru=subaru, aao=aao, vbo=vbo, apo=apo,
                             keck=keck, kpno=kpno, lapalma=lapalma, asiago=asiago)

        if string.lower() in observatories:
            return observatories[string.lower()]
        else:
            return
