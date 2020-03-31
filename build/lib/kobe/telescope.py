#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/telescope.py
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
from kobe import cookbook

__all__ = ('telescope')

class telescope(object):
    
    # Static version info
    version = 1.0

    # default key args
    defkwargs = {'location'           : None,
                 'time'               : None,
                 'timeu'              : 'sec',
                 'name'               : None,
                 'latitude'           : None,
                 'longitude'          : None,
                 'elevation'          : 0*u.m,
                 'solarobj'           : None,
                 'pressure'           : 0.0,
                 'temperature'        : 100.0,
                 'relative_humidity'  : 0.0,
                 'obswl'              : 1.0,
                 'iers'               : True,
                 'clear_iers'         : False,
                 'iers_timeout'       : 10,
                 'clobber'            : True}
    
    def __init__(self, logger=None, **kwargs):
        """
        initialize
        """

        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger
           
        # ----- set keys ----- #        
        for _key in self.defkwargs:                        
            kwargs.setdefault(_key, self.defkwargs[_key])
        self.kwargs = kwargs
        
        # ----- default runs----- #                
        if not 'location' in self.__dict__: self.location = None
        if not 'obstime' in self.__dict__: self.obstime = None
        self._check()
        
    def _check(self, **kwargs):        
        self._check_iers(**kwargs)
        self._check_loc(**kwargs)
        self._check_time(**kwargs)

    def _check_iers(self, **kwargs):

        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

        from astropy.utils import iers
        iers.remote_timeout = kwargs['iers_timeout']
        
        if not kwargs['iers']:
            # turn off downloading            
            iers.conf.auto_download = False
            
        if kwargs['clear_iers']:
            from astropy.utils.data import clear_download_cache
            clear_download_cache()  


    def _check_loc(self, **kwargs):
                
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        if not self.location is None and not kwargs['clobber']: return
        
        if not kwargs['name'] is None:
            self.location = astropy.coordinates.EarthLocation.of_site(kwargs['name'])
        
        elif kwargs['location'] is None and (kwargs['latitude'] is \
                                    not None and kwargs['longitude'] is not None):
            self.location = astropy.coordinates.EarthLocation.from_geodetic(kwargs['longitude'], \
                                    kwargs['latitude'], kwargs['elevation'])
            
        elif isinstance(kwargs['location'], astropy.coordinates.EarthLocation):
            self.location = kwargs['location']

        else:
            self.logger.info ('Observatory location must be specified with '
                              'either (1) an instance of '
                              'astropy.coordinates.EarthLocation or (2) '
                              'latitude and longitude in degrees as '
                              'accepted by astropy.coordinates.Latitude and '
                              'astropy.coordinates.Latitude.')
            return

    def _check_time(self, **kwargs):
        """ define observing time

        Parameters
        ----------   
        time :   astropy.time, `float`, `int`, None
          observing time.
          None for now; 
          `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; 
          `astropy.time` set specific observing time

        returns
        ----------   
        obstime :   `astropy.time`
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not self.obstime is None and not kwargs['clobber']: return        
        if kwargs['time'] is None: return
        
        if kwargs['time'] == 'now':
            self.obstime = astropy.time.Time.now()
            
        elif type(kwargs['time']) in [int,float]:
            self.obstime = astropy.time.Time.now() + astropy.time.TimeDelta(kwargs['time'], format=kwargs['timeu'])

        elif isinstance(kwargs['time'], astropy.time.Time):
            self.obstime = kwargs['time']
        
        else:
            try:
                self.obstime = astropy.time.Time(t, scale='utc')
            except:
                self.logger.info ('time keyword should be a string, or an '
                                  'instance of astropy.time.Time')

    def altaz(self, ra=None, dec=None, theta=None, phi=None, **kwargs):
        """calculate altazimuth

        Parameters
        ----------   
        lon :         `float` or None
          longtitude of telescope, defaul: la palma chile
        lat :         `float` or None
          latitude of telescope, defaul: la palma chile
        alt :         `float` or None
          altitude of telescope, defaul: la palma chile
        obstime :     `float` or None or astropy.time obj
          observing time.
          None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; `astropy.time` set specific observing time
        cls :         `list`
          list of confidence level, default: [.5, .9]
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        returns
        ----------
        altaz :         `list`
          altaz list

        obstime :      `astropy.time`
          observing time, 

        observaroty :  `astropy.coordinate`
          observaroty location
        """

        self._check(**kwargs)
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.location is None or self.obstime is None:            
            self.logger.info ('Error: wrong telescope location and obstime setting')
            return
        
        if not ra is None and not dec is None:            
            if cookbook.is_seq(ra) and cookbook.is_seq(dec):
                ra, dec = np.array(ra), np.array(dec)
            
        elif not theta is None and not phi is None:
            if cookbook.is_seq(theta) and cookbook.is_seq(phi):
                theta, phi = np.array(thata), np.array(phi)
            dec, ra = -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
            
        else:
            self.logger.info ('Error: wrong altaz input')
            return
        
        radecs = astropy.coordinates.SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        frame = astropy.coordinates.AltAz(obstime=self.obstime, location=self.location,
                    pressure=kwargs['pressure']*u.hPa, temperature=kwargs['temperature']*u.deg_C,
                    relative_humidity=kwargs['relative_humidity'], obswl=kwargs['obswl']*u.micron)        
        altaz = radecs.transform_to(frame)        
        return altaz
    
    def calc_airmass(self, **kwargs):
        
        """calculate airmass for tilings
        
        Parameters
        ----------   
            lon :         `float` or None
                longtitude of telescope, defaul: la palma chile
            lat :         `float` or None
                latitude of telescope, defaul: la palma chile
            alt :         `float` or None
                altitude of telescope, defaul: la palma chile
            obstime :     `float` or None or astropy.time obj
                observing time.
                None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; `astropy.time` set specific observing time               

        returns
        ----------
            table :         `astropy.tab`
                  aismass of each tilings

        Examples
        -------- 
        >>> from kobe.pipeline.KBGetTilings import KBGetTilings
        >>> a = KBGetTilings()
        >>> a.generate(limdec=[-20,90]) 
        >>> b=a.calc_airmass()
        >>> b      
        """
        altaz = self.altaz(**kwargs)
        if not altaz is None: return altaz.secz
            
    def calc_solar(self, ra=None, dec=None, theta=None, phi=None, **kwargs):
        """calculate distance of tilings to solar object

        Parameters
        ----------   
        sobj :        `string`
                options:  `earth`, `sun`, `earth-moon-barycenter`, `moon`, `mercury`, 
                `venus`, `mars`, `jupiter`, `saturn`, `uranus`, `neptune`
        lon :         `float` or None
                longtitude of telescope, defaul: la palma chile
        lat :         `float` or None
                latitude of telescope, defaul: la palma chile
        alt :         `float` or None
                altitude of telescope, defaul: la palma chile
        obstime :     `float` or None or astropy.time obj
                observing time.
                None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; `astropy.time` set specific observing time               

        returns
        ----------
        table :         `astropy.tab`
                distance of each galaxies to selected solar object, unit in deg

        Examples
        -------- 
        >>> from kobe.pipeline.KBGetTilings import KBGetTilings
        >>> a = KBGetTilings()
        >>> a.generate(limdec=[-20,90]) 
        >>> b=a.calc_solar('moon')
        >>> b                
        """
        self._check(**kwargs)
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        _list = ['earth', 'sun', 'earth-moon-barycenter', 'moon',
                 'mercury','venus','mars','jupiter', 'saturn',
                 'uranus','neptune']
        if not kwargs['solarobj'] in _list:
            self.logger.info ('Error: solar obj select from %s' % _list)
            return
                
        altaz = self.altaz(ra=ra, dec=dec, the=theta, phi=phi, **kwargs)        
        if not ra is None and not dec is None:            
            if cookbook.is_seq(ra) and cookbook.is_seq(dec):
                ra, dec = np.array(ra), np.array(dec)
            
        elif not theta is None and not phi is None:
            if cookbook.is_seq(theta) and cookbook.is_seq(phi):
                theta, phi = np.array(thata), np.array(phi)
            dec, ra = -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
            
        else:
            self.logger.info ('Error: wrong altaz input')
            return
        
        radecs = astropy.coordinates.SkyCoord(ra=ra*u.deg, dec=dec*u.deg)             
        _sloc = astropy.coordinates.get_body(kwargs['solarobj'], self.obstime, self.location).transform_to(altaz)
        return radecs.separation(_sloc)

    def calc_sun(self, **kwargs):        

        """calculate sun height

        Parameters
        ----------   
            lon :         `float` or None
                  longtitude of telescope, defaul: la palma chile
            lat :         `float` or None
                  latitude of telescope, defaul: la palma chile
            alt :         `float` or None
                  altitude of telescope, defaul: la palma chile
            obstime :     `float` or None or astropy.time obj
                  observing time.
                  None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; `astropy.time` set specific observing time               

        returns
        ----------
        table :         `astropy`
                sun height

        Examples
        -------- 
        >>> from kobe.pipeline.KBGetTilings import KBGetTilings
        >>> a = KBGetTilings()
        >>> a.generate(limdec=[-20,90]) 
        >>> b=a.calc_sun()
        >>> b
        <Latitude -33.6547778 deg>               
        """
        
        altaz = self.altaz(**kwargs)
        if not altaz is None:
            sun_altaz = astropy.coordinates.get_sun(self.obstime).transform_to(altaz)
            return sun_altaz.alt

    def at_night(self):
        return
        
    def calc_vis(self, limalt=2.5, limsun=-18, obstime=None,
                 cls=None, nest=None, lat=None, lon=None, alt=None):               
        _hp = self.checkhp()
        if not _hp:
            self.logger.info ('### Warning: no hpmap found')
            return
        
        from kobe.cookbook import is_seq, is_seq_of_seq
        if is_seq_of_seq(self.data['hpmap']):
            (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        elif is_seq(self.data['hpmap']):
            hpx = self.data['hpmap']
        else:
            return

        _alist = {}
        
        if limalt is None:
            self.logger.info ('no airmass constrain')
            return       
             
        nside = hp.get_nside(hpx)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2        
        altazl, obstime, observatory = self.altaz(obstime=obstime,
                        cls=cls, nest=nest, lat=lat, lon=lon, alt=alt)               
        for cc in altazl:
            altaz = altazl[cc]
            sun_altaz = astropy.coordinates.get_sun(obstime).transform_to(altaz)            
            if limsun is None:
                self.logger.info ('no sun constrain')
                cond = np.logical_and(altaz.secz <= limalt, altaz.secz >= 1)
            else:               
                cond = np.logical_and(np.logical_and(altaz.secz <= limalt, altaz.secz >= 1),
                                      sun_altaz.alt <= limsun*u.deg)           
            _alist[cc] = len(altaz[cond])*areapix
        return _alist
    
def _generate_24hr_grid(t0, start, end, n_grid_points, for_deriv=False):
    """
    Generate a nearly linearly spaced grid of time durations.
    The midpoints of these grid points will span times from ``t0``+``start``
    to ``t0``+``end``, including the end points, which is useful when taking
    numerical derivatives.
    Parameters
    ----------
    t0 : `~astropy.time.Time`
        Time queried for, grid will be built from or up to this time.
    start : float
        Number of days before/after ``t0`` to start the grid.
    end : float
        Number of days before/after ``t0`` to end the grid.
    n_grid_points : int (optional)
        Number of grid points to generate
    for_deriv : bool
        Generate time series for taking numerical derivative (modify
        bounds)?
    Returns
    -------
    `~astropy.time.Time`
    """

    if for_deriv:
        time_grid = np.concatenate([[start - 1 / (n_grid_points - 1)],
                                    np.linspace(start, end,
                                                n_grid_points)[1:-1],
                                    [end + 1 / (n_grid_points - 1)]]) * u.day
    else:
        time_grid = np.linspace(start, end, n_grid_points) * u.day

    # broadcast so grid is first index, and remaining shape of t0
    # falls in later indices. e.g. if t0 is shape (10), time_grid
    # will be shape (N, 10). If t0 is shape (5, 2), time_grid is (N, 5, 2)
    while time_grid.ndim <= t0.ndim:
        time_grid = time_grid[:, np.newaxis]
    # we want to avoid 1D grids since we always want to broadcast against targets
    if time_grid.ndim == 1:
        time_grid = time_grid[:, np.newaxis]
    return t0 + time_grid
