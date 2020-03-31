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
from kobe import telescope

__all__ = ('observatory')

class observatory(telescope):
    """
    observaroty: define an observaroty.
    inherient from kobe.telescope

    Parameters
    ----------
    location           : `class`
        astropy.coordinates.EarthLocation object
    strategy           : `str` 
        options: `T`, `G` 
        `T` generate a tiling network
        `G` generate a list of galaxies
    time               : `class` or now or `str` or `float` or `int` or sequence
        define one or a list of observing time, so that for a calculation of 
                 airmass, sun height, etc
        options: 1. astropy.time.Time class
                 2. now: present current time
                 3. `str`, e.g. 1999-01-01T00:00:00.123456789
                 4. `float`, `int`: how much time before or later than currently
                    e.g. when timeu=sec, time=-3600 means 1 hour before
                 5. sequence of above options: will run tasks for a list of observing times
    timeu              : `str`
        unit of TimeDelta
        options: 
           jd:   float, long, decimal, str, bytes
           sec:  float, long, decimal, str, bytes
    name               : `str`,
        observatory name, if None will construct with its location
    latitude           : `int`, `float`, `str`
        latitude of observatory
    longitude          : `int`, `float`, `str`
        longitude of observatory
    elevation          : `int`, `float`, `str`
        elevation of observatory
    solarobj           : `str`
        name of solar object
        options: earth, sun, earth-moon-barycenter, moon,
                 mercury,venus,mars,jupiter, saturn,
                 uranus,neptune
    pressure           : `int`, `float`
        pressure in observatory, unit in hpa
    temperature        : `int`, `float`
        temperature in observatory, unit in C deg
    relative_humidity  : `int`, `float`
        humidity in observatory
    obswl : `~astropy.units.Quantity`
        Wavelength of the observation used in the calculation.
    iers               : `bool`
        if False, turn off iers auto downloading
    clear_iers         : `bool`
        if True, clear iers cached files
    iers_timeout       : `int`, `float`
        maximum time for iers query timeout
    clobber            : `bool`
       if any products already available in self, clobber=True will redo it

    See Also
    --------
    KBtilings, KBgalaxies

    Examples
    --------
    see https://github.com/saberyoung/kobe/blob/master/notebook/test_obs.ipynb
    """
        
    version = 1.0
    '''
    class kobe.KBplotter version
    '''
    
    defkwargs = {'location'           : None,                                
                 'obsname'            : None,
                 'latitude'           : None,
                 'longitude'          : None,
                 'elevation'          : 0*u.m,
                 'time'               : None,
                 'timeu'              : 'sec',     
                 'solarobj'           : None,
                 'pressure'           : 0.0,
                 'temperature'        : 100.0,
                 'relative_humidity'  : 0.0,
                 'obswl'              : 1.0,
                 'iers'               : False,
                 'clear_iers'         : False,
                 'iers_timeout'       : 20,                 
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
    location of observaroty
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
        self.parse_location(**kwargs)
        self.set_obsname(**kwargs)
        
    def set_obsname(self, **kwargs):
        """set observaroty name

        Parameters
        ----------
        name :      `string` 
                  if not telescope name given, construct with telescope coordinate:
                  lon-lat-alt
        """
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  
            
        if not kwargs['obsname'] is None:
            self.obsname = kwargs['obsname']

        elif kwargs['location'] is None and (kwargs['latitude'] is \
                                             not None and kwargs['longitude'] is not None):
            self.obsname = 'T%.2f-%.2f' % (kwargs['longitude'], kwargs['latitude'])
            
        elif isinstance(kwargs['location'], astropy.coordinates.EarthLocation):
            self.obsname = 'T%.2f-%.2f' % (kwargs['location'].lon.deg,\
                                           kwargs['location'].lat.deg)
            
        else:
            self.logger.info ('Warning: Observatory name set failed')
            return
        
    def _check_iers(self, **kwargs):

        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  
        
        # check iers
        from astropy.utils import iers
        iers.remote_timeout = kwargs['iers_timeout']
        
        if not kwargs['iers']:
            # turn off downloading            
            iers.conf.auto_download = False
            
        if kwargs['clear_iers']:
            from astropy.utils.data import clear_download_cache
            clear_download_cache()

    def add_telescope(self, **kwargs):

        kwargs = self.setkeys(kwargs)
        if kwargs is None: return        
        self.set_telname(**kwargs)
        
        for _name, _data in zip(['pointings', 'telname'],
                                [self.data, self.telname]):
            if _data is None:
                self.logger.info ('Warning: %s not parsed'%_name)
                return
        if self.telname in self.telescopes and not kwargs['clobber']:
            self.logger.info ('telescope %s already added'%self.telname)
        else:
            self.telescopes[self.telname] = self.data            
            self.logger.info ('telescope %s added'%self.telname)        
            
    def del_telescope(self, **kwargs):
        
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return
        self.set_telname(**kwargs)
        
        for _name, _data in zip(['pointings', 'telname'],
                                [self.data, self.telname]):
            if _data is None:
                self.logger.info ('Warning: %s not parsed'%_name)
                return            
        del self.telescopes[self.telname]
    
    def parse_location(self, **kwargs):

        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  
        
        if not self.location is None and not kwargs['clobber']:
            return
        
        if not kwargs['obsname'] is None:
            _loc = self.EarthLocation_from_kobe(kwargs['obsname'])
            if _loc is None:
                self.location = astropy.coordinates.EarthLocation.of_site(kwargs['obsname'])
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
            
    def frame(self, **kwargs):
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
        
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return
        
        if self.location is None:            
            self.logger.info ('Error: set telescope location first')
            return
        _times = self.parse_time(**kwargs)
        if _times is None:return
        frame = astropy.coordinates.AltAz(obstime=_times, location=self.location,
                    pressure=kwargs['pressure']*u.hPa, temperature=kwargs['temperature']*u.deg_C,
                    relative_humidity=kwargs['relative_humidity'], obswl=kwargs['obswl']*u.micron)
        return frame

    def parse_time(self, **kwargs):
        '''
        parse time for astronomical calculation
        '''
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return
        
        if kwargs['time'] is None:
            return
        
        elif self.is_seq(kwargs['time']) and not type(kwargs['time']) is str:
            _times = []
            for _time in kwargs['time']:
                _t = self.parse_time_single(_time)
                if not _t is None: _times.append(_t)
            return np.array(_times)

        else:
            return self.parse_time_single(kwargs['time'])

    def parse_time_single(self, _time, **kwargs):
        '''
        support for parse_time
        '''
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return
        
        if _time == 'now':
            return astropy.time.Time.now()
        
        elif isinstance(_time, astropy.time.Time):
            return _time

        try:            
            return astropy.time.Time.now() + \
                astropy.time.TimeDelta(float(_time), format=kwargs['timeu'])
        except:
            pass
        
        try:
            return astropy.time.Time(_time, scale='utc')
        except:
            self.logger.info ('Warning: time %s format wrong'%_time)
            
    def calc_sun_height(self, **kwargs):        

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
        
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return  
        
        frame = self.frame(**kwargs)
        _times = self.parse_time(**kwargs)
        if _times is None: return
                
        if isinstance(_times, astropy.time.Time):            
            sun_loc = astropy.coordinates.get_sun(_times)
            altaz = astropy.coordinates.AltAz(obstime=_times, location=sun_loc)     
            sun_altaz = sun_loc.transform_to(frame)
            sunalt = sun_altaz.alt
        else:
            sunalt = []            
            for _ii in range(len(_times)):                                 
                sun_loc = astropy.coordinates.get_sun(_times[_ii])
                altaz = astropy.coordinates.AltAz(obstime=_times[_ii], location=sun_loc)     
                sun_altaz = sun_loc.transform_to(frame[_ii])                
                sunalt.append(sun_altaz.alt)                
        return sunalt

    def calc_moon_illumination(self, ephemeris=None, **kwargs):        
        """
        Calculate fraction of the moon illuminated.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time of observation

        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            `~astropy.coordinates.solar_system_ephemeris` (which is
            set to 'builtin' by default).

        Returns
        -------
        k : float
            Fraction of moon illuminated
        """
        
        _times = self.parse_time(**kwargs)
        if _times is None: return

        # Calculate lunar orbital phase in radians. 
        if isinstance(_times, astropy.time.Time):            
            sun = astropy.coordinates.get_sun(_times)
            moon = astropy.coordinates.get_moon(_times, ephemeris=ephemeris)
            elongation = sun.separation(moon)
            i = np.arctan2(sun.distance*np.sin(elongation),
                           moon.distance - sun.distance*np.cos(elongation))
            k = (1 + np.cos(i))/2.0
            mill = k.value
    
        else:
            mill = []            
            for _ii in range(len(_times)):                
                sun = astropy.coordinates.get_sun(_times[_ii])
                moon = astropy.coordinates.get_moon(_times[_ii], ephemeris=ephemeris)
                elongation = sun.separation(moon)
                i = np.arctan2(sun.distance*np.sin(elongation),
                               moon.distance - sun.distance*np.cos(elongation))
                k = (1 + np.cos(i))/2.0
                mill.append(k.value)
        return mill
    
    def at_night(self, twilight1='astronomical', twilight2='astronomical',
                 time_start=0, time_end=24*3600, time_resolution=600, **kwargs):
        
        if time_start > time_end:
            self.logger.info ('Error: time_start > time_end')
            return

        kwargs = self.setkeys(kwargs)
        if kwargs is None: return
        
        _timediff = np.arange(time_start, time_end, time_resolution)
        _times = self.parse_time(**kwargs)        
        if _times is None: return
        
        elif isinstance(_times, astropy.time.Time):
            kwargs['time'] = _times + astropy.time.TimeDelta(_timediff, format=kwargs['timeu'])            
            _sunalt = self.calc_sun_height(**kwargs)
            _res = self.at_night_single(_sunalt, twilight1, twilight2)
            if _res is None:
                return
            else:
                return kwargs['time'][_res[0]], kwargs['time'][_res[1]]
        else:
            _resl = []
            for _ii in range(len(_times)):
                kwargs['time'] = _times[_ii] + astropy.time.TimeDelta(_timediff, format=kwargs['timeu'])
                _sunalt = self.calc_sun_height(**kwargs)
                _res = self.at_night_single(_sunalt, twilight1, twilight2)
                if _res is None:
                    continue
                else:
                    _nl = kwargs['time'][_res[0]], kwargs['time'][_res[1]]
                    _resl.append(_nl)
            return _resl
        
    def at_night_single(self, _sunalt, twilight1, twilight2):
        
        def get_twilight(_type):        
            if _type == 'civil': return -6*u.deg
            elif _type == 'nautical': return -12*u.deg
            elif _type == 'astronomical': return -18*u.deg
            else: return

        if get_twilight(twilight1) is None or \
           get_twilight(twilight1) is None:
            self.logger.info ('Error: twilight setting wrong')
            return
        
        _start, _end = None, None
        for _ii in range(len(_sunalt)):
            _isunalt = _sunalt[_ii]
            if not _start is None and not _end is None:
                continue
            if _isunalt <= get_twilight(twilight1) and \
               _isunalt <= get_twilight(twilight2) and _start is None:
                _start = _ii
                continue
            if _isunalt >= get_twilight(twilight2) and \
               not _start is None and _end is None:
                _end = _ii
                
        if _start is None or _end is None:
            self.logger.info ('Error: night setting wrong')
            return
        else:                
            return _start, _end
            
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

        observatories = dict(lco=lco, subaru=subaru, aao=aao, vbo=vbo, apo=apo,
                             keck=keck, kpno=kpno, lapalma=lapalma)

        if string.lower() in observatories:
            return observatories[string.lower()]
        else:
            return
