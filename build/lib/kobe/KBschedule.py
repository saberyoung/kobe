#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/schedule.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2020
# Last Modified Date: 01.03.2020
# Last Modified By  : syang <saberyoung@gmail.com>

import time
import logging
import astropy.coordinates
import astropy.units as u
import astropy.time
import numpy as np
from kobe import trigger, observatory, candidates

__all__ = ['frame', 'schedule']


class frame(observatory, candidates):
    """define an observaroty.        
    """    

    frame=None
    '''
    define observatory frame
    '''
    
    obstime = None
    '''
    define obstime for telescope
    '''
    
    def parse_frame(self, time=None, timeu='sec',
                    pressure=0.0, temperature=100.0,
                    relative_humidity= 0.0, obswl=1.0):
        """calculate altazimuth

        Parameters
        ----------   
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
        pressure           : `int`, `float`
        pressure in observatory, unit in hpa
        temperature        : `int`, `float`
        temperature in observatory, unit in C deg
        relative_humidity  : `int`, `float`
        humidity in observatory
        obswl : `~astropy.units.Quantity`
        Wavelength of the observation used in the calculation. 
        """                 
        assert not self.location is None, 'set telescope location first'          
        self.parse_time(time=time, timeu=timeu)
        if self.obstime is None: return
        self.frame = astropy.coordinates.AltAz(obstime=self.obstime,
                location=self.location, pressure=pressure*u.hPa,
                temperature=temperature*u.deg_C, obswl=obswl*u.micron,
                relative_humidity=relative_humidity)                

    def parse_time(self, time=None, timeu='sec'):
        '''parse time for astronomical calculation

        Parameters
        ----------    
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
        '''        
        if time is None:
            self.logger.info('obstime is not set')
            return        
        if self.is_seq(time) and not type(time) is str:
            _times = []
            for _time in time:
                _t = self.parse_time_single(_time, timeu)
                if not _t is None: _times.append(_t)
            self.obstime = np.array(_times)
            
        else:
            self.obstime = self.parse_time_single(time, timeu)

    def parse_time_single(self, _time, _timeu):
        '''support to parse_time
        '''        
        if _time == 'now':
            return astropy.time.Time.now()
        
        elif isinstance(_time, astropy.time.Time):
            return _time

        try:            
            return astropy.time.Time.now() + \
                astropy.time.TimeDelta(float(_time), format=_timeu)
        except:
            pass
        
        try:
            return astropy.time.Time(_time, scale='utc')
        except:
            self.logger.info ('Warning: time %s format wrong'%_time)
            
    def calc_sun_height(self, time=None, timeu='sec', pressure=0.0,
                        temperature=100.0, relative_humidity= 0.0, obswl=1.0,
                        show=False, **kwargs):        
        """calculate sun height

        Parameters
        ----------   
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
        pressure           : `int`, `float`
        pressure in observatory, unit in hpa
        temperature        : `int`, `float`
        temperature in observatory, unit in C deg
        relative_humidity  : `int`, `float`
        humidity in observatory
        obswl : `~astropy.units.Quantity`
        Wavelength of the observation used in the calculation. 

        returns
        ----------
        table :         `astropy`
                sun height        
        """ 
        self.parse_frame(time=time, timeu=timeu,
                pressure=pressure, temperature=temperature,
                relative_humidity=relative_humidity, obswl=obswl)
        if self.frame is None:return
        
        if isinstance(self.obstime, astropy.time.Time):            
            sun_loc = astropy.coordinates.get_sun(self.obstime)            
            altaz = astropy.coordinates.AltAz(obstime=self.obstime, location=sun_loc)     
            sun_altaz = sun_loc.transform_to(self.frame)
            sunalt = sun_altaz.alt
            
        else:            
            sunalt = []            
            for _ii in range(len(self.obstime)):                                 
                sun_loc = astropy.coordinates.get_sun(self.obstime[_ii])
                altaz = astropy.coordinates.AltAz(obstime=self.obstime[_ii], location=sun_loc)     
                sun_altaz = sun_loc.transform_to(self.frame[_ii])                
                sunalt.append(sun_altaz.alt)
                
        if show: self.obsplot(self.obstime, sunalt, **kwargs)
        return sunalt

    def calc_moon_illumination(self, ephemeris=None, time=None,
                               timeu='sec', show=False, **kwargs):        
        """
        Calculate fraction of the moon illuminated.

        Parameters
        ----------        
        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            `~astropy.coordinates.solar_system_ephemeris` (which is
            set to 'builtin' by default).
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

        Returns
        -------
        k : float
            Fraction of moon illuminated
        """        
        self.parse_time(time=time, timeu=timeu)
        if self.obstime is None: return
        
        # Calculate lunar orbital phase in radians. 
        if isinstance(self.obstime, astropy.time.Time):            
            sun = astropy.coordinates.get_sun(self.obstime)
            moon = astropy.coordinates.get_moon(self.obstime,
                                        ephemeris=ephemeris)
            elongation = sun.separation(moon)
            i = np.arctan2(sun.distance*np.sin(elongation),
                           moon.distance - sun.distance*np.cos(elongation))
            k = (1 + np.cos(i))/2.0
            mill = k.value
    
        else:
            mill = []            
            for _ii in range(len(self.obstime)):                
                sun = astropy.coordinates.get_sun(self.obstime[_ii])
                moon = astropy.coordinates.get_moon(self.obstime[_ii],
                                                    ephemeris=ephemeris)
                elongation = sun.separation(moon)
                i = np.arctan2(sun.distance*np.sin(elongation),
                               moon.distance - sun.distance*np.cos(elongation))
                k = (1 + np.cos(i))/2.0
                mill.append(k.value)

        if show: self.obsplot(self.obstime, mill, **kwargs)
        return mill
    
    def at_night(self, twilight1='astronomical', twilight2='astronomical',
                 ephemeris=None, time=None, time_start=0, time_end=24*3600,
                 time_resolution=600, timeu='sec', pressure=0.0,
                 temperature=100.0, relative_humidity= 0.0, obswl=1.0):        
        """
        Calculate fraction of the moon illuminated.

        Parameters
        ----------        
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
        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            `~astropy.coordinates.solar_system_ephemeris` (which is
            set to 'builtin' by default).
        pressure           : `int`, `float`
        pressure in observatory, unit in hpa
        temperature        : `int`, `float`
        temperature in observatory, unit in C deg
        relative_humidity  : `int`, `float`
        humidity in observatory
        obswl : `~astropy.units.Quantity`
        Wavelength of the observation used in the calculation. 
        """
        
        assert time_start <= time_end, 'time_start > time_end'          
        _timediff = np.arange(time_start, time_end, time_resolution)
        self.parse_time(time=time, timeu=timeu)       
        if self.obstime is None:return
        
        if isinstance(self.obstime, astropy.time.Time):
            _time = self.obstime + astropy.time.TimeDelta(_timediff, format=timeu)            
            _sunalt = self.calc_sun_height(time=_time, timeu=timeu,
                        pressure=pressure, temperature=temperature,
                        relative_humidity=relative_humidity, obswl=obswl)            
            _res = self.at_night_single(_sunalt, twilight1, twilight2)
            if _res is None:
                return
            else:
                return self.obstime[_res[0]], self.obstime[_res[1]]
        else:
            _resl = []
            for _ii in range(len(self.obstime)):
                _time = self.obstime[_ii] + astropy.time.TimeDelta(_timediff, format=timeu)
                _sunalt = self.calc_sun_height(time=_time, timeu=timeu,
                        pressure=pressure, temperature=temperature,
                        relative_humidity=relative_humidity, obswl=obswl)
                _res = self.at_night_single(_sunalt, twilight1, twilight2)
                if _res is None:
                    continue
                else:
                    _nl = self.obstime[_res[0]], self.obstime[_res[1]]
                    _resl.append(_nl)
            return _resl
    
    def at_night_single(self, _sunalt, twilight1, twilight2):
        
        def get_twilight(_type):        
            if _type == 'civil': return -6*u.deg
            elif _type == 'nautical': return -12*u.deg
            elif _type == 'astronomical': return -18*u.deg
            else: return

        assert not get_twilight(twilight1) is None, 'twilight1 setting wrong'
        assert not get_twilight(twilight1) is None, 'twilight2 setting wrong'
        
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
        
    def altaz(self, ra=None, dec=None, theta=None, phi=None,
              pressure=0.0, temperature=100.0, relative_humidity= 0.0,
              obswl=1.0, time=None, timeu='sec'):
        """calculate altazimuth        
        """
        self.parse_frame(time=time, timeu=timeu, obswl=obswl,
                         pressure=pressure, temperature=temperature,
                         relative_humidity=relative_humidity)
        if self.frame is None:return
        
        radecs = self.radecs(ra=ra, dec=dec, theta=theta, phi=phi)        
        if radecs is None: return        
        radecs = astropy.coordinates.SkyCoord(ra=radecs[0]*u.deg, dec=radecs[1]*u.deg)        
        
        if self.is_seq(self.frame.obstime.value) and \
           not type(self.frame.obstime.value) is str:
            altaz = []            
            for _ii in range(len(self.frame.obstime.value)):
                altaz.append(radecs.transform_to(self.frame[_ii]))
                
        else:                      
            altaz = radecs.transform_to(self.frame)            
        return altaz

    def altaz_pointings(self, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        assert not self.pointings is None
        assert not self.pointings.data is None
        return self.altaz(ra=self.pointings.data['ra'],
                dec=self.pointings.data['dec'], pressure=pressure,
                temperature=temperature, relative_humidity=relative_humidity,
                time=time, timeu=timeu)
    
    def altaz_candidates(self, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        assert 'candidates' in self.__dict__
        assert not self.candidates is None
        return self.altaz(ra=self.candidates['ra'], dec=self.candidates['dec'],
                          pressure=pressure, temperature=temperature,
                          relative_humidity=relative_humidity, time=time, timeu=timeu)
    
    def calc_airmass(self, ra=None, dec=None, theta=None, phi=None,
            pressure=0.0, temperature=100.0, relative_humidity= 0.0,
            obswl=1.0, time=None, timeu='sec', show=False, **kwargs):
        
        """calculate airmass for tilings        
        """
        altaz = self.altaz(ra=ra, dec=dec, theta=theta, phi=phi,
                pressure=pressure, temperature=temperature,
                relative_humidity=relative_humidity, time=time,
                timeu=timeu)
        
        if altaz is None:
            return

        elif type(altaz) is list:
            airmass = []
            for _altaz in altaz: airmass.append(_altaz.secz)            

        else:
            airmass = altaz.secz

        if show: self.cooplot(self.obstime, airmass,
                        ids=['%.2f %.2f'%(i,j) for i,j in zip(ra,dec)], **kwargs)
        return airmass

    def calc_airmass_pointings(self, pressure=0.0, temperature=100.0,
                relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec',
                show = False, **kwargs):
        assert not self.pointings is None
        assert not self.pointings.data is None
        return self.calc_airmass(ra=self.pointings.data['ra'],
                dec=self.pointings.data['dec'], pressure=pressure,
                temperature=temperature, relative_humidity=relative_humidity,
                time=time, timeu=timeu, show=show, **kwargs)
    
    def calc_airmass_candidates(self, pressure=0.0, temperature=100.0,
                relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec',
                show=False, **kwargs):
        assert 'candidates' in self.__dict__
        assert not self.candidates is None
        return self.calc_airmass(ra=self.candidates['ra'], dec=self.candidates['dec'],
                        pressure=pressure, temperature=temperature,
                        relative_humidity=relative_humidity, time=time,
                        timeu=timeu, show=show, **kwargs)
    
    def calc_solar(self, solarobj=None, ra=None, dec=None,
            theta=None, phi=None, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        """calculate distance of tilings to solar object        

        Parameters
        ----------    
        solarobj           : `str`
        name of solar object
        options: earth, sun, earth-moon-barycenter, moon,
                 mercury,venus,mars,jupiter, saturn,
                 uranus,neptune       

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
        _list = ['earth', 'sun', 'earth-moon-barycenter', 'moon',
                 'mercury','venus','mars','jupiter', 'saturn',
                 'uranus','neptune']
        if not solarobj in _list:
            self.logger.info ('Error: solar obj select from %s' % _list)
            return
        radecs = self.radecs(ra=ra, dec=dec, theta=theta, phi=phi)
        if radecs is None: return        
        radecs = astropy.coordinates.SkyCoord(ra=radecs[0]*u.deg, dec=radecs[1]*u.deg)
        
        altaz = self.altaz(ra=ra, dec=dec, theta=theta, phi=phi,
                pressure=pressure, temperature=temperature,
                relative_humidity=relative_humidity, time=time,
                timeu=timeu)
        if altaz is None: return
        
        if not isinstance(self.obstime, astropy.time.Time):
            _sdist = []
            for _ii in range(len(self.obstime)):
                _sloc = astropy.coordinates.get_body(solarobj, self.obstime[_ii],
                                    self.location).transform_to(altaz[_ii])
                _sdist.append(radecs.separation(_sloc))                
        else:            
            _sloc = astropy.coordinates.get_body(solarobj, self.obstime,
                                        self.location).transform_to(altaz)
            _sdist = radecs.separation(_sloc)
        return _sdist

    def calc_solar_pointings(self, solarobj=None, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        assert not self.pointings is None
        assert not self.pointings.data is None
        return self.calc_solar(solarobj=solarobj, ra=self.pointings.data['ra'],
                dec=self.pointings.data['dec'], pressure=pressure,
                temperature=temperature, relative_humidity=relative_humidity,
                time=time, timeu=timeu)
    
    def calc_solar_candidates(self, solarobj=None, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        assert 'candidates' in self.__dict__
        assert not self.candidates is None
        return self.calc_solar(solarobj=solarobj,
                    ra=self.candidates['ra'], dec=self.candidates['dec'],
                    pressure=pressure, temperature=temperature,
                    relative_humidity=relative_humidity, time=time, timeu=timeu)
    
    def calc_galactic(self, ra=None, dec=None,
            theta=None, phi=None, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        
        altaz = self.altaz(ra=ra, dec=dec, theta=theta, phi=phi,
                pressure=pressure, temperature=temperature,
                relative_humidity=relative_humidity, time=time,
                timeu=timeu)

        if altaz is None:
            return

        elif type(altaz) is list:
            gal_altaz = []
            for _altaz in altaz:
                gal_altaz.append(_altaz.transform_to(astropy.coordinates.Galactic).b)
                
        else:        
            gal_altaz = abs(altaz.transform_to(astropy.coordinates.Galactic).b)            
        return gal_altaz

    def calc_galactic_pointings(self, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        assert not self.pointings is None
        assert not self.pointings.data is None
        return self.calc_galactic(ra=self.pointings.data['ra'],
                dec=self.pointings.data['dec'], pressure=pressure,
                temperature=temperature, relative_humidity=relative_humidity,
                time=time, timeu=timeu)
    
    def calc_galactic_candidates(self, pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0, time=None, timeu='sec'):
        assert 'candidates' in self.__dict__
        assert not self.candidates is None
        return self.calc_galactic(ra=self.candidates['ra'], dec=self.candidates['dec'],
                            pressure=pressure, temperature=temperature,
                            relative_humidity=relative_humidity, time=time, timeu=timeu)    

    def rankp_visibility(self, mode=1, nside=64, nest=False,
            threshold=1, sort=1, limairmass=2.5, limgalactic=10,
            limsolar={'sun':30, 'moon':None},  exptime=60, ovhtime=100,
            filters='r', twilight1='astronomical', twilight2='astronomical',
            ephemeris=None, time=None, time_start=0, time_end=24*3600,
            time_resolution=600, timeu='sec', pressure=0.0, temperature=100.0,
            relative_humidity= 0.0, obswl=1.0):
        """rank pointings (conservative approach): rank highest airmass

        Parameters
        ----------           
        mode     : `int`
                  1. strict approach: follow above strategy strictly
                  2. adjacent approach: start from the westest point 
                     or highest probability point and arange the next 
                     pointing either adjacent or closest to the previous one
        """                                  
        # define start point        
        _idstart = np.argmin(_data['ra'])
        
        if _res is None: return
        assert len(_res) == 2
        _tstart, _tend = _res
        assert isinstance(_tstart, astropy.time.Time)
        assert isinstance(_tend, astropy.time.Time)

        # define telescope list
        _tellist = {}
        if 'observatories' in self.__dict:
            
            self.location = location
            
             # define time period
            _res = self.at_night(self, twilight1=twilight1, twilight2=twilight2,
                        ephemeris=ephemeris, time=time, time_start=time_start,
                        time_end=time_end, time_resolution=time_resolution,
                        timeu=timeu, pressure=pressure, temperature=temperature,
                        relative_humidity= relative_humidity, obswl=obswl)
        
        elif 'telescopes' in self.__dict:
            assert 'location' in self.__dict__
            
        else:
            assert 'location' in self.__dict__
            if 'pointings' in self.__dict__:            
                assert not self.pointings.data is None, 'Warning: tiling data not parsed'
                _data = self.pointings.data
            else:
                assert 'location' in self.__dict__
                assert not self.pointings.data is None, 'Warning: tiling data not parsed'
                _data = self.data
                
        for _nt,_tt in enumerate(_glist):

            _tellist[_nt] = {}

            # define observatory and time per frame
            _obs, _tf, [_exp,_rot,_numfield,_numexp] = \
                                self.def_obs(_tt)

            # define time
            _timenow, _jd = self.def_time(_tt)

            # time last
            _deltat = float(_tt['observe']['timelast'])            
            _deltat = astropy.time.TimeDelta(_deltat*3600., format='sec') 

            # time list
            _timel = np.arange(_timenow, _timenow+_deltat, _tf)

            # store them
            # innight remove time, when the sun is high
            _tellist[_nt]['timelist'] = \
                    pst.innight(_timel,_obs,\
                    float(_tt['observe']['limsun']))
            _tellist[_nt]['tf'] = _tf
            _tellist[_nt]['name'] = _tt['telescope']['name']
            _tellist[_nt]['ra'] = []
            _tellist[_nt]['dec'] = []
            if _glist[_nt]['pointings']['scheduler'] == 'T':
                _tellist[_nt]['fovw'] = []
                _tellist[_nt]['fovh'] = []
            _tellist[_nt]['obs'] = _obs

        # define time list
        _timelist = []
        for zz in _tellist:
            _timelist.extend(_tellist[zz]['timelist'])

        _timelist = np.unique(_timelist)
        _timelist = np.sort(_timelist)

        # For slew angle case, start from latest fields
        ralatest, declatest = {}, {}
        for ntel,tel in enumerate(_glist):
            ralatest[ntel], declatest[ntel] = None, None

        # field assignment start
        for _nt,_time in enumerate(_timelist):

            # stop till: 
            # - 1. the end of time;
            # - 2. used up all available pointings

            _timedif = _time - _timelist[0]
            if self.verbose:
                print ('\n\t+ %.2f hs'%(_timedif.value*24))


                
        if mode == 1:
            ipx = np.argsort(_data['ra'])
            
        elif mode == 2:            
            idx = np.arange(len(_data['ra']))                                
            ipx, _idr = np.array([], dtype=int), _idstart
        
            while len(ipx) < len(idx)-1:
                ipx = np.append(ipx, _idr)  
                _ra, _dec = _data['ra'][_idr],_data['dec'][_idr]            
                _idl = np.array([], dtype=int)
                for ii in idx:                
                    if not ii in ipx: _idl = np.append(_idl, ii)            
                _ral, _decl = _data['ra'][_idl], _data['dec'][_idl]            
                _id = self.adjacent_move(_ra,_dec,_ral,_decl,
                                         threshold=threshold,sort=sort)
                _idr = idx[np.logical_and(_data['ra'] == _ral[_id],
                                          _data['dec'] == _decl[_id])][0]
                np.append(ipx, _idr) 

        else:
            self.logger.info ('wrong mode')
            return
            
        self.pointings.data = _data[ipx]
        self.pointings.data['n'] = np.arange(len(self.pointings.data['n']), dtype=int)

        
class schedule(trigger, frame):

    obs = None
    '''
    define telescope observing blocks
    '''
    
    observatories = {}
    '''
    define a list of observatories
    '''

    defkwargs = {**trigger.defkwargs, **frame.defkwargs}
    '''
    define default keys
    '''
    
    def __init__(self, logger=None, **kwargs):        
        super(schedule, self).__init__(logger=logger, **kwargs)

    def process_gcn(self, payload, root):        
        # parse trigger
        self.root(root)

    def add_observatory(self, *args, **kwargs):
                
        #self.set_obsname(*args, **kwargs)
        self.parse_location(*args, **kwargs)
            
        kwargs = self.setkeys(kwargs)
        self.check_obs()
        self.telescopes = {}
        if self.obsname in self.observatories and not kwargs['clobber']:
            self.logger.info ('observatory %s already added'%self.obsname)
        else:
            self.observatories[self.obsname] = [self.telescopes, self.location]
            self.logger.info ('observatory %s added'%self.obsname)
            
    def del_observatory(self, obsname):
                        
        del self.observatories[obsname]

    def to_observatory(self, obsname):        
        assert obsname in self.observatories
        self.obsname = obsname        
        self.telescopes, self.location = self.observatories[obsname]

    def update_observatory(self, obsname):        
        assert obsname in self.observatories
        self.observatories[obsname] = [self.telescopes, self.location]

    def genobs(self, filters='r', ovhtime=120, exptime=45, mode=1):
        """set telescope name

        Parameters
        ----------
        filter:
                telescope filter
                if multiple, split with comma
        exptime:        
                estimated time for telescope exposure per frame + field slewing
                unit in second               
        ovhtime:       
                estimated time for telescope overhat
                unit in second
        """
        try: self.pointings.data
        except: return
        assert self.obs is None
                
        ff=[]
        for _nf, _filt in enumerate(filters.split(',')):
            if _nf == 0:                                
                _data = self.pointings.data                                
            else:                                
                _data = vstack([_data, self.pointings.data], join_type='outer')
            ff += [_filt]*len(self.pointings.data)                        
        ffc = Column(ff, name='filter')
                
        self.obs = _data
        self.obs.add_column(ffc)

        # arrange times
        tt = Column(np.zeros(len(self.obs['n'])), name='time')
        self.obs.add_column(tt)

        _time = 0
        for _nn in np.unique(self.obs['n']):
            _idx = np.where(self.obs['n']==_nn)
            if mode==1:                                
                ra, dec = self.obs[_idx]['ra'][0], \
                    self.obs[_idx]['dec'][0]
                ral, decl = self.obs[_idx]['ra'][1:], \
                    self.obs[_idx]['dec'][1:]
                idl = self.spherical_angle(ra,dec,ral,decl,sort=True)
                try:_tmptab.add_row(self.obs[_idx][0])
                except:_tmptab = Table(self.obs[_idx][0])
                
                _time += ovhtime
                _tmptab[-1]['time'] += _time
                                
                for _tt in self.obs[_idx][1:][idl]:
                    _time += exptime
                    _tt['time'] += _time
                    _tmptab.add_row(_tt)                                
            else:
                for _ff in np.unique(self.obs[_idx]['filter']):
                    _idy =np.where(self.obs[_idx]['filter']==_ff)
                    ra, dec = self.obs[_idx][_idy]['ra'][0], \
                        self.obs[_idx][_idy]['dec'][0]
                    ral, decl = self.obs[_idx][_idy]['ra'][1:], \
                        self.obs[_idx][_idy]['dec'][1:]
                    idl = self.spherical_angle(ra,dec,ral,decl,sort=True)
                    try:_tmptab.add_row(self.obs[_idx][_idy][0])
                    except:_tmptab = Table(self.obs[_idx][_idy][0])
                    _time += ovhtime
                    _tmptab[-1]['time'] += _time
                    for _tt in self.obs[_idx][_idy][1:][idl]:
                        _time += exptime
                        _tt['time'] += _time
                    _tmptab.add_row(_tt)   
        self.obs = _tmptab
                
    def prioritization(self):
        
        start_time = time.time()        
        assert len(self.observatories) > 0
        
        _tellist = {}
        for _nt,_tt in enumerate(_glist):

            _tellist[_nt] = {}

            # define observatory and time per frame
            _obs, _tf, [_exp,_rot,_numfield,_numexp] = \
                                self.def_obs(_tt)

            # define time
            _timenow, _jd = self.def_time(_tt)

            # time last
            _deltat = float(_tt['observe']['timelast'])            
            _deltat = astropy.time.TimeDelta(_deltat*3600., format='sec') 

            # time list
            _timel = np.arange(_timenow, _timenow+_deltat, _tf)

            # store them
            # innight remove time, when the sun is high
            _tellist[_nt]['timelist'] = \
                    pst.innight(_timel,_obs,\
                    float(_tt['observe']['limsun']))
            _tellist[_nt]['tf'] = _tf
            _tellist[_nt]['name'] = _tt['telescope']['name']
            _tellist[_nt]['ra'] = []
            _tellist[_nt]['dec'] = []
            if _glist[_nt]['pointings']['scheduler'] == 'T':
                _tellist[_nt]['fovw'] = []
                _tellist[_nt]['fovh'] = []
            _tellist[_nt]['obs'] = _obs

        # define time list
        _timelist = []
        for zz in _tellist:
            _timelist.extend(_tellist[zz]['timelist'])

        _timelist = np.unique(_timelist)
        _timelist = np.sort(_timelist)

        # For slew angle case, start from latest fields
        ralatest, declatest = {}, {}
        for ntel,tel in enumerate(_glist):
            ralatest[ntel], declatest[ntel] = None, None

        # field assignment start
        for _nt,_time in enumerate(_timelist):

            # stop till: 
            # - 1. the end of time;
            # - 2. used up all available pointings

            _timedif = _time - _timelist[0]
            if self.verbose:
                print ('\n\t+ %.2f hs'%(_timedif.value*24))

            for ntel,tel in enumerate(_glist):

                # check telescope available or not
                if _time in _tellist[ntel]['timelist']:
                    if self.verbose:
                        print ('\n\t - tel: %s avalable now'%\
                               tel['telescope']['name'])
                else: # not available
                    if self.verbose:
                        print ('\n\t - tel: %s not avalable now'%\
                               tel['telescope']['name'])
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # if reached pointing limit
                try:
                    nf = float(_glist[ntel]['observe']['limfields'])
                    if sum([bb is not None for bb in \
                            _tellist[ntel]['ra']]) >= nf:
                        if self.verbose:
                            print ('\t - tel: %s reached pointing limit'%\
                                   tel['telescope']['name'])
                        index = np.where(_tellist[ntel]['timelist']==_time)
                        _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                        continue
                except: 
                    pass

                # if used up all pointings
                if sum([bb is not None for bb in \
                    _tellist[ntel]['ra']]) >= \
                    len(_glist[ntel]['ra']):
                    if self.verbose:
                        print ('\t - tel: %s finished all pointings'%\
                               tel['telescope']['name'])
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # remove duplicated fields
                _tralist,_tdeclist = [],[]
                for _ra,_dec in zip(tel['ra'],tel['dec']):
                    
                    _use = True
                    if tel['pointings']['scheduler'] == 'T':
                        # remain pointings for followup
                        # index different for different telescopes
                        # check if 2 pointings have overlap

                        # current pointing
                        fovw, fovh = float(tel['telescope']['fovw']),\
                                     float(tel['telescope']['fovh'])

                        # remove input ra,dec which should be high priority
                        _idx = pst.overlapregioncut(_ra,_dec,_cralist,_cdeclist)
                        _cov = pst.overlapregion(_ra,_dec,fovw,fovh,\
                                                 np.array(_cralist)[_idx],\
                                                 np.array(_cdeclist)[_idx],\
                                                 np.array(_cfovw)[_idx],\
                                                 np.array(_cfovh)[_idx])
                        if not _cov is None:
                            _frac = _cov/fovw/fovh
                            if _frac > float(tel['pointings']['uarea']):
                                _use = False

                        # remove already done for the same priority
                        # different telescopes
                        for _ntel in _tellist:
                            _idx = pst.overlapregioncut(_ra,_dec,\
                                    _tellist[_ntel]['ra'],_tellist[_ntel]['dec'])
                            _cov = pst.overlapregion(_ra,_dec,fovw,fovh,\
                                    np.array(_tellist[_ntel]['ra'])[_idx],\
                                    np.array(_tellist[_ntel]['dec'])[_idx],\
                                    np.array(_tellist[_ntel]['fovw'])[_idx],\
                                    np.array(_tellist[_ntel]['fovh'])[_idx])
                            if not _cov is None:
                                _frac = _cov/fovw/fovh
                                if _frac > float(tel['pointings']['uarea']):
                                    _use = False

                    elif _glist[ntel]['pointings']['scheduler'] == 'G':
                        # remain galaxies for followup
                        # index same for different telescopes
                        # check if 2 galaxies are same

                        # remove input ra,dec which should be high priority
                        for _ra0,_dec0 in zip(_cralist,_cdeclist):
                            if _ra == _ra0 and _dec == _dec0:
                                _use = False

                        # remove already done for the same priority
                        # different telescopes
                        for _ntel in _tellist:
                            for _ra0,_dec0 in zip(_tellist[_ntel]['ra'],\
                                                  _tellist[_ntel]['dec']):
                                if _ra == _ra0 and _dec == _dec0:
                                    _use = False

                    if _use:
                        _tralist.append(_ra)
                        _tdeclist.append(_dec)

                if self.verbose:
                    print ('\t - %s fields remained for tel %s, after removing duplication'%\
                           (len(_tralist),tel['telescope']['name']))
                if len(_tralist)==0:
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # ra dec
                radecs = astropy.coordinates.SkyCoord(ra=_tralist*u.deg, \
                                                      dec=_tdeclist*u.deg)

                # read observatory
                observatory = _tellist[ntel]['obs']

                # Alt/az reference frame at observatory, now
                frame = astropy.coordinates.AltAz(obstime=_time, location=observatory)

                # Transform grid to alt/az coordinates at observatory, now
                altaz = radecs.transform_to(frame)

                # Where is the sun, now?
                sun_altaz = astropy.coordinates.get_sun(_time).transform_to(altaz)

                # sun, airmass constrains                      
                _limalt,_limsun = float(tel['observe']['limalt']), \
                                  float(tel['observe']['limsun'])               

                _cond = np.logical_and(np.logical_and(altaz.secz <= _limalt, \
                                altaz.secz >= 1),sun_altaz.alt <= _limsun*u.deg)
                gradecs = radecs[_cond]

                # go on or not
                if self.verbose:
                    print ('\t - %s fields visible for tel %s, after sun+airmass constrain'%\
                           (len(gradecs),tel['telescope']['name']))
                if len(gradecs) == 0: 
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # moon constrains
                # define seperation to the moon
                _limmoon = tel['observe']['limmoon']
                if _limmoon == 'auto':
                    if 'T' in str(_time):
                        _year, _month, _day = str(_time.value).split('T')[0].split('-')
                    else:
                        _year, _month, _day = str(_time.value).split()[0].split('-')
                    _date, _status, _light = moon_phase(int(_month),int(_day),int(_year))

                    ''' limitation for the moon / TBD '''
                    if _light>=80: _limmoon = 30
                    elif _light<80 and _light>=40: _limmoon = 20
                    elif _light<40 and _light>=10: _limmoon = 10
                    else:_limmoon = 1
                else:
                    try: _limmoon = float(_limmoon)
                    except: 
                        print ('!!! Warning: moon limitation wrong... set default')
                        _limmoon = 30

                _limmoon *= 3600 # to deg
                _loc = astropy.coordinates.get_body('moon', \
                        _time, observatory).transform_to(altaz)
                sep = gradecs.separation(_loc)              
                _length = len(gradecs[np.where(sep.arcsecond<_limmoon)])
                if _length>0:
                    if self.verbose:
                        print('\t - remove %s sources due to moon'%(_length))                   
                    gradecs = gradecs[sep.arcsecond>_limmoon]

                # go on or not
                if self.verbose:
                    print ('\t - %s fields visible for tel %s, after moon constrain'%\
                           (len(gradecs),tel['telescope']['name']))
                if len(gradecs) == 0: 
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # Where is the other sources inside solor system
                _solorobj = tel['observe']['limsobj'].split(',')  
                try:_limsolor = float(tel['observe']['limsobjr'])*3600 # to deg
                except:_limsolor=False
                if _limsolor:
                    for _source in _solorobj:
                        # astropy.coordinates.solar_system_ephemeris.bodies
                        if _source in ['mercury','venus','mars',\
                                       'jupiter','saturn','uranus','neptune']:
                            _loc = astropy.coordinates.get_body(_source, \
                                    _time, observatory).transform_to(altaz)
                            sep = gradecs.separation(_loc)
                            _length = len(gradecs[np.where(sep.arcsecond<_limsolor)])
                            if _length>0:
                                if self.verbose:
                                    print('\t - remove %s sources due to %s'%(_length,_source))
                                gradecs = gradecs[sep.arcsecond>_limsolor]
                        else:
                            if self.verbose:
                                print ('\t!!! Warning: unknown source %s'%_source)

                # go on or not
                if self.verbose:
                    print ('\t - %s fields visible for tel %s, after solor obj constrain'%\
                           (len(gradecs),tel['telescope']['name']))
                if len(gradecs) == 0: 
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # ranking fields
                _order = int(tel['observe']['order'])

                # rank with prob
                if _order == 1: pass
                elif _order == 2: # from west to east
                    _idrank = np.argsort(gradecs.ra)
                    gradecs = gradecs[_idrank]

                elif _order in [3,4]: # consider slewing angle
                    if ralatest[ntel] is None and declatest[ntel] is None:
                        # initialize start field
                        if _order==3: # start from west
                            _ids = np.argmin(gradecs.ra)
                        elif _order==4:
                            # start from high ranking
                            # Note: the input radecs are already ranked with ranking
                            _ids = 0
                        ras,decs = gradecs.ra[_ids].deg,\
                                   gradecs.dec[_ids].deg
                    else:
                        ras,decs = ralatest[ntel], declatest[ntel]
                    index = pst.slew_angle(ras,decs,gradecs.ra.deg,\
                                           gradecs.dec.deg)
                    gradecs = gradecs[index]
               
                if self.verbose:
                    print ('\t - %s fields ranked for tel %s, with order=%i'%\
                           (len(gradecs),tel['telescope']['name'],_order))

                # append fields
                for _nn,(_ra,_dec) in enumerate(zip(gradecs.ra.deg, gradecs.dec.deg)):
                    if _nn<_numfield: # append numfield fields
                        _tellist[ntel]['ra'].append(_ra)
                        _tellist[ntel]['dec'].append(_dec)
                        ralatest[ntel], declatest[ntel] = _ra, _dec
                if self.verbose:
                    print ('\t%i fields for %s\t'%\
                           (len(_tellist[ntel]['ra']),\
                            tel['telescope']['name']))                
                if _glist[ntel]['pointings']['scheduler'] == 'T':
                    for ii in range(_numfield):
                        _tellist[ntel]['fovw'].append(float(_glist[ntel]['telescope']['fovw']))
                        _tellist[ntel]['fovh'].append(float(_glist[ntel]['telescope']['fovh']))
        if self.verbose:print ('\t### time finished')
        return _tellist
