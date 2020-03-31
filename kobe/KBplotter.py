#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : kobe/KBPlotter.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

# Standard library
import logging
import os

# third-party module
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import astropy.coordinates
import astropy.time
import astropy.units as u

# self import
from kobe import utils

__all__ = ['visualization', 'vtilings', 'vgalaxies',
           'vtrigger', 'vcandidates']

class visualization(utils):
    """provide various visualization functions.
    
    Parameters
    ----------    
    rot_theta  : int, float
       *healpix.visualization* rotation angle.
       theta = pi/2 - decl (unit in deg), default: 0.
       setted to rotate localization map along longtitude direction
    rot_phi    : int, float
       *healpix.visualization* rotation angle.
       phi = ra (unit in deg), default: 0.
       rotate localization map along latitude direction
    rot_psi    :  int, float
       *healpix.visualization* rotation angle.
       the point at longitude and latitude will be at the center.   
       An additional rotation of angle psi around this direction is applied.
       default: 0.
    num          :   int, str
        *matplotlib* or *healpy* figure number.
        default: 1.
    wdir      : str   
       working directory.
       default: ./
    clobber   : bool  
       if any product has already been parsed, 
       set clobber as True will remove and redo the calculation.
       default: False

    Notes
    ----------   
    *visualization* inherients from **utils**,
    and inheriented by *vtrigger*/*vgalaxies*/*vtilings*/*vcandidates*.
    """      
        
    def _init_figure(self, **kwargs):
        """initialize/switch to a figure
        
        Parameters
        ----------        
        num          :   int or  str
            matplotlib figure number
        Matplotlib.figure Parameters              
        """                
        # init plot
        kwargs = self.setkeys(kwargs)        
        return plt.figure(num=kwargs['num'], **self.getkeys(kwargs, 'matplotlib'))    
    
    def savefig(self, filepath='tmp_kobe.png', **kwargs):
        """save a figure
        
        Parameters
        ----------
        filepath:  str
        absolute path that would be used when one want to store a figure       
        num          :   int or  str
        matplotlib figure number
        wdir      : str   
        working directory  
        clobber   : bool  
        if filepath exists, remove and redo the plot or not       
        Matplotlib.savefig Parameters
        Matplotlib.figure Parameters     
        """
        kwargs = self.setkeys(kwargs)        
        filename = '%s/%s' % (kwargs['wdir'], filepath)
        if os.path.exists(filename) and not kwargs['clobber']:
            self.logger.info ('Warning: %s already exists'%filename)
            return        
        fig = self._init_figure(**kwargs)        
        fig.savefig(filename, **self.getkeys(kwargs, 'savefig'))

    def hpmapproj(self, hpmap, showtype='m', **kwargs):
        """healpy Map projections
        
        Parameters
        ---------- 
        hpmap  :
           healpix map           
        fig          :   int or  str
           matplotlib figure number
        showtype :
           function to call        
        Healpix Parameters
        Matplotlib.figure Parameters
        """
        kwargs = self.setkeys(kwargs)
        
        # Choose color map and set background to white
        cmap = cm.YlOrRd
        cmap.set_under("w")

        # Plot GW skymap
        assert showtype in ['m','g','c','o'], 'wrong option for showtype'
        if showtype == 'm': _func = hp.mollview
        elif showtype == 'g': _func = hp.gnomview
        elif showtype == 'c': _func = hp.cartview
        elif showtype == 'o': _func = hp.orthview        
        _func(map=hpmap, fig=kwargs['num'], rot=[kwargs['rot_theta'],
                kwargs['rot_phi'], kwargs['rot_psi']], cmap=cmap,
                return_projected_map=False, **self.getkeys(kwargs, 'healpixview'))
        
    def projplot(self, theta, phi, **kwargs):
        """healpy projplot

        Parameters
        ----------          
        theta, phi
        Matplotlib.figure Parameters
        Healpix projplot Parameters
        """              
        hp.projplot(theta, phi, **self.getkeys(kwargs, 'projplot'))
        
    def closefig(self, num):
        """close a figure
        
        Parameters
        ---------- 
        Matplotlib.figure Parameters
        """        
        plt.close(num)

    def _grid(self, **kwargs):
        """generate background grids for figure

        Parameters
        ----------        
        Healpix graticule Parameters
        """                       
        hp.graticule(**self.getkeys(kwargs, 'graticule'))

    def zoomin(self, zoomx=None, zoomy=None):
        """zoomin figure

        Parameters
        ----------        
        zoomx       :  list
        zoomin in x direction
        zoomy      :  list
        zoomin in y direction
        """
        for _zoom, _func in zip([zoomx, zoomy], [plt.xlim, plt.ylim]):
            if not _zoom is None:
                assert len(_zoom) == 2, 'zoomx/y should be a 2 element list'            
                _func(_zoom[0], _zoom[1])                
            
    def plot_sky(self, showgrid=True, **kwargs):
        """plot: horizon, galactic plane, sun and moon

        Parameters
        ----------        
        Healpix graticule Parameters
        """
        kwargs = self.setkeys(kwargs)
        r = hp.Rotator(deg=True, rot=[kwargs['rot_theta'],
                        kwargs['rot_phi'], kwargs['rot_psi']])
                        
        # create figure        
        fig = self._init_figure(**kwargs)
        if showgrid:
            self._grid(**kwargs)                        
            self.plot_coord(rot = [kwargs['rot_theta'],
                        kwargs['rot_phi'], kwargs['rot_psi']], **kwargs)
            
        # plot the galactic plane        
        ra = np.arange(0,360,10)
        dec = np.zeros(len(ra))

        thetal, phil = [], []
        for _ra,_dec in zip(ra,dec):
        
            # from galactic coordinates to equatorial
            _radecs = astropy.coordinates.SkyCoord(l=_ra*u.deg,
                                        b=_dec*u.deg, frame='galactic')                

            # transform to theta phi
            _ra,_dec = _radecs.icrs.ra.deg, _radecs.icrs.dec.deg
            _theta,_phi = np.pi/2.-np.radians(_dec),np.radians(_ra)            
            thetal.append(_theta)
            phil.append(_phi)
            
        kwargs['color'] = 'k'
        kwargs['marker'] = 'x'
        kwargs['label'] = 'galactic plane'
        thetal, phil = r(thetal, phil)
        self.projplot(thetal, phil, **kwargs)        
            
        if 'obstime' in self.__dict__:
            assert isinstance(self.obstime, astropy.time.Time)
            
            # plot the sun
            _sra,_sdec = astropy.coordinates.get_sun(self.obstime).ra.deg,\
                astropy.coordinates.get_sun(self.obstime).dec.deg
            _stheta,_sphi = np.pi/2.-np.radians(_sdec),np.radians(_sra)
            
            # visualization
            kwargs['color'] = 'y'
            _stheta,_sphi = r(_stheta,_sphi)
            hp.projtext(_stheta,_sphi, '$\odot$', lonlat=False,
                        **self.getkeys(kwargs, 'projtext'))
            
            # plot the moon
            for _nd in np.arange(0,5,1):
                _timeplus = _nd*24
                timelater = self.obstime + astropy.time.TimeDelta(_timeplus*3600, format='sec')
                _mra,_mdec = astropy.coordinates.get_moon(timelater).ra.deg,\
                    astropy.coordinates.get_moon(timelater).dec.deg
                _mtheta,_mphi = np.pi/2.-np.radians(_mdec),np.radians(_mra)

                # visualization
                kwargs['color'] = 'b'
                _mtheta,_mphi = r(_mtheta,_mphi)
                hp.projtext(_mtheta,_mphi,'$\oplus$', lonlat=False,
                            **self.getkeys(kwargs, 'projtext'))
                hp.projtext(_mtheta,_mphi-.1,'%id'%_nd, lonlat=False,
                            **self.getkeys(kwargs, 'projtext'))
            
            # write labels
            xx,yy = -1.5,1.
            plt.text(xx, .9, 'sun: $\odot$',fontsize=20,\
                     ha="center", va="center", color='y')
            plt.text(xx, 1., 'moon: $\oplus$',fontsize=20,\
                     ha="center", va="center", color='b')

        # plot the horizon for each observatories
        if not 'observatories' in self.__dict__: return
        for obs in self.observatories:
            
            for tel in self.observatories[obs]:
                
                name, lat, lon, alt = tel.name, tel.conf['lat'], tel.conf['lon'], tel.conf['alt']
                # define observatory
                observatory = astropy.coordinates.EarthLocation(lat=lat*u.deg,
                                                    lon=lon*u.deg, height=alt*u.m)
                _smlabel=True
                for _timedelta in np.arange(0,360,1):                               
                    AltAzcoordiantes = astropy.coordinates.SkyCoord(alt = 0*u.deg, 
                                    az = 0*u.deg + _timedelta*u.deg, obstime = obstime,
                                    frame = 'altaz', location = observatory)
                    ra,dec = AltAzcoordiantes.icrs.ra.deg, AltAzcoordiantes.icrs.dec.deg
                    theta, phi = np.pi/2.-np.radians(dec),np.radians(ra)                
                    if _smlabel:
                        hp.projplot(theta,phi,'.', color = _colorlist[nt], \
                                    coord=coord, ms = 2, label='%s horizon now'%name)
                        _smlabel=False
                    else:
                        hp.projplot(theta,phi,'.', color = _colorlist[nt], \
                                    coord=coord, ms = 2)
        
    def plot_coord(self, rot=[0,0,0], bkgcolor='k', **kwargs):
        """show specific coordinates in skymap

        Parameters
        ----------        
        rot   :  list
           rotion scheme
        Healpix graticule Parameters
        """        
        kwargs = self.getkeys(kwargs, 'projtext')
        kwargs['color'] = bkgcolor
        for _t in [60,120,180,240,300,360]:
            # deg to hms
            c= astropy.coordinates.SkyCoord(ra=_t*u.degree,
                                dec=0*u.degree, frame='icrs')          
    
            # visualization            
            hp.projtext(_t-rot[0], 0-rot[1], '%ih'%c.ra.hms.h, lonlat=True, **kwargs)

        for _p in [30,60,-30,-60]:       
        
            # visualization            
            hp.projtext(0-rot[0], _p-rot[1], '%.f$^\circ$'%_p, lonlat=True, **kwargs)

    @staticmethod
    def compute_contours(proportions, samples, nside=64):
        ''' Compute containment contour around desired level.
        '''
        try:  import meander
        except:  return
        
        # binnned map to lower resolution in order to save time
        samples = hp.pixelfunc.ud_grade(samples,nside)
        samples = samples/np.sum(samples)
            
        levels = []
        sorted_samples = list(reversed(list(sorted(samples))))
        nside = hp.pixelfunc.get_nside(samples)
        sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T        
        for proportion in proportions:
            level_index = (np.cumsum(sorted_samples) > \
                           proportion).tolist().index(True)
            level = (sorted_samples[level_index] + \
                     (sorted_samples[level_index+1] \
                      if level_index+1 < len(samples) else 0)) / 2.0
            levels.append(level)
            
        contours_by_level = meander.spherical_contours(sample_points, samples, levels)
        theta_list = {}; phi_list={}
        for cc,contours in enumerate(contours_by_level):            
            _cnt = proportions[cc]
            try:theta_list[_cnt]; phi_list[_cnt]
            except: theta_list[_cnt] = []; phi_list[_cnt] = []
            for contour in contours:            
                theta, phi = contour.T
                phi[phi<0] += 2.0*np.pi
                theta_list[_cnt].append(theta)
                phi_list[_cnt].append(phi)
        return theta_list, phi_list
                   
    def plot_lines(self, ra, dec, hh, ww, n=None, rot=None,
                   colors=['r','k','b','y','g','m','c'], **kwargs):
        '''plot a series of vertices.
        '''
        for _id in range(len(ra)):
            _ra,_dec,_hh,_ww = ra[_id], dec[_id], hh[_id], ww[_id]
            v1_ra,v2_ra,v3_ra,v4_ra,v1_dec,v2_dec,v3_dec,v4_dec = self.vertices(_ra,_dec,_hh,_ww)
            ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra, v1_ra], \
                                         [v1_dec, v2_dec, v4_dec, v3_dec, v1_dec])
            theta, phi = np.pi/2.-np.radians(dec_vertices),np.radians(ra_vertices)
            if not rot is None:                
                r = hp.Rotator(deg=True, rot=rot)
                theta, phi = r(theta, phi)
            if not n is None:
                kwargs['color'] = colors[n[_id] % len(colors)]            
            self.projplot(theta, phi, **kwargs)
            kwargs['label'] = None

    def plot_points(self, ra, dec, n=None, rot=None,
                   colors=['r','k','b','y','g','m','c'], **kwargs):
        '''plot a series of dots.      
        '''                          
        for _id in range(len(ra)):
            _ra, _dec = ra[_id], dec[_id]
            theta, phi = np.pi/2.-np.radians(_dec),np.radians(_ra)            
            if not rot is None:                
                r = hp.Rotator(deg=True, rot=rot)
                theta, phi = r(theta, phi)
            if not n is None:
                kwargs['color'] = colors[n[_id] % len(colors)]  
            self.projplot(theta, phi, **kwargs)
            kwargs['label'] = None        

    def plot_point_data(self, _data, showgrid=True, n=None,                      
                    colors=['r','k','b','y','g','m','c'], **kwargs):
   
        if _data is None:return
        if len(_data['ra']) != len(_data['dec']): return

        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)        
        if showgrid:
            self._grid(**kwargs)            
            self.plot_coord(rot = [kwargs['rot_theta'],
                kwargs['rot_phi'], kwargs['rot_psi']], **kwargs)
            
        # plot galaxies        
        self.plot_points(_data['ra'], _data['dec'], rot=
                [kwargs['rot_theta'], kwargs['rot_phi'],
                 kwargs['rot_psi']], n=n, colors=colors, **kwargs)                

    def plot_box_data(self, _data, showgrid=True, n=None,                      
                      colors=['r','k','b','y','g','m','c'], **kwargs):
        
        if _data is None: return
        for _key in ['ra','dec','fovra','fovdec']:            
            if not _key in _data.keys():return            
            if len(_data[_key]) != len(_data['ra']):return
            
        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)        
        if showgrid:
            self._grid(**kwargs)
            self.plot_coord(rot = [kwargs['rot_theta'],
                    kwargs['rot_phi'], kwargs['rot_psi']], **kwargs)
            
        # plot tilings        
        self.plot_lines(_data['ra'], _data['dec'], _data['fovra'], _data['fovdec'],
                        rot=[kwargs['rot_theta'], kwargs['rot_phi'], kwargs['rot_psi']],
                        n=n, colors=colors, **kwargs)

    def locshow_tel(self, **kwargs):        
        self.pointings.locshow(**kwargs)
            
    def locshow_obs(self, showgrid=True, show_label=False,                 
            colors=['r','k','b','y','g','m','c'], **kwargs):
        
        assert type(self.telescopes) is dict
        for _nt, _tel in enumerate(self.telescopes):
            _data = self.telescopes[_tel]
            if _nt > 0: showgrid=False                                     
            if show_label: kwargs['label'] = '%s' %  _tel
            else: kwargs['label'] = None            
            kwargs['color'] = colors[_nt % len(colors)]
            kwargs['showgrid'] = showgrid            
            _data.locshow(**kwargs)            

    def locshow_all(self, showgrid=True, show_label=False,
            colors=['r','k','b','y','g','m','c'], **kwargs):
        """show tilings in current self.observatories
        """
        assert self.observatories      
        for _no, _obs in enumerate(self.observatories):
            for _nt, _tel in enumerate(self.observatories[_obs]):
                _data = self.observatories[_obs][_tel]
                if _nt > 0 or _no > 0: showgrid=False                                        
                if show_label:
                    kwargs['label'] = '%s tilings (%s)' %  (_tel, _obs)
                else:
                    kwargs['label'] = None
                kwargs['color'] = colors[(_nt+_no) % len(colors)]
                kwargs['showgrid'] = showgrid 
                _data.locshow(**kwargs)                
                
    def plot_route(self, _data, rot=[0,0,0], showgrid=True, **kwargs):
        
        if _data is None: return
        assert 'ra' in _data.keys()
        assert 'dec' in _data.keys()        
        assert len(_data['ra']) == len(_data['dec'])
        
        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)        
        if showgrid:
            self._grid(**kwargs)
            self.plot_coord(rot = rot, **kwargs)
            
        # start point
        _rot = hp.Rotator(deg=True, rot=rot)        
        _theta, _phi = np.pi/2.-np.radians(_data['dec'][0]),\
            np.radians(_data['ra'][0])
        _theta, _phi = _rot(_theta, _phi)
        kwargs['marker'] = 'x'
        self.projplot(_theta,_phi,**kwargs)                
        
        # routes
        theta, phi = np.pi/2.-np.radians(_data['dec']),\
            np.radians(_data['ra'])
        theta, phi = _rot(theta, phi)
        kwargs['marker'] = '.'
        self.projplot(theta,phi,**kwargs) 

    def routeshow_tel(self, rot=[0,0,0], showgrid=True, **kwargs):
        assert self.pointings.data
        assert not self.pointings.data is None
        self.plot_route(self.pointings.data, rot=rot, showgrid=showgrid, **kwargs)

    def routeshow(self, rot=[0,0,0], showgrid=True, **kwargs):
        assert self.data
        assert not self.data is None
        self.plot_route(self.data, rot=rot, showgrid=showgrid, **kwargs)
        
class vtilings(visualization):
    """visualize one tiling network.
       
    Notes
    ----------   
    *vtilings* inherients from **visualization**,
    and inheriented by *tilings*.
    """    
    def locshow(self, data=None, showgrid=True, shown=False,
                colors=['r','k','b','y','g','m','c'], **kwargs):
        """show one tiling network

        Parameters
        ----------
        showgrid :          bool
           show the grids or not          
        Matplotlib.figure Parameters
        Healpix graticule Parameters        
        """
        if data is None:
            data = self.data
        # plot tilings
        if shown:
            assert 'n' in data
            n=data['n']
        else:
            n=None        
        self.plot_box_data(data, showgrid=showgrid, n=n, colors=colors, **kwargs)
        
class vgalaxies(visualization):
    """visualize one galaxy network.
    
    Notes
    ----------   
    *vgalaxies* inherients from **visualization**,
    and inheriented by *galaxies*.
    """
    def locshow(self, data=None, showgrid=True, shown=False,
                colors=['r','k','b','y','g','m','c'], **kwargs):
        """show one galaxy network

        Parameters
        ----------        
        showgrid :          bool
           show the grids or not
        Matplotlib.figure Parameters 
        Healpix graticule Parameters
        """
        if data is None:
            data = self.data
        if shown:
            assert 'n' in data.keys()
            n=data['n']
        else:
            n=None        
        self.plot_point_data(data, showgrid=showgrid, n=n, colors=colors, **kwargs)       

    def hplocshow(self, cls=None, showhp=True, showtype='m', showgrid=True, **kwargs):
        
        if self.hpmapm is None:
            self.logger.info ('Warning: hpmapm not parsed')
            return
        
        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)
                        
        # 2d trigger healpix map
        if showhp:                   
            self.hpmapproj(self.hpmapm, showtype=showtype, **kwargs)
        else:
            self._grid(**kwargs) 
            self.logger.info ('skip showing healpix map')
                    
        # show trigger contours        
        if cls is None:
            self.logger.info ('skip showing contours')            
        
        elif self.is_seq(cls):
            contours = self.compute_contours(cls, self.hpmapm)
            if contours is None:
                self.logger.info ('Warning: no meander installed')
                return
            
            theta_contours, phi_contours = contours
            for _ii, _jj in enumerate(theta_contours):                
                for _theta, _phi in zip(theta_contours[_jj], phi_contours[_jj]):
                    if len(_theta)==0:continue                                       
                    if showhp:
                        self.projplot(_theta,_phi,**kwargs)
                    else:
                        for _kk in ['rot_theta', 'rot_phi', 'rot_psi']:
                            kwargs.setdefault(_kk, 0)                        
                        rot= [kwargs['rot_theta'], kwargs['rot_phi'], kwargs['rot_psi']]
                        r = hp.Rotator(deg=True, rot=rot)
                        _theta, _phi = r(_theta, _phi)
                        self.projplot(_theta,_phi,**kwargs)
                        
        # show grids        
        if showgrid:            
            self._grid(**kwargs)            
            if showhp:
                self.plot_coord(**kwargs)
            else:
                self.plot_coord(rot = [kwargs['rot_theta'],
                    kwargs['rot_phi'], kwargs['rot_psi']], **kwargs)
    
    def distshow(self, histtype='stepfilled',
                 drange=None, nbin=100, **kwargs):
        
        # check if galaxy exists
        if self.data is None: return        
        if not 'dist' in self.data.keys(): return        
        _len = []
        for _key in ['ra','dec','name','dist','mag']:
            _len.append(len(self.data[_key]))
        if len(np.unique(_len)) != 1 or np.unique(_len)[0]==0:
            return

        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)        
                         
        # cut dist
        if not drange is None:
            assert len(drange)==2                    
            dist0 = dist0[np.logical_and(dist0>min(drange), dist0<max(drange))]
        else:
            dist0 = self.data['dist']

        # plot
        kwargs = self.getkeys(kwargs, 'matplotlib')
        plt.hist(dist0,nbin,histtype=histtype,**kwargs)
        plt.xlabel('Distance (Mpc)')
        plt.ylabel('Nuber of galaxies')

    def lumshow(self, nbin=1, color1='k', color2='r', **kwargs):
        # check if galaxy exists
        if self.data is None: return        
        if not 'dist' in self.data.keys(): return        
        _len = []
        for _key in ['ra','dec','name','dist','mag']:
            _len.append(len(self.data[_key]))
        if len(np.unique(_len)) != 1 or np.unique(_len)[0]==0:
            return

        # set key        
        kwargs = self.setkeys(kwargs)
        
        # plot
        fig = self._init_figure(**kwargs)                        
        kwargs = self.getkeys(kwargs, 'matplotlib')        
        Lums = 10**((4.74-self.data['mag'])/2.5)
        dist = self.data['dist']
        ticks,Lcum,Lbin = [],[],[]
        for ii in np.arange(min(dist),max(dist),nbin):
            ticks.append((ii+.5)*nbin)
            Lbin.append(sum(Lums[np.logical_and(dist<ii,dist>ii-nbin)]))
            Lcum.append(sum(Lums[np.logical_and(dist<ii,dist>min(dist))]))            
        plt.plot(ticks,Lbin,drawstyle="steps",\
                label='binned luminosity',color=color1)
        plt.plot(ticks,Lcum,drawstyle="steps",\
                label='cumulative luminosity',color=color2)
        plt.fill_between(ticks,Lbin,step="pre", alpha=1,color=color1)        
        plt.xlabel('Distance (Mpc)')
        plt.ylabel('$L_{sun}$')
        plt.legend()

class vcandidates(visualization):
    """visualize candidates.

    Notes
    ----------   
    *vcandidates* inherients from **visualization**,
    and inheriented by *candidates*.   
    """
    def candshow(self, showgrid=True, **kwargs):
        """show candidates

        Parameters
        ----------        
        showgrid :          bool
           show the grids or not
        Matplotlib.figure Parameters 
        Healpix graticule Parameters
        """
        if 'candidates' in self.__dict__:
            self.plot_point_data(self.candidates, showgrid=showgrid, **kwargs) 

    def lcshow(self, filts=None, xlim=None, ylim=None,ls='',
               vert=None, hori=None, showlim=True, showlegend=True,
               **kwargs): 
        
        if self.lc is None:return
        assert 'time' in self.lc.keys()
        assert 'magnitude' in self.lc.keys()
        assert 'band' in self.lc.keys()

        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)        
                    
        # plot lc                        
        ax = fig.add_subplot(1, 1, 1)
        if not xlim is None: ax.set_xlim(xlim)
        if not ylim is None: ax.set_ylim(ylim)
        ax.invert_yaxis()
        ax.set_xlabel('mjd')
        ax.set_ylabel('mag')
        if not vert is None and not xlim is None:
            for _vert in vert:                                
                ax.plot([_vert, _vert],[min(xlim)-1,max(xlim)+1],'k--')
        if not hori is None and not ylim is None:
            for _hori in hori:
                ax.plot([min(ylim)-1,max(ylim)+1],[_hori, _hori],'k--')
                                                        
        # filters
        if filts is None:
            filts = np.unique(self.lc['band'])

        colors = cm.rainbow(np.linspace(0, 1, len(filts)))
        markers = 'ov^<>*hH+xXDd12348.,spP|_'
        
        for _mm, _filter in enumerate(filts):
                
            # obtain pps and uls
            _pps, _uls = self.pps(pps=True, band=_filter), \
                self.pps(pps=False, band=_filter)
            if len(_pps)==0: continue
                
            # obtain error bars
            _dat = self.errors(pps=True, band=_filter)
            if _dat is None: xerrdat, yerrdat = None, None
            else: xerrdat, yerrdat = _dat
                
            _lmt = self.errors(pps=False, band=_filter)
            if _lmt is None: xerrlim, yerrlim = None, None
            else: xerrlim, yerrlim = _lmt
            
            # data points
            if not _pps is None:
                ax.errorbar(_pps['time'], _pps['magnitude'],
                            xerr=xerrdat, yerr=yerrdat, ls=ls,
                            marker=markers[_mm % len(markers)],
                            color=colors[_mm], label=_filter)
                        
                # upper limits
                if not _uls is None and showlim:
                    ax.errorbar(_uls['time'], _uls['magnitude'],
                                yerr=1, lolims=True, ls=ls,
                                marker=markers[_mm % len(markers)],
                                color=colors[_mm])
        if showlegend: plt.legend()
        return fig
    

    def obsplot(self, _time, _deg, **kwargs):
        '''                
        '''
        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)
        
        # judge        
        if isinstance(_time, astropy.time.Time):
            _time = [_time]       
            _deg = [_deg]        
                                   
        for _x, _y in zip(_time, _deg):
            try: _y=_y.deg
            except:pass
            plt.plot(_x.mjd, _y, **self.getkeys(kwargs, 'projplot'))

    def cooplot(self, _time, _deg, ids=False,
                  colors=['r','k','b','y','g','m','c'], **kwargs):
        '''                
        '''        
        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)
        
        # judge        
        if isinstance(_time, astropy.time.Time):
            _time = [_time]       
            _deg = [_deg]
        
        _x, _y = [], {}                                
        for _ii, (_time0, _deg0) in enumerate(zip(_time, _deg)):            
            _x.append(_time0.mjd)
            for _nn, _deg00 in enumerate(_deg0):                
                if ids: _nn = ids[_nn]                
                try: _y[_nn]
                except: _y[_nn]=[]
                _y[_nn].append(_deg00)
        
        #plot
        for _mm, _m in enumerate(_y):
            kwargs['color'] = colors[_mm % len(colors)]
            if ids: kwargs['label'] = _m            
            plt.plot(_x, _y[_m], **self.getkeys(kwargs, 'projplot'))     
                
class vtrigger(visualization):
    """visualize a trigger map and its contours

    Notes
    ----------   
    *vtrigger* inherients from **visualization**,
    and inheriented by *trigger*.
    """        
    
    def locshow(self, cls=None, showhp=True, showtype='m', showgrid=True, **kwargs):
        """show trigger via healpy.mollview/gnormview....
        parse trigger first so that self.data is not None
        
        Parameters
        ----------                       
        showgrid :          bool
           show the background grids or not
        rot_theta  : int, float
           healpix.visualization rotation angle
           pi/2 - decl, unit in deg, rotate localization map along longtitude direction
        rot_phi    : int, float
           phi = ra, unit in deg, rotate localization map along latitude direction
        rot_psi    :  int, float
           the point at longitude lon and latitude lat will be at the center.   
           An additional rotation of angle psi around this direction is applied.
        num          :   int or  str
           matplotlib/healpy figure number
        showhp :    bool
           show healpix map or not
        showtype  :       string
           if show healpix map, which healpix function was called to use:
           options: [m]ollview, [g]nomview, [c]artview, and [o]rthview
           default: m
        cls:   list
           confidence levels for contours to show
           set to None, will not show contours
           e.g. cls=[.5, .9] would show 50% and 90% CL contours                  
        
        Healpix Parameters                      
        Matplotlib.figure Parameters 
        Healpix graticule Parameters (if showgrid=True)        
        Healpix projplot Parameters (if show contours)
        """
        
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        
        # set key        
        kwargs = self.setkeys(kwargs)
        
        # create figure
        fig = self._init_figure(**kwargs)
                        
        # 2d trigger healpix map
        if showhp:                   
            self.hpmapproj(self.hpmap, showtype=showtype, **kwargs)
        else:
            self._grid(**kwargs) 
            self.logger.info ('skip showing healpix map')
                    
        # show trigger contours        
        if cls is None:
            self.logger.info ('skip showing contours')            
        
        elif self.is_seq(cls):
            contours = self.compute_contours(cls, self.hpmap)
            if contours is None:
                self.logger.info ('Warning: no meander installed')
                return
            
            theta_contours, phi_contours = contours
            for _ii, _jj in enumerate(theta_contours):                
                for _theta, _phi in zip(theta_contours[_jj], phi_contours[_jj]):
                    if len(_theta)==0:continue                                       
                    if showhp:
                        self.projplot(_theta,_phi,**kwargs)
                    else:
                        for _kk in ['rot_theta', 'rot_phi', 'rot_psi']:
                            kwargs.setdefault(_kk, 0)                        
                        rot= [kwargs['rot_theta'], kwargs['rot_phi'], kwargs['rot_psi']]
                        r = hp.Rotator(deg=True, rot=rot)
                        _theta, _phi = r(_theta, _phi)
                        self.projplot(_theta,_phi,**kwargs)
                        
        # show grids        
        if showgrid:            
            self._grid(**kwargs)            
            if showhp:
                self.plot_coord(**kwargs)
            else:
                self.plot_coord(rot = [kwargs['rot_theta'],
                    kwargs['rot_phi'], kwargs['rot_psi']], **kwargs)


            

    ''' '''    
    def cumshow(pparams):

        # arg params
        full=pparams['full']
        fignum=pparams['fignum']
        select=pparams['select']    
        _figsize = (int(pparams['figsize'].split(',')[0]),\
                    int(pparams['figsize'].split(',')[1]))

        # opt params
        try: nameloc=float(pparams['nameloc'])
        except: 
            print ('### Error: nameloc should be a float')
            return
        try: number=int(pparams['number'])
        except: 
            print ('### Error: number should be an interger')
            return
        try: showname=eval(pparams['showname'])
        except: 
            print ('### Error: showname should be True/False')
            return

        # initialize plot
        fig = plt.figure(fignum, figsize=_figsize)
        ax1=pl.axes([.1,.10,.4,.3])
        ax2=pl.axes([.1,.40,.4,.55])
        ax11=pl.axes([.5,.10,.35,.3])
        ax22=pl.axes([.5,.40,.35,.55])
        ax1.set_xlim([0,number+1])
        ax2.set_xlim([0,number+1])
        ax1.set_ylim([10**(-8),.3])
        ax1.set_yticks([10**(-8),.3])
        ax2.set_ylim([.01,1])
        ax1.set_yscale('log')
        ax1.set_xlabel(' '*50+'$N_{gal}$')      
        ax1.set_ylabel('Score')
        ax2.set_ylabel('Cumulative Score')
        ax2.tick_params(labelbottom=False) 
        ax11.set_ylim([10**(-8),.3])   
        ax22.set_ylim([.01,1])   
        ax11.set_yscale('log')
        ax11.set_xscale('log')    
        ax22.set_xscale('log')
        ax22.tick_params(labelbottom=False) 
        ax22.tick_params(labelleft=False) 
        ax11.tick_params(labelleft=False) 

        # plot
        _cc = 0
        _color = ['b','g','y','c','d']
        color1 = 'r'
        color2 = 'k'
        for _nt,_tt in enumerate(full):
       
            ##### full part
            scorei = np.array(pst.decomposit(_tt['score']))
            rai = np.array(pst.decomposit(_tt['ra']))
            deci = np.array(pst.decomposit(_tt['dec']))
            try: namei = np.array(pst.decomposit(_tt['name']))
            except: namei = np.array(pst.decomposit(_tt['ra']))

            # sort score
            idx = np.argsort(np.asarray(scorei))[::-1]

            scorei = scorei[idx]
            rai = rai[idx]
            deci = deci[idx]
            namei = namei[idx]

            # show only interesting fields
            _nm = 5
            score = scorei[:number*_nm]
            ra = rai[:number*_nm]
            dec = deci[:number*_nm]
            name = namei[:number*_nm]

            ax1.plot(np.arange(len(score))+1,score,color1+'.')
            ax2.plot(np.arange(len(score))+1,\
                [sum(score[:y]) for y in range(1, len(score) + 1)],\
                 color1+'.')   
            if showname:
                for i, txt in enumerate(name):
                    ax2.annotate(txt, (i+1, \
                                sum(score[:i+1])+nameloc),\
                                 fontsize=6,rotation=90)

            ax11.plot(np.arange(len(score))+1,score,color1+'.')
            ax22.plot(np.arange(len(score))+1,[sum(score[:y]) \
                                for y in range(1, len(score) + 1)],color1+'.') 
            ax11.set_xlim([number+1,len(score)])    
            ax22.set_xlim([number+1,len(score)])
            ax11.set_xticks([number+1,len(score)])  
            
            ##### for selected dots
            tellist = None
            for _weight in select:
                for _ntt in select[_weight]:
                    if select[_weight][_ntt]['name'] ==\
                       _tt['telescope']['name']:
                        tellist = select[_weight][_ntt]
            if not tellist is None:
                # for one telescope
                _title = '%i pointings with %.2f %% probability'%\
                    (number, sum(score[:number])*100)
                for _ra,_dec in zip(tellist['ra'],tellist['dec']):
                    s,t = np.where(rai == _ra),\
                        np.where(deci == _dec)
                    ids = np.intersect1d(s, t)                
                    ax1.plot(ids+1,scorei[ids],_color[_cc]+'x')
                    ax2.plot(ids+1,np.cumsum(scorei)[ids],_color[_cc]+'x')
                _title += '\n%s: %s'%\
                    (tellist['name'],_color[_cc])
                _cc+=1
        plt.title(_title)
        return fig
