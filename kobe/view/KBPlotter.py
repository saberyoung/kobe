#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/view/PstPlotter.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import logging
import astropy.coordinates
import astropy.time
import astropy.units as u

__all__ = ('PstPlotter')

class PstPlotter():
    """PstGetGalaxies: Generate, read and process galaxies

    Parameters
    ----------
    defconf :    `dict`           
      default configure, if any PstGetGalaxies parameter was included in defconf dictionary, 
      then its default value would be applied
    logger :     `class`
      logging object
    """        

    # Static version info
    version = 1.0
        
    def __init__(self, defconf=None, logger = None):
        """
        initialize
        """

        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # ----- define default parameters ----- #
        self.run_config(defconf)
        
        self.colors_contour = ['b','g','k','y','c','m']
        self.colors_field = ['b','g','k','y','c','m']       
    
    def run_config(self, defconf):
        """read default value for parameters

        Parameters
        ----------
        defconf :      `dict`
          input parameter settings, if None, use default

        Examples
        --------    
        >>> from pst.view.PstPlotter import PstPlotter
        >>> a = PstPlotter()
        >>> dict = {'theta': '90', 'obstime': 3600*2}
        >>> a.run_config(dict)
        """
         
        self.conf = {
            'contours':    [.5,.9],     # 2D contours to show for triggers
            'obstime':     None,        # observing time
            'theta':       0,           # longtitude, theta = pi/2 - decl, unit in deg
            'phi':         0,           # latitude, phi = ra, unit in deg
            'psi':         0,           # 
            'nest':        False,       # healpix ordering
            'cbar':        False,       # healpix color bar
            'norm':        'None',      # normalization, options:
                                        # hist: histagram / log: logarithmic / None: linear
            'coord':       'C',         # coordinate system, options: C, E, G               
            'vmin':        0,           # minimum range value for 2d normalized map
            'vmax':        1e-4,        # maximum range value for 2d normalized map
            'ptype':       'm',         # healpix plot type
            'nside':       64,          # healpix resolution for contour view
            'figsize':     (15, 10)     # default figure size
        }
        
        if defconf is None: return        
        for _k in self.conf.keys():
            if _k in defconf.keys():
                self.conf[_k] = defconf[_k]
            else:
                self.logger.info ('### Warning: use default value for %s'%_k)
    
    def locshow(self, psttrigger=None, psttilings=None, pstgalaxies=None,
                tellist=None, theta=None, phi=None, psi=None, figsize=None,
                contours=None, nest=None, coord=None, norm=None, vmin=None,
                vmax=None, fignum=1, ptype=None, cbar=None, obstime=None,
                nside=None, title='sky localization'):
        """show 2d localization plot for trigger, galaxies, tilings

        Parameters
        ----------
        psttrigger :       `class`
           PstParseTriggers object
        psttilings :       `class`
           PstGetTilings object
        pstgalaxies :      `class`
           PstGetGalaxies object
        tellist :          sequence of `class`
           list of PstGetTilings obj or PstGetGalaxies obj
           to be used to define the observatory so that plot the horizon lines        
        theta :            `class` between 0 and 360
           rotate along the healpix theta directions, unit in deg
        phi :              `class` between 0 and 360
           rotate along the healpix phi directions, unit in deg
        psi :              `class` between 0 and 360
           rotate along the healpix psi directions, unit in deg
        figsize :          `tuple`
           define the size of figures
        contours :         `list`
           list of confidence level, e.g. [.5, .9]
        nest :             `bool`
           healpix ordering options: 
           if True, healpix map use `nest` ordering, otherwise, use `ring` instead 
        coord :             sequence of character
           Either one of ‘G’, ‘E’ or ‘C’ to describe the coordinate system of the map, 
           or a sequence of 2 of these to rotate the map from the first to the 
           second coordinate system.
        norm :              `hist`, `log`, None
           Color normalization, hist= histogram equalized color mapping, log= logarithmic color mapping, 
           default: None (linear color mapping)
        vmin :              `float`
           The minimum range value for mollview/gnomview/etc
        vmax :              `float`
           The maximum range value for mollview/gnomview/etc
        fignum :            `int`
           The figure number of generated matplotlib plot
        ptype :             `string`
           type of trigger plots: `m` = mollview, `g` = gnomview, `c` = cartview, and `o` = orthview
        cbar :              `bool`
           Display the colorbar
        obstime :           `float` or None or astropy.time obj
           observing time.
           None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; 
           `astropy.time` set specific observing time  
        nside :             `int`
           healpix nside parameter, must be a power of 2, less than 2**30
        title :             `string`
           figure title of the plot

        Examples
        --------    
        >>> from pst.view.PstPlotter import PstPlotter
        >>> a = PstPlotter()       
        """
        
        # create figure
        if not figsize is None: self.conf['figsize'] = figsize               
        try:
            fig = plt.figure(fignum, figsize=self.conf['figsize'])
        except:
            self.logger.info ('### Error: wrong figsize or fignum')
            return        

        # plot
        if not psttrigger is None:
            self.triggershow(psttrigger, theta=theta, phi=phi, psi=psi, title=title,
                        contours=contours, nest=nest, coord=coord, norm=norm, vmin=vmin,
                        vmax=vmax, fignum=fignum, ptype=ptype, cbar=cbar, nside=nside)

        hp.graticule(dpar=None, dmer=None, coord=coord) 
        '''
        if not psttilings is None:
            self.tilingshow(psttilings)
                    
        if not pstgalaxies is None:
            self.galaxyshow(pstgalaxies)

        # others
        self.notes()
        
        # Set grid lines
        hp.graticule(dpar=None, dmer=None, coord=coord)
        plt.legend()
        '''
        
    def triggershow(self, psttrigger, contours=None, nside=None,
            theta=None, phi=None, psi=None, nest=None, coord=None,
            norm=None, vmin=None, vmax=None, fignum=1, ptype=None,
            cbar=None, title='sky localization'):
        """ 
        2d trigger healpix map
        """                
        from pst.cookbook import is_seq, is_seq_of_seq
        if is_seq_of_seq(psttrigger.data['hpmap']):
            (hpx, hpd1, hpd2, hpd3) = psttrigger.data['hpmap']
        elif is_seq(psttrigger.data['hpmap']):
            hpx = psttrigger.data['hpmap']
        else:
            self.logger.info ('### Warning: failed to load healpix map')
            return
        
        # read parameters        
        if theta is None:    theta    =  float(self.conf['theta'])
        if phi is None:      phi      =  float(self.conf['phi'])
        if psi is None:      psi      =  float(self.conf['psi'])
        if nest is None:     nest     =  self.conf['nest']
        if cbar is None:     cbar     =  self.conf['cbar']    
        if coord is None:    coord    =  self.conf['coord']
        if norm is None:     norm     =  self.conf['norm']        
        if vmin is None:     vmin     =  self.conf['vmin']
        if vmax is None:     vmax     =  self.conf['vmax']        
        if ptype is None:    ptype    =  self.conf['ptype'] 
        
        # Choose color map and set background to white
        cmap = cm.YlOrRd
        cmap.set_under("w")

        # Plot GW skymap
        if ptype == 'm':
            hp.mollview(map=hpx, fig=fignum, rot=[theta,phi,psi], coord=coord, \
                    unit=r'Probability', xsize=800, title=title, \
                    nest=nest, min=vmin, max=vmax, flip='astro', \
                    remove_dip=False, remove_mono=False, gal_cut=0, \
                    format='%g', format2='%g', cbar=cbar, cmap=cmap, \
                    notext=False, norm=norm, hold=True, margins=None, \
                    sub=None, return_projected_map=False)      
        elif ptype == 'g':
            hp.gnomview(map=hpx, fig=fignum, rot=[theta,phi,psi], coord=coord,\
                    unit='', xsize=5000, ysize=None, reso=1.5, title=title, \
                    nest=nest, remove_dip=False, remove_mono=False, gal_cut=0, \
                    min=vmin, max=vmax, flip='astro', format='%.3g', cbar=False, \
                    cmap=cmap, badcolor='gray', bgcolor='white', norm=norm, \
                    hold=True, sub=None, margins=None, notext=False, \
                    return_projected_map=False, no_plot=False)
        elif ptype == 'c':
            hp.cartview(map=hpx, fig=fignum, rot=[theta,phi,psi], \
                    zat=None, coord=coord, unit='', xsize=800, ysize=None, \
                    lonra=None, latra=None, title=title, nest=nest, \
                    remove_dip=False, remove_mono=False, gal_cut=0, \
                    min=vmin, max=vmax, flip='astro', format='%.3g', \
                    cbar=False, cmap=cmap, badcolor='gray', bgcolor='white', \
                    norm=norm, aspect=None, hold=False, sub=None, margins=None, \
                    notext=False, return_projected_map=False)
        elif ptype == 'o':
            hp.orthview(map=hpx, fig=fignum, rot=[theta,phi,psi], coord=coord, \
                    unit='', xsize=800, half_sky=False, title=title, nest=nest, \
                    min=vmin, max=vmax, flip='astro', remove_dip=False, \
                    remove_mono=False, gal_cut=0, format='%g', format2='%g', \
                    cbar=False, cmap=cmap, badcolor='gray', bgcolor='white', \
                    notext=False, norm=norm, hold=False, margins=None, sub=None, \
                    return_projected_map=False)
        else:
            self.logger.info ('### Error: wrong option for ptype [m, g, c, o]')
            return

        """ 
        show trigger contours
        """               
        # read parameters                
        if contours is None: contours =  self.conf['contours']
        if nside is None:    nside    =  self.conf['nside']
        if is_seq(contours):
            theta_contour, phi_contour = self.compute_contours(contours,hpx,nside=nside)
            for ndd,dd in enumerate(theta_contour):                
                _theta_contour, _phi_contour = theta_contour[dd], phi_contour[dd]        
                for i in range(len(_theta_contour)):                    
                    if len(_theta_contour[i])==0:continue                    
                    hp.projplot(_theta_contour[i],_phi_contour[i], coord=coord, linewidth=1,
                                c=self.colors_contour[ndd])

    def tilingshow(self, psttilings, theta=None, phi=None, psi=None,
                   color=None, label=None, nest=None, coord=None,
                   fignum=1):

        n, ra, dec, fovra, fovdec = psttilings.data['n'], psttilings.data['ra'], \
            psttilings.data['dec'], psttilings.data['fovra'], psttilings.data['fovdec']
        if len(ra) >0: pass
        else:
            self.logger.info ('### Warning: skipping tilingshow since tilings are not parsed')
            return
        
        if color is None:    color    =  self.colors_field[0]
        if label is None:    label    =  psttilings.name
        if coord is None:    coord    =  self.conf['coord']
        if fignum is None:   fignum   =  1               

        fig = plt.figure(fignum)               
        self.plot_lines(ra, dec, fovra, fovdec,
                        coord=coord, color=color,label=label)        
        return fig        
    
    def notes(self, psttrigger=None, tellist=None,
              coord=None, obstime=None, fignum=1,
              contours=None):
        """
        show specific coordiantes in healpy plots
        """        
        self.plot_coord(coord)        

        """
        plot the sky: 
        sun, moon, horizon, galactic plane,
        """
        if obstime is None:    obstime    =  self.conf['obstime'] 
        obstime = self.obstime(obstime)
        self.plot_sky(tellist,coord,obstime,fignum)

        # write labels
        xx,yy = -2.,1.
        plt.text(xx, .9, 'sun: $\odot$',fontsize=20,\
                 ha="center", va="center", color='y')
        plt.text(xx, 1., 'moon: $\oplus$',fontsize=20,\
                 ha="center", va="center", color='b')        

        # contou informations
        from pst.cookbook import is_seq, is_seq_of_seq
        if is_seq(contours):
            area = psttrigger.calc_area(contours)            
            for ii,cc in enumerate(contours):
                yy+=.1
                info = '%.2f%% area: %.2f deg^2' % (cc*100, area[cc])
                plt.text(xx, yy, info, ha="center", va="center",
                         color=self.colors_contour[ii], fontsize=20)        
   
    def galaxies(self, pstgalaxies, nest=None, coord=None, colors=None, fignum=1):
        """ galaxy dot plots:
        pstgalaxies: 
        1. PstGetGalaxies obj
        2. PstGetGalaxies obj list
        color: list
        """
        from pst.cookbook import is_seq, is_seq_of_seq
        if coord is None:    coord    =  self.conf['coord']
        if nest is None:     nest     =  self.conf['nest']        
        if colors is None:
            colors = self.colors_field
        else:
            if not is_seq(colors):
                self.logger.info('### Error: colors should be a list')
                return
            if len(colors) < len(pstgalaxies):
                self.logger.info('### Error: not enough colors')
                return 
        if is_seq(pstgalaxies):
            for pstgalaxy, color in zip(pstgalaxies, colors):
                self.galaxy(pstgalaxy, nest=nest, coord=coord,
                            color=color, fignum=fignum)
        else:
            self.galaxy(pstgalaxies, nest=nest, coord=coord,
                        color=colors[0], fignum=fignum)

    def galaxy(self, pstgalaxy, nest=None, coord=None, color=None, fignum=1):

        if coord is None:    coord    =  self.conf['coord']
        if nest is None:     nest     =  self.conf['nest']     
        if pstgalaxy.data['ra'] is None or pstgalaxy.data['dec'] is None:                    
            self.logger.info('### Warning: fields not found, skipped for %s'%pstgalaxy.name)
            return
        ra  =  pstgalaxy.data['ra']
        dec =  pstgalaxy.data['dec']        
        ms=4

        # switch figure
#        plt.figure(fignum)

        # Set grid lines
        hp.graticule(dpar=None, dmer=None, coord=coord)

        # points plot
        theta, phi = np.pi/2.-np.radians(dec),np.radians(ra)                
        hp.projplot(theta[0],phi[0],'x', color =color,
                    coord=coord, ms = ms, label='%s: %i galaxies'%(pstgalaxy.name, len(ra)))
        hp.projplot(theta,phi,'x', color =color,
                    coord=coord, ms = ms)
        plt.legend()

    @staticmethod
    def plot_coord(coord):

        for _t in [60,120,180,240,300,360]:
            # deg to hms
            c= astropy.coordinates.SkyCoord(ra=_t*u.degree,
                                dec=0*u.degree, frame='icrs')          
    
            # visualization
            hp.projtext(_t, 0, '%ih'%c.ra.hms.h, lonlat=True, coord=coord)

        for _p in [30,60,-30,-60]:       
        
            # visualization
            hp.projtext(0, _p, '%.f$^\circ$'%_p, lonlat=True, coord=coord)

    @staticmethod
    def plot_sky(tellist,coord,obstime=None,fignum=1):
        '''
        tellist: list of PstGetTilings obj or PstGetGalaxies obj
        '''
        
        _colorlist = ['b','g','k','y','c','m']
#        plt.figure(fignum)

        if not tellist is None:
            # plot the horizon for each telescopes
            for nt, tel in enumerate(tellist):
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

        # plot the galactic plane        
        ra = np.arange(0,360,10)
        dec = np.zeros(len(ra))
        _smlabel=True
        for _ra,_dec in zip(ra,dec):
        
            # from galactic coordinates to equatorial
            _radecs = astropy.coordinates.SkyCoord(l=_ra*u.deg,
                                        b=_dec*u.deg, frame='galactic')                

            # transform to theta phi
            _ra,_dec = _radecs.icrs.ra.deg, _radecs.icrs.dec.deg
            _theta,_phi = np.pi/2.-np.radians(_dec),np.radians(_ra)            

            # plot
            if _smlabel:
                hp.projplot(_theta,_phi,'x', color = 'k', 
                            coord=coord, ms = 10, label='galactic plane')
                _smlabel=False
            else:
                hp.projplot(_theta,_phi,'x', color = 'k', 
                            coord=coord, ms = 10)

        # plot the sun
        _sra,_sdec = astropy.coordinates.get_sun(obstime).ra.deg,\
            astropy.coordinates.get_sun(obstime).dec.deg
        _stheta,_sphi = np.pi/2.-np.radians(_sdec),np.radians(_sra)        
        hp.projtext(_stheta,_sphi,'$\odot$', color = 'y',
                    coord=coord,fontsize=20)     

        # plot the moon lasting 5 days
        for _nd in np.arange(0,5,1):
            _timeplus = _nd*24
            timelater = obstime + astropy.time.TimeDelta(_timeplus*3600, format='sec')
            _mra,_mdec = astropy.coordinates.get_moon(timelater).ra.deg,\
                astropy.coordinates.get_moon(timelater).dec.deg
            _mtheta,_mphi = np.pi/2.-np.radians(_mdec),np.radians(_mra)            
            hp.projtext(_mtheta,_mphi,'$\oplus$', color = 'b',
                        coord=coord,fontsize=20)
            hp.projtext(_mtheta,_mphi-.1,'%id'%_nd, color = 'b', \
                        coord=coord,fontsize=12)        

    @staticmethod
    def obstime(t=None):
        """ define observing time
        t: None, or float, unit in sec, or astropy.time
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

    @staticmethod
    def compute_contours(proportions, samples, nside=64):
        r''' Plot containment contour around desired level.
        E.g 90% containment of a PDF on a healpix map
        '''
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

        import meander        
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

    def distview(pparams):
        
        _distmin=float(pparams['distmin'])
        _distmax=float(pparams['distmax'])
        dist0=pparams['dist']
        color1=pparams['color1']
        color2=pparams['color2']    
        try:nbin=int(pparams['nbin'])
        except:
            return ('!!! Warning: nbin should be an integer')            
        _figsize = (int(pparams['figsize'].split(',')[0]),\
                    int(pparams['figsize'].split(',')[1]))    
        label=pparams['label']
        fignum=pparams['fignum']
        scale=pparams['scale']

        # cut dist
        dist0 = dist0[np.logical_and(dist0>_distmin,
                                 dist0<_distmax)]
        fig = plt.figure(fignum, _figsize)
        for nii,ii in enumerate(np.arange(_distmin,_distmax,abs(_distmax-_distmin)/nbin)):                 
            cum = len(dist0[np.logical_and(dist0<ii,dist0>_distmin)])
            pl.text(ii,-2,cum,color=color1)
        pl.hist(dist0,nbin,label=label,histtype='stepfilled',color=color2)
        if scale=='log': plt.yscale('log')
        plt.xlabel('Distance (Mpc)')
        plt.ylabel('Nuber of galaxies')
        plt.legend()
        return fig

    def lumsview(pparams):
        
        _distmin=float(pparams['distmin'])
        _distmax=float(pparams['distmax'])    
        dist0=pparams['dist']   
        mag0=pparams['mag']
        color1=pparams['color1']  
        color2=pparams['color2']      
        nbin=int(pparams['nbin'])   
        label=pparams['label']
        fignum=pparams['fignum']
        scale=pparams['scale']    
        _figsize = (int(pparams['figsize'].split(',')[0]),\
                    int(pparams['figsize'].split(',')[1]))
        fig = plt.figure(fignum, _figsize)
        Lums = 10**((4.74-mag0)/2.5)
        ticks,Lcum,Lbin = [],[],[]
        for ii in np.arange(_distmin,_distmax,nbin):
            ticks.append((ii+.5)*nbin)
            Lbin.append(sum(Lums[np.logical_and(dist0<ii,dist0>ii-nbin)]))
            Lcum.append(sum(Lums[np.logical_and(dist0<ii,dist0>_distmin)]))            
        pl.plot(ticks,Lbin,drawstyle="steps",\
                label='binned luminosity',color=color1)
        pl.plot(ticks,Lcum,drawstyle="steps",\
                label='cumulative luminosity',color=color2)
        pl.fill_between(ticks,Lbin,step="pre", alpha=1,color=color1)
        pl.title(label)   
        if scale=='log': plt.yscale('log')
        plt.xlabel('Distance (Mpc)')
        plt.ylabel('$L_{sun}$')
        plt.legend()     
        return fig

    def routeview(pparams):
        
        ralist=pparams['ra']
        declist=pparams['dec']
        _time=pparams['time']
        color=pparams['color']
        label=pparams['label']
        rot_phi=float(pparams['phi'])
        rot_theta=float(pparams['theta'])
        fignum=pparams['fignum']
        coord=pparams['coord']
        ms = 4

        if len(ralist)>0:pass
        else:return

        fig = plt.figure(fignum)
        hp.graticule()

        # start point
        _rot = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])
        _theta,_phi = pst.RadecToThetaphi(ralist[0],declist[0])
        hp.projplot(_rot(_theta,_phi),'*', color=color, \
                    coord=coord, ms = ms)

        # routes
        theta,phi = pst.RadecToThetaphi(ralist,declist)
        hp.projplot(_rot(theta,phi),color=color,\
                    coord=coord) #,label=label)

        plt.legend()
        return fig

    def plot_lines(self, ra, dec, hh, ww, color='k',coord='C',label=None):
        '''
        '''
        
        from pst.pipeline.PstGetTilings import PstGetTilings        
        if not label is None: _smlabel=True
        else: _smlabel=False
        for _ra,_dec,_hh,_ww in zip(ra,dec,hh,ww):
            v1_ra,v2_ra,v3_ra,v4_ra,v1_dec,v2_dec,v3_dec,v4_dec = PstGetTilings.vertices(_ra,_dec,_hh,_ww)
            ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra, v1_ra], \
                                         [v1_dec, v2_dec, v4_dec, v3_dec, v1_dec])
            theta, phi = np.pi/2.-np.radians(dec_vertices),np.radians(ra_vertices)
            if _smlabel:
                hp.projplot(theta, phi, color=color, coord=coord,label=label)
                _smlabel=False
            else:
                hp.projplot(theta, phi,coord=coord, color=color)

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

def interactive_show(func,_pm,_opt):   
    """
    Usage: intractive(plotfunc,dict,list)
    """                
    answ = False
    while not answ:
        _fig = func(_pm)   
        _show = '    current options:'
        for ii in _opt:
            try: _show += ' %s=%.2e '%(ii,float(_pm[ii])) 
            except: _show += ' %s=%s '%(ii,_pm[ii])
        print(_show)
        answ = True
        good_answ=False
        while not good_answ: 
            good_answ = True
            try:
                _answ = input('>>> good options (y) or quit (q) or change'+\
                          ' options (eg. %s=%.2f) [y]: '%\
                          (_opt[0],float(_pm[_opt[0]])))    
            except:
                print ('!!! Error: parameter wrong')
                _answ = ''
                good_answ = False
            if len(_answ) > 0: 
                if _answ=='y': answ= True
                elif _answ=='q':sys.exit('quit...')
                else:
                    try:
                        if _answ[:_answ.index('=')] in _opt:
                            _pm[_answ[:_answ.index('=')]] = \
                                _answ[_answ.index('=')+1:].strip() 
                            try:_fig.clf()
                            except:print ('!!! Warning: no return from inter func')
                        else: 
                            print("!!! Warning; option "+_answ+\
                                     " not available !!!" )
                            good_answ = False
                    except: 
                        print("!!! Error: wrong input: retry")
                        good_answ = False
                    answ = False
    return _pm,_fig
