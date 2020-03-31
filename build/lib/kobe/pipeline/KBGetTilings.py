#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/KBGetTilings.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import numpy as np
import healpy as hp
import os
import random
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
import logging

__all__ = ('KBGetTilings')

class KBGetTilings():
        """KBGetTilings: Generate, read and process tilings

        Parameters
        ----------
        name :       `string`
          telescope name, if none, will be constructed with coordinates
        ra :         `list`           
          tiling ra list, default: None
        dec :        `list`           
          tiling dec list, default: None
        fovra :      `list`           
          tiling field of view list in ra direction, 
             default: None
        fovdec :     `list`           
          tiling field of view list in dec direction,
             default: None
        index :      `list`           
          tiling index list, default: None
        defconf :    `dict`           
          default configure, if any KBGetTilings parameter was included in defconf dictionary, 
          then its default value would be applied
        logger :     `class`
          logging object

        See Also
        --------
        KBParseTriggers, KBGetGalaxies

        Examples
        --------
        see also https://github.com/saberyoung/kobe/blob/master/notebook/test_tilings.ipynb
        
        >>> from kobe.pipeline.KBGetTilings import KBGetTilings
        >>> a = KBGetTilings()
        >>> a.data
        """    

        # Static version info
        version = 1.0
   
        def __init__(self, name=None, ra=None, dec=None,
                     fovra=None, fovdec=None,
                     index=None, defconf=None,
                     logger=None):
                """                                        
                initialize        
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set tiling properties ----- #
                self.data   =   {'ra':ra, 'dec':dec,'n':index,
                                 'fovra':fovra, 'fovdec':fovdec}

                # ----- define default parameters ----- #
                self.run_config(defconf)

                # ----- set name for telescope ----- #
                self.set_name(name)
                
        def set_name(self, name):
                """set telescope name

                Parameters
                ----------
                name :      `string` 
                  if not telescope name given, construct with telescope coordinate:
                  lon-lat-alt
                """
                if name is None:
                        self.name = '%.2f-%.2f-%i' % (self.conf['lon'],
                                self.conf['lat'], self.conf['alt'])
                else:
                        self.name = name
                
        def run_config(self, defconf):
                """read default value for parameters

                Parameters
                ----------
                defconf :      `dict`
                  input parameter settings, if None, use default

                Examples
                --------    
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> dict = {'limdec': [-20,90]}
                >>> a.run_config(dict)
                """
                      
                self.conf = {
                        'limra':      [0,360.],  # `range`   [0,360.]
                        'limdec':     [-89,89],  # `range`   [-89.,89.]
                        'fovra':      1.,        # `float`   > 0
                        'fovdec':     1.,        # `float`   > 0
                        'shiftra':    0.,        # `float`   abs(shiftra) < fovra
                        'shiftdec':   0.,        # `float`   abs(shiftdec) < fovdec
                        'obra':       1,         # `int`
                        'obdec':      1,         # `int`
                        'skipfile':   None,      # `str`     options: `txt`, `npz`
                        'skipfrac':   0.,        # `float`   [0, 1]
                        'nside':      512,       # `int`     default healpix resolution
                        'nest':       False,     # `bool`  healpix map defualt ordering: nest or ring
                        'wdir':       '/tmp/',   # `str`     working directory
                        'filetype':   'npz',     # `str`     options: `txt`, `npz`
                        'filename':   'kobe_tilings', # `str`     tilings file name
                                                     #           relative path, use wdir set directory
                                                     #           without suffix
                        'obstime':    None,          #
                        'lat':        24.625,        # la parma
                        'lon':        70.403,        # la parma
                        'alt':        2635           # la parma
                }
        
                if defconf is None: return        
                for _k in self.conf.keys():
                        if _k in defconf.keys():
                                self.conf[_k] = defconf[_k]
                        else:
                                self.logger.info ('### Warning: use default value for %s'%_k)
                
        def generate(self,limra=None,limdec=None,
                     fovra=None,fovdec=None,
                     shiftra=None,shiftdec=None):
                """create pointings by tiling sky

                Parameters
                ----------
                limra :         `range`
                  tilings: ra range, default: [0, 360]
                limdec :        `range`           
                  tilings: dec range, default: [-90, 90]
                fovra :         `float`           
                  tiling field of view in ra direction
                fovdec :        `float`           
                  tiling field of view in dec direction
                shiftra :       `float`                  
                  ra of initial tiling will be 0 + shiftra
                shiftdec :      `float`           
                  dec of initial tiling will be 0 + shiftdec

                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])
                >>> a.data
                {'n': array([    0,     1,     2, ..., 27788, 27789, 27790]), 'ra': array([  1.06418,   2.12836,   3.19253, ..., 286.53708, 315.19079, 343.8445 ]), 'dec': array([-20., -20., -20., ...,  88.,  88.,  88.]), 'fovra': array([1., 1., 1., ..., 1., 1., 1.]), 'fovdec': array([1., 1., 1., ..., 1., 1., 1.])}
                >>> a.astrotab()
                <Table length=27791>
                n       ra      dec    fovra   fovdec
                int64  float64  float64 float64 float64
                ----- --------- ------- ------- -------
                0   1.06418   -20.0     1.0     1.0
                1   2.12836   -20.0     1.0     1.0
                2   3.19253   -20.0     1.0     1.0
                3   4.25671   -20.0     1.0     1.0
                4   5.32089   -20.0     1.0     1.0
                5   6.38507   -20.0     1.0     1.0
                6   7.44924   -20.0     1.0     1.0
                7   8.51342   -20.0     1.0     1.0
                8    9.5776   -20.0     1.0     1.0
                9  10.64178   -20.0     1.0     1.0
                10  11.70596   -20.0     1.0     1.0
                11  12.77013   -20.0     1.0     1.0
                12  13.83431   -20.0     1.0     1.0
                13  14.89849   -20.0     1.0     1.0
                ...       ...     ...     ...     ...
                27776 305.71716    87.0     1.0     1.0
                27777 324.82448    87.0     1.0     1.0
                27778 343.93181    87.0     1.0     1.0
                27779  28.65371    88.0     1.0     1.0
                27780  57.30742    88.0     1.0     1.0
                27781  85.96113    88.0     1.0     1.0
                27782 114.61483    88.0     1.0     1.0
                27783 143.26854    88.0     1.0     1.0
                27784 171.92225    88.0     1.0     1.0
                27785 200.57596    88.0     1.0     1.0
                27786 229.22967    88.0     1.0     1.0
                27787 257.88338    88.0     1.0     1.0
                27788 286.53708    88.0     1.0     1.0
                27789 315.19079    88.0     1.0     1.0
                27790  343.8445    88.0     1.0     1.0
                """
                
                _hp = self.checkdata()
                if _hp: return
                
                if not limra is None: self.conf['limra'] = limra
                if not limdec is None: self.conf['limdec'] = limdec
                if not fovra is None: self.conf['fovra'] = fovra
                if not fovdec is None: self.conf['fovdec'] = fovdec
                if not shiftra is None: self.conf['shiftra'] = shiftra
                if not shiftdec is None: self.conf['shiftdec'] = shiftdec                
        
                # limit shift to proper range
                if abs(self.conf['shiftra']) >= self.conf['fovra']:
                        self.logger.info ('# Warning: abs(shiftra) < fovra')
                        return
                
                if abs(self.conf['shiftdec']) >= self.conf['fovdec']:
                        self.logger.info ('# Warning: abs(shiftdec) < fovdec')
                        return
                
                # cut ra,dec
                ramin, ramax = min(self.conf['limra'])+self.conf['shiftra'], \
                        max(self.conf['limra'])+self.conf['shiftra']
                decmin, decmax = min(self.conf['limdec'])+self.conf['shiftdec'], \
                        max(self.conf['limdec'])+self.conf['shiftdec']
                ramin = max(ramin, 0)
                ramax = min(ramax, 360)
                decmin = max(decmin, -89)
                decmax = min(decmax, 89)
                
                # dec range
                decrange= np.arange(decmin,decmax,self.conf['fovdec'])

                # get network, i.e. ra,dec list
                ralist,declist=[],[]
                fovralist,fovdeclist=[],[]
                for _dec in decrange:
                        npoint = 360*np.cos(_dec*np.pi/180)/self.conf['fovra']
                        for nn in np.arange(0,npoint,1):  
                                _ra = 360./npoint*nn+self.conf['shiftra']
                                if _ra < ramax and _ra > ramin:
                                        ralist.append(float('%.5f'%_ra))
                                        declist.append(float('%.5f'%_dec))
                                        fovralist.append(float('%.5f'%self.conf['fovra']))
                                        fovdeclist.append(float('%.5f'%self.conf['fovdec']))
                                        
                self.data = {'n':       np.arange(len(ralist)),
                             'ra':      np.array(ralist),
                             'dec':     np.array(declist),                             
                             'fovra':   np.array(fovralist),
                             'fovdec':  np.array(fovdeclist)}

        def generate_mc(self, skymap, num,
                        limra=None,limdec=None,
                        fovra=None,fovdec=None):
                """monte carlo approach on shiftra and shiftdec to maxmize trigger probability with a number of pointings

                Parameters
                ----------
                skymap :        `class`
                  KBParseTriggers object
                num :           `int`           
                  number of pointings
                limra :         `range`
                  tilings: ra range, default: [0, 360]
                limdec :        `range`           
                  tilings: dec range, default: [-90, 90]
                fovra :         `float`           
                  tiling field of view in ra direction
                fovdec :        `float`
                  tiling field of view in dec direction
                
                """
                _hp = self.checkdata()
                if _hp: return
                
                if not limra is None: self.conf['limra'] = limra
                if not limdec is None: self.conf['limdec'] = limdec
                if not fovra is None: self.conf['fovra'] = fovra
                if not fovdec is None: self.conf['fovdec'] = fovdec                

                # monte carlo for tiling
                _log, _nloop = [0.], 3
                shifthi, shiftwi = 0, 0
                for nn in [5.,10.,20.]:
                        if verbose: 
                                print(' - searching in fovh/%i fovw/%i'%(nn,nn))
                        shifth, shiftw=fovh/nn, fovw/nn

                        nloop = 0
                        answ1 = False
                        while not answ1:  # if angle OK: loop 100 times
                                angle = random.uniform(0,2*np.pi)
                                if verbose: print('\t %i with angle: %.2f'%(nloop,angle))
                                _shifth,_shiftw = np.sqrt(shifth**2+shiftw**2)*np.sin(angle),\
                                        np.sqrt(shifth**2+shiftw**2)*np.cos(angle) 
            
                                answ2 = False
                                while not answ2:  # if OK, go on, if no, change
                                        shifthi += _shifth
                                        shiftwi += _shiftw
                                        
                                        # generate pointings
                                        _ral,_decl = kobe.gen_pointings(limdec=limdec,limra=limra,\
                                                fovh=fovh,fovw=fovw,shifth=shifthi,shiftw=shiftwi)
                                        # cal prob for tiling list
                                        t = kobe.calprob_tile(skymap,_ral,_decl,fovh,fovw)  

                                        # judge converge or not
                                        if sum(t)>_log[-1]:  # accept, direction correct
                                                _log.append(sum(t))
                                                print ('\t\tcovered %.5e probs'%_log[-1])
                                        else:  # reject, change direction
                                                nloop+=1
                                                answ2 = True
                                                if nloop>=_nloop: answ1=True
                                        _ral,_decl = kobe.gen_pointings(limdec=limdec,limra=limra,\
                                                fovh=fovh,fovw=fovw,shifth=shifthi,shiftw=shiftwi)
                                        _ral,_decl = kobe.remove_fields(skipfile,_ral,_decl,\
                                                                np.zeros(len(limra))+fovw,\
                                                                np.zeros(len(limra))+fovh,verbose)
                                        


        def checkdata(self):
                """check if KBGetTilings is parsed or not
        
                returns
                ----------
                res :        `bool`                  
                """
                 
                if not self.data['ra'] is None and \
                   not self.data['dec'] is None and \
                   not self.data['fovra'] is None and \
                   not self.data['fovdec'] is None:
                        self.logger.info ('tilings has already been parsed')
                        return True
                else:
                        return False
                
        def read(self, filename=None, filetype=None, wdir=None):
                """read pointings from file

                Parameters
                ------------
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file, options: `txt` (plain text), `npz` (pickle file)
                wdir :             `string`
                  working directory for file

                Examples
                --------        
        
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.read(filename='kobe_tilings', filetype='npz', wdir='./')
                """
                
                _res = self.readfile(filename=filename, filetype=filetype, wdir=wdir)
                if not _res is None:                        
                        self.data = {'ra':_res['ra'], 'dec':_res['dec'],'n':_res['n'],
                                     'fovra':_res['fovra'], 'fovdec':_res['fovdec']}
                
        def readfile(self, filename=None, filetype=None, wdir=None):
                """read pointings from file

                Parameters
                ----------               
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file, options: `txt` (plain text), `npz` (pickle file)
                wdir :             `string`
                  working directory for file

                Returns
                --------        
                data :          `dictionary`
                  Tilings: `ra`, `dec`, `fovra`, `fovdec`, `n`
                """
 
                if filename is None: filename = self.conf['filename']   
                if filetype is None: filetype = self.conf['filetype']                
                if not wdir is None: self.conf['wdir'] = wdir

                if filename is None:
                        self.logger.info ('### Warning: filename not defined')
                        return

                _data = {}                
                if filetype is None:     return
                
                elif filetype == 'npz':
                        cachefile = '%s/%s.npz' % (self.conf['wdir'], filename)
                        if os.path.exists(cachefile):
                                for _k in ['ra','dec','fovra','fovdec','n']:
                                        try:
                                                _data[_k] = np.load(cachefile)[_k]
                                        except:
                                                self.logger.info ('### Warning: missing keyword %s'%_k)
                                                return                               
                        else:
                                self.logger.info ('### Warning: %s not found'%cachefile)
                                return                        
                
                elif filetype == 'txt':
                        cachefile = '%s/%s.txt' % (self.conf['wdir'], filename)
                        if os.path.exists(cachefile):
                                _datat, _keys = {}, {}
                                _skipl = 1
                                for ll in open(cachefile).readlines():
                                        if ll[0] == '#':
                                                # read header
                                                for _nk,_k in enumerate(ll.replace('#','').split()):
                                                        _keys[_nk] = _k
                                                        _datat[_k] = []
                                        else:
                                                # read data
                                                for _nk,_k in enumerate(ll.split()):
                                                        if _nk in _keys:
                                                                try:
                                                                        _datat[_keys[_nk]].append(float(_k))
                                                                except:
                                                                        self.logger.info ('### Warning: skip line %i'%_skipl)
                                                                        _skipl += 1
                                for _k in ['ra','dec','fovra','fovdec','n']:
                                        if _k in _datat.keys():
                                                _data[_k] = _datat[_k]
                                        else:
                                                self.logger.info ('### Warning: %s not found'%_k)
                                                
                        else:
                                self.logger.info ('### Warning: %s not found'%cachefile)
                                return                        
                else:
                        self.logger.info ('### Warning: filetype %s unknown'%filetype)
                        return                

                return _data        
                        
        def astrotab(self):
                """output astropy table (https://docs.astropy.org/en/stable/table/)                

                Returns
                --------        
                table :          `astropy.tab`
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])
                >>> a.data
                {'n': array([    0,     1,     2, ..., 27788, 27789, 27790]), 'ra': array([  1.06418,   2.12836,   3.19253, ..., 286.53708, 315.19079, 343.8445 ]), 'dec': array([-20., -20., -20., ...,  88.,  88.,  88.]), 'fovra': array([1., 1., 1., ..., 1., 1., 1.]), 'fovdec': array([1., 1., 1., ..., 1., 1., 1.])}
                >>> a.astrotab()
                <Table length=27791>
                n       ra      dec    fovra   fovdec
                int64  float64  float64 float64 float64
                ----- --------- ------- ------- -------
                0   1.06418   -20.0     1.0     1.0
                1   2.12836   -20.0     1.0     1.0
                2   3.19253   -20.0     1.0     1.0
                3   4.25671   -20.0     1.0     1.0
                4   5.32089   -20.0     1.0     1.0
                5   6.38507   -20.0     1.0     1.0
                6   7.44924   -20.0     1.0     1.0
                7   8.51342   -20.0     1.0     1.0
                8    9.5776   -20.0     1.0     1.0
                9  10.64178   -20.0     1.0     1.0
                10  11.70596   -20.0     1.0     1.0
                11  12.77013   -20.0     1.0     1.0
                12  13.83431   -20.0     1.0     1.0
                13  14.89849   -20.0     1.0     1.0
                ...       ...     ...     ...     ...
                27776 305.71716    87.0     1.0     1.0
                27777 324.82448    87.0     1.0     1.0
                27778 343.93181    87.0     1.0     1.0
                27779  28.65371    88.0     1.0     1.0
                27780  57.30742    88.0     1.0     1.0
                27781  85.96113    88.0     1.0     1.0
                27782 114.61483    88.0     1.0     1.0
                27783 143.26854    88.0     1.0     1.0
                27784 171.92225    88.0     1.0     1.0
                27785 200.57596    88.0     1.0     1.0
                27786 229.22967    88.0     1.0     1.0
                27787 257.88338    88.0     1.0     1.0
                27788 286.53708    88.0     1.0     1.0
                27789 315.19079    88.0     1.0     1.0
                27790  343.8445    88.0     1.0     1.0                
                """

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                
                return Table([self.data['n'],
                              self.data['ra'],
                              self.data['dec'],
                              self.data['fovra'],
                              self.data['fovdec']],
                             names=('n', 'ra', 'dec', 'fovra', 'fovdec'))                 
               

        def remove_fileds_coo(self, ra, dec, fovra, fovdec,
                              nside=None, nest=None, skipfrac=None):
                """remove fields by input coordinate
                
                Parameters
                ------------
                ra :         `list`           
                  ra list to be skipped
                dec :        `list`           
                  dec list to be skipped
                fovra :      `list`           
                  fovra list to be skipped            
                fovdec :     `list`           
                  fovdec list to be skipped
                nside :      `int`
                  healpix nside parameter, must be a power of 2, less than 2**30
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead
                skipfrac :   `float` between 0 and 1
                  skipping fraction
                  tiling would be skipped when a fraction of it was already covered               
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> a.astrotab()
                <Table length=27791>
                n       ra      dec    fovra   fovdec
                int64  float64  float64 float64 float64
                ----- --------- ------- ------- -------
                0   1.06418   -20.0     1.0     1.0
                1   2.12836   -20.0     1.0     1.0
                2   3.19253   -20.0     1.0     1.0
                3   4.25671   -20.0     1.0     1.0
                4   5.32089   -20.0     1.0     1.0
                5   6.38507   -20.0     1.0     1.0
                6   7.44924   -20.0     1.0     1.0
                7   8.51342   -20.0     1.0     1.0
                8    9.5776   -20.0     1.0     1.0
                9  10.64178   -20.0     1.0     1.0
                10  11.70596   -20.0     1.0     1.0
                11  12.77013   -20.0     1.0     1.0
                12  13.83431   -20.0     1.0     1.0
                13  14.89849   -20.0     1.0     1.0
                ...       ...     ...     ...     ...
                27776 305.71716    87.0     1.0     1.0
                27777 324.82448    87.0     1.0     1.0
                27778 343.93181    87.0     1.0     1.0
                27779  28.65371    88.0     1.0     1.0
                27780  57.30742    88.0     1.0     1.0
                27781  85.96113    88.0     1.0     1.0
                27782 114.61483    88.0     1.0     1.0
                27783 143.26854    88.0     1.0     1.0
                27784 171.92225    88.0     1.0     1.0
                27785 200.57596    88.0     1.0     1.0
                27786 229.22967    88.0     1.0     1.0
                27787 257.88338    88.0     1.0     1.0
                27788 286.53708    88.0     1.0     1.0
                27789 315.19079    88.0     1.0     1.0
                27790  343.8445    88.0     1.0     1.0   
                >>> a.remove_fileds_coo([10],[10],[10],[10],skipfrac=0)
                <Table length=27670>
                n       ra      dec    fovra   fovdec
                int64  float64  float64 float64 float64
                ----- --------- ------- ------- -------
                0   1.06418   -20.0     1.0     1.0
                1   2.12836   -20.0     1.0     1.0
                2   3.19253   -20.0     1.0     1.0
                3   4.25671   -20.0     1.0     1.0
                4   5.32089   -20.0     1.0     1.0
                5   6.38507   -20.0     1.0     1.0
                6   7.44924   -20.0     1.0     1.0
                7   8.51342   -20.0     1.0     1.0
                8    9.5776   -20.0     1.0     1.0
                9  10.64178   -20.0     1.0     1.0
                10  11.70596   -20.0     1.0     1.0
                11  12.77013   -20.0     1.0     1.0
                12  13.83431   -20.0     1.0     1.0
                13  14.89849   -20.0     1.0     1.0
                ...       ...     ...     ...     ...
                27776 305.71716    87.0     1.0     1.0
                27777 324.82448    87.0     1.0     1.0
                27778 343.93181    87.0     1.0     1.0
                27779  28.65371    88.0     1.0     1.0
                27780  57.30742    88.0     1.0     1.0
                27781  85.96113    88.0     1.0     1.0
                27782 114.61483    88.0     1.0     1.0
                27783 143.26854    88.0     1.0     1.0
                27784 171.92225    88.0     1.0     1.0
                27785 200.57596    88.0     1.0     1.0
                27786 229.22967    88.0     1.0     1.0
                27787 257.88338    88.0     1.0     1.0
                27788 286.53708    88.0     1.0     1.0
                27789 315.19079    88.0     1.0     1.0
                27790  343.8445    88.0     1.0     1.0
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                
                if not skipfrac is None: self.conf['skipfrac'] = skipfrac
                if not nside is None: self.conf['nside'] = nside
                if not nest is None: self.conf['nest'] = nest 
                areasingle =  hp.nside2pixarea(self.conf['nside'], degrees=True)

                _data = {'n':np.array([]), 'ra':np.array([]), 'dec':np.array([]),
                         'fovra':np.array([]), 'fovdec':np.array([])}
                idtiles = KBGetTilings.ipix_in_box(self.data['ra'], self.data['dec'],
                                self.data['fovra'],self.data['fovdec'],self.conf['nside'],
                                self.conf['nest'])                        
                idhpx = KBGetTilings.ipix_in_box(np.array(ra), np.array(dec), np.array(fovra),
                                np.array(fovdec), self.conf['nside'], self.conf['nest'])
                ele = self.common_ele(idtiles, idhpx)                        
                nn = 0
                for nii,ii in enumerate(idtiles):
                        _ele = ele[nn: nn+len(ii)]
                        _areasi = np.count_nonzero(_ele)*areasingle
                        nn += len(ii)
                        _frac = _areasi/self.data['fovra'][nii]/self.data['fovra'][nii]                                
                        if _frac <= self.conf['skipfrac']:
                                _data['n']=np.append(_data['n'], self.data['n'][nii])
                                _data['ra']=np.append(_data['ra'], self.data['ra'][nii])
                                _data['dec']=np.append(_data['dec'], self.data['dec'][nii])
                                _data['fovra']=np.append(_data['fovra'], self.data['fovra'][nii])
                                _data['fovdec']=np.append(_data['fovdec'], self.data['fovdec'][nii])                
                self.data = _data
        
        def remove_fields_file(self, filename=None, filetype=None, wdir=None,
                               nside=None, skipfrac=None, nest=None):
                """remove fields by input file
                
                Parameters
                ------------         
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file, options: `txt` (plain text), `npz` (pickle file)
                wdir :             `string`
                  working directory for file       
                nside :      `int`
                  healpix nside parameter, must be a power of 2, less than 2**30
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead
                skipfrac :   `float` between 0 and 1
                  skipping fraction
                  tiling would be skipped when a fraction of it was already covered               
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> a.remove_fields_file()
                """
                
                if filename is None: filename = self.conf['filename']  
                if filetype is None: filetype = self.conf['filetype']                
                if wdir is None: wdir = self.conf['wdir']
                if nside is None: nside = self.conf['nside']  
                if skipfrac is None: skipfrac = self.conf['skipfrac']
                if nest is None: nest = self.conf['nest']
                
                _res = self.readfile(filename=filename, filetype=filetype, wdir=wdir)
                if _res is None: return

                # tilings should be skipped
                ra,dec,fovra,fovdec = _res['ra'], _res['dec'], _res['fovra'], _res['fovdec']

                # remove coolist
                self.remove_fileds_coo(ra, dec, fovra, fovdec,
                        nest=nest, nside=nside, skipfrac=skipfrac)
        
        def save(self, filetype=None, filename=None, wdir=None):
                """save tilings to a file
                
                Parameters
                ------------         
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file, options: `txt` (plain text), `npz` (pickle file)
                wdir :             `string`
                  working directory for file                        
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> a.save()
                """
 
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return

                if not filetype is None: self.conf['filetype'] = filetype
                if not filename is None: self.conf['filename'] = filename
                if not wdir is None: self.conf['wdir'] = wdir
                
                if self.conf['filetype'] is None:
                        return
                
                elif self.conf['filetype'] == 'npz':
                        cachefile = '%s/%s.npz' % (self.conf['wdir'], self.conf['filename'])
                        if os.path.exists(cachefile): os.remove(cachefile)
                        
                        # store to cachefile                        
                        np.savez(cachefile,ra=self.data['ra'],dec=self.data['dec'],
                                 fovra=self.data['fovra'],fovdec=self.data['fovdec'],
                                 n=self.data['n'])
                                
                elif self.conf['filetype'] == 'txt':
                        cachefile = '%s/%s.txt' % (self.conf['wdir'], self.conf['filename'])
                        if os.path.exists(cachefile): os.remove(cachefile)

                        # store to cachefile
                        ww = open(cachefile,'w')                       
                        ww.write('# ra dec fovra fovdec n \n')
                        for _ra,_dec,_fra,_fdec,_n in zip(self.data['ra'], self.data['dec'],
                                        self.data['fovra'],self.data['fovdec'],self.data['n']):
                                ww.write('%i %.5f %.5f %.2f %.2f \n'%(_n,_ra,_dec,_fra,_fdec))
                        ww.close()
                else:
                        self.logger.info ('### Warning: filetype %s unknown'%self.conf['filetype'])

        def calc_prob_loc(self, triggerobj, nest=None):                
                """calculate each tilings the probability that can cover targeting sources
                
                Parameters
                ------------         
                triggerobj :       `class`
                  KBParseTriggers object
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead                   
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
                >>> b = KBParseTriggers()   
                >>> b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
                >>> c = a.calc_prob_loc(b)
                >>> c
                <Table length=27791>
                n            prob
                int64        float64
                ----- ----------------------
                0  7.762505447792172e-05
                1 0.00012162028911602872
                2 0.00018488661247951556
                3 0.00027625315863596554
                4 0.00038066692835706515
                5  0.0004793012342172087
                6  0.0005629543081437528
                7  0.0006256609931332212
                8  0.0006782250883419683
                9  0.0005599687530534132
                10 0.00048668849730236496
                11  0.0004166983917616087
                12  0.0003430322368362478
                13  0.0002687973618959957
                ...                    ...
                27776 1.6135793843685656e-05
                27777   1.53210855578347e-05
                27778 1.5388783387778718e-05
                27779   9.30311429730731e-06
                27780  8.821616723970638e-06
                27781  7.765044427002605e-06
                27782 1.0692238193097823e-05
                27783 1.2511435412197865e-05
                27784 1.3244990370456746e-05
                27785 1.3703530269145885e-05
                27786 1.3517874076051426e-05
                27787 1.2687796993406217e-05
                27788 1.6181818795770013e-05
                27789 1.5537190545016493e-05
                27790 1.5394664113960722e-05
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if nest is None: nest = self.conf['nest']
                
                from kobe.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                elif is_seq(triggerobj.data['hpmap']):
                        hpx = triggerobj.data['hpmap']
                else: return
                
                probs, ns = [], []
                nside = hp.get_nside(hpx)
                for _ra,_dec,_fovw,_fovh,_n in zip(self.data['ra'], self.data['dec'],
                                self.data['fovra'], self.data['fovdec'], self.data['n']):
                        ipix_poly=(self.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,nest))
                        _probs = hpx[ipix_poly].sum()
                        probs.append(_probs)
                        ns.append(_n)
                return Table([ns, probs], names=('n', 'prob'))

        def set_contour(self, triggerobj, cls=None, nest=None, frac=0):
                '''
                '''
                if cls is None:      cls  = triggerobj.conf['cls']
                from kobe.cookbook import is_seq
                if is_seq(cls):  cls = [cls[0]]
                _data = self.cut_contours(triggerobj, cls=cls, nest=nest, frac=frac)
                self.data = _data[cls[0]]
                
        def cut_contours(self, triggerobj, cls=None, nest=None, frac=0):                
                """remove tilings the probability that can cover targeting sources
                
                Parameters
                ------------         
                triggerobj :       `class`
                  KBParseTriggers object
                cls :         `list`
                  list of confidence level, default: [.5, .9]
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead                   
                frac :   `float` between 0 and 1
                  fraction, tiling would be remained when its fraction that covered by a CL region is larger than `float`

                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
                >>> b = KBParseTriggers()   
                >>> b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
                >>> c = a.cut_contours(b, cls=[.9])
                >>> c
                {0.9: {'n': array([0.0000e+00, 1.0000e+00, 2.0000e+00, ..., 2.7781e+04, 2.7782e+04, 2.7783e+04]), 'ra': array([  1.06418,   2.12836,   3.19253, ...,  85.96113, 114.61483, 143.26854]), 'dec': array([-20., -20., -20., ...,  88.,  88.,  88.]), 'fovra': array([1., 1., 1., ..., 1., 1., 1.]), 'fovdec': array([1., 1., 1., ..., 1., 1., 1.])}}
                """
                
                if cls is None:      cls  = triggerobj.conf['cls']
                if nest is None:     nest = self.conf['nest'] 
                
                idlist = triggerobj.calc_contours(cls=cls)
                from kobe.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                elif is_seq(triggerobj.data['hpmap']):
                        hpx = triggerobj.data['hpmap']
                else: return
                
                nside = hp.get_nside(hpx)
                areasingle =  hp.nside2pixarea(nside, degrees=True)
                
                _data =  {}                                                        
                for cc in idlist:
                        _data[cc] =  {'n':np.array([]), 'ra':np.array([]), 'dec':np.array([]),
                                      'fovra':np.array([]), 'fovdec':np.array([])}
                        idhpx = idlist[cc]                        
                        idtiles = KBGetTilings.ipix_in_box(self.data['ra'], self.data['dec'],
                                        self.data['fovra'],self.data['fovdec'],nside,nest)
                        ele = self.common_ele(idtiles, idhpx)                        
                        nn = 0
                        for nii,ii in enumerate(idtiles):
                                _ele = ele[nn: nn+len(ii)]
                                _areasi = np.count_nonzero(_ele)*areasingle
                                nn += len(ii)
                                _frac = _areasi/self.data['fovra'][nii]/self.data['fovra'][nii]                                
                                if _frac  > frac:
                                        _data[cc]['n']=np.append(_data[cc]['n'], self.data['n'][nii])
                                        _data[cc]['ra']=np.append(_data[cc]['ra'], self.data['ra'][nii])
                                        _data[cc]['dec']=np.append(_data[cc]['dec'], self.data['dec'][nii])
                                        _data[cc]['fovra']=np.append(_data[cc]['fovra'], self.data['fovra'][nii])
                                        _data[cc]['fovdec']=np.append(_data[cc]['fovdec'], self.data['fovdec'][nii])
                return _data
                
        
        def calc_prob_dis(self, triggerobj, nest=None, limdist=400):               
                """calculate each tilings the probability that can reach targeting sources
                
                Parameters
                ------------         
                triggerobj :       `class`
                  KBParseTriggers object
                nest :             `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead                   
                limdist :          `float`
                  limiting distance of telescope on targeting sources, e.g. kilonovae, unit in Mpc
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
                >>> b = KBParseTriggers()   
                >>> b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
                >>> c = a.calc_prob_dis(b)
                >>> c
                <Table length=27791>
                n            prob
                int64        float64
                ----- ----------------------
                0 2.3419805886408296e-06
                1  2.500055445348643e-06
                2  2.544504009265499e-06
                3  2.719092839650723e-06
                4 2.8069783715468555e-06
                5 2.9320193609482796e-06
                6  2.987128541519283e-06
                7  3.063352038540268e-06
                8  3.097282271848335e-06
                9 3.0991167082713078e-06
                10 3.0892570423443215e-06
                11 3.0842840851790526e-06
                12 3.0677284821807573e-06
                13  3.012545385329722e-06
                ...                    ...
                27776  6.781067862926155e-06
                27777  6.553070919138916e-06
                27778  6.553070919138916e-06
                27779 5.4498666141403325e-06
                27780 5.4498666141403325e-06
                27781 4.6269040182437215e-06
                27782  5.144928039011793e-06
                27783  5.144928039011793e-06
                27784  4.688104407494185e-06
                27785  5.341106408103707e-06
                27786  5.628845287057135e-06
                27787  5.388979196318181e-06
                27788  6.781067862926155e-06
                27789  6.181369738911302e-06
                27790  6.553070919138916e-06
                """
                                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if nest is None: nest = self.conf['nest']
                
                from kobe.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                else:
                        return

                from scipy.stats import norm                             
                r = np.linspace(0, limdist)
                nside = hp.get_nside(hpx)
                pixarea = hp.nside2pixarea(nside, degrees=True)
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=nest)                
                dmu, dsigma, dnorm = hpd1[ipix], hpd2[ipix], hpd3[ipix]
                probl = [dnorm * norm(dmu, dsigma).pdf(rr)/pixarea for rr in r]
                probs = [max(ii) for ii in list(map(list, zip(*probl)))]               
                return Table([self.data['n'], probs], names=('n', 'prob'))  

        def calc_prob_mass(self, galaxyobj, nside=None, nest=None):               
                """calculate probability of local luminous mass (galaxies) for each tilings
                
                Parameters
                ------------         
                galaxyobj :       `class`
                  KBGetGalaxies object
                nside :      `int`
                  healpix nside parameter, must be a power of 2, less than 2**30
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead                   
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> from kobe.pipeline.KBGetGalaxies import KBGetGalaxies
                >>> b = KBGetGalaxies()   
                >>> b.generate(limdec=[-20,90],catalog='GLADE',limdist=[0,40])
                >>> c = a.calc_prob_mass(b)
                >>> c
                <Table length=27791>
                n            prob
                int64        float64
                ----- ----------------------
                0                    0.0
                1                    0.0
                2                    0.0
                3                    0.0
                4                    0.0
                5                    0.0
                6                    0.0
                7                    0.0
                8                    0.0
                9                    0.0
                10                    0.0
                11                    0.0
                12                    0.0
                13                    0.0
                ...                    ...
                27776                    0.0
                27777                    0.0
                27778                    0.0
                27779                    0.0
                27780                    0.0
                27781                    0.0
                27782                    0.0
                27783                    0.0
                27784                    0.0
                27785                    0.0
                27786                    0.0
                27787                    0.0
                27788                    0.0
                27789 0.00024173661941324885
                27790                    0.0
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if nside is None: nside = self.conf['nside']
                if nest is None: nest = self.conf['nest']
                
                hpx = galaxyobj.hpmap(nside=nside)
                probs, ns = [], []               
                for _ra,_dec,_fovw,_fovh,_n in zip(self.data['ra'], self.data['dec'],
                                self.data['fovra'], self.data['fovdec'], self.data['n']):
                        ipix_poly=(self.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,nest))
                        _probs = hpx[ipix_poly].sum()
                        probs.append(_probs)
                        ns.append(_n)
                return Table([ns, probs], names=('n', 'prob'))                

        def altaz(self, obstime=None, lat=None, lon=None, alt=None):
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

                returns
                ----------
                altaz :         `list`
                  altaz list

                obstime :      `astropy.time`
                  observing time, 

                observaroty :  `astropy.coordinate`
                  observaroty location
                """
                
                # observiting time
                if obstime is None: obstime = self.conf['obstime']                
                obstime = self.obstime(obstime)

                # observatory
                if not lat is None: self.conf['lat'] = lat
                if not lon is None: self.conf['lon'] = lon
                if not alt is None: self.conf['alt'] = alt
                observatory = astropy.coordinates.EarthLocation(lat=self.conf['lat']*u.deg,
                                lon=self.conf['lon']*u.deg, height=self.conf['alt']*u.m)
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg,
                                                      dec=self.data['dec']*u.deg)               

                # Alt/az reference frame at observatory, now
                frame = astropy.coordinates.AltAz(obstime=obstime, location=observatory)

                # Transform grid to alt/az coordinates at observatory, now
                altaz = radecs.transform_to(frame)
                return altaz, obstime, observatory
                
        def calc_airmass(self, obstime=None, lat=None, lon=None, alt=None):

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
                <Table length=27791>
                n        airmass
                int64      float64
                ----- ------------------
                0  2.082351877381677
                1 2.0366523500208884
                2 1.9937839280277412
                3 1.9535296936928823
                4  1.915696477565739
                5 1.8801104138666542
                6  1.846615265170311
                7  1.815069130505551
                8 1.7853445670112462
                9 1.7573257899791146
                10 1.7309076785530657
                11  1.705994857597817
                12 1.6824997774744856
                13  1.660343245357579
                ...                ...
                27777  2.361500041443534
                27778 2.2815119508872708
                27779 2.2437920691844204
                27780 2.2419237250450905
                27781 2.2787374614856644
                27782 2.3484040877350263
                27783  2.437960365979658
                27784  2.526884196120685
                27785  2.590320402338957
                27786  2.607713411858955
                27787  2.573007412223912
                27788 2.4981055056098107
                27789  2.406259305188414
                27790 2.3215526907548627 
                """
                  
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                return Table([self.data['n'], altaz.secz], names=('n', 'airmass'))  
        
        def calc_sun(self, obstime=None, lat=None, lon=None, alt=None):

                """calculate sun height for tilings

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
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the sun, now?
                sun_altaz = astropy.coordinates.get_sun(obstime).transform_to(altaz)
                return sun_altaz.alt

        def calc_solar(self, sobj, obstime=None, lat=None, lon=None, alt=None):
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
                <Table length=27791>
                n          dist
                deg
                int64      float64
                ----- ------------------
                0 111.48230429252158
                1 110.52744709508559
                2 109.57391448903458
                3  108.6217290640065
                4 107.67094152725195
                5 106.72159490501019
                6 105.77374246664934
                7 104.82741217070645
                8 103.88266016978994
                9  102.9395352333374
                10 101.99808772943797
                11 101.05837851768933
                12 100.12044371734643
                13  99.18434776772013
                ...                ...
                27777  69.84262000949794
                27778  69.11775063255271
                27779  67.07074034783025
                27780  66.17721327085206
                27781   65.6006729976947
                27782  65.48707626852809
                27783  65.86541699110153
                27784  66.63943402915353
                27785  67.61480362745114
                27786  68.55123667202297
                27787  69.22232875770946
                27788  69.46812448873096
                27789   69.2304760529254
                27790  68.56560022596703
                >>> b=a.calc_solar('sun')
                >>> b
                <Table length=27791>
                n          dist
                deg
                int64      float64
                ----- ------------------
                0  122.6472102176067
                1 121.67138532219919
                2  120.6964468468499
                3 119.72240300080793
                4 118.74929022624694
                5  117.7771365754901
                6 116.80598003086068
                7 115.83583204441057
                8 114.86673232977677
                9 113.89871238597732
                10  112.9318046668423
                11 111.96605169115686
                12 111.00146979487849
                13 110.03810362674838
                ...                ...
                27777  71.86662386557985
                27778   71.3055620381576
                27779  69.19218545294606
                27780  68.22709702079628
                27781   67.4870335327892
                27782  67.15772949901837
                27783  67.32279011585291
                27784  67.94023255598934
                27785   68.8543514227436
                27786  69.83825250213873
                27787  70.65206222160572
                27788  71.10037152954308
                27789  71.07656783320597
                27790  70.58629443695105
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if not sobj in ['earth', 'sun', 'earth-moon-barycenter', 'moon',
                                'mercury','venus','mars','jupiter', 'saturn',
                                'uranus','neptune']:
                        self.logger.info ('### Error: sobj wrong')
                        return
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg,
                                                      dec=self.data['dec']*u.deg)
                
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the solar object, now
                moonloc = astropy.coordinates.get_body(sobj, obstime, observatory).transform_to(altaz)
                return Table([self.data['n'], radecs.separation(moonloc)], names=('n', 'dist'))
        
        def divide_OB(self, nobw=None, nobh=None):
                """divide each tiling to a serious of sub-tilings

                Parameters
                ----------  
                nobw :         `int`
                  number of pointings in OB, in ra direction
                nobh :         `int`
                  number of pointings in OB, in dec direction

                Examples
                -------- 
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90]) 
                >>> a.divide_OB(3,3)
                >>> a.astrotab()
                <Table length=250119>
                n       ra      dec    fovra   fovdec
                int64  float64  float64 float64 float64
                ----- --------- ------- ------- -------
                1  -0.00696   -21.0 0.33333 0.33333
                1   1.06418   -21.0 0.33333 0.33333
                1   2.13532   -21.0 0.33333 0.33333
                1       0.0   -20.0 0.33333 0.33333
                1   1.06418   -20.0 0.33333 0.33333
                1   2.12836   -20.0 0.33333 0.33333
                1   0.00656   -19.0 0.33333 0.33333
                1   1.06418   -19.0 0.33333 0.33333
                1    2.1218   -19.0 0.33333 0.33333
                2   1.05722   -21.0 0.33333 0.33333
                2   2.12836   -21.0 0.33333 0.33333
                2    3.1995   -21.0 0.33333 0.33333
                2   1.06418   -20.0 0.33333 0.33333
                2   2.12836   -20.0 0.33333 0.33333
                ...       ...     ...     ...     ...
                27790 286.53708    88.0 0.33333 0.33333
                27790 315.19079    88.0 0.33333 0.33333
                27790  343.8445    88.0 0.33333 0.33333
                27790  257.8921    89.0 0.33333 0.33333
                27790 315.19079    89.0 0.33333 0.33333
                27790 372.48948    89.0 0.33333 0.33333
                27791 324.73718    87.0 0.33333 0.33333
                27791  343.8445    87.0 0.33333 0.33333
                27791 362.95182    87.0 0.33333 0.33333
                27791 315.19079    88.0 0.33333 0.33333
                27791  343.8445    88.0 0.33333 0.33333
                27791 372.49821    88.0 0.33333 0.33333
                27791 286.54581    89.0 0.33333 0.33333
                27791  343.8445    89.0 0.33333 0.33333
                27791 401.14319    89.0 0.33333 0.33333
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                
                if not nobw is None: self.conf['obra'] = nobw
                if not nobh is None: self.conf['obdec'] = nobh                        

                ralist, declist, fovralist, fovdeclist, oblist = [], [], [], [], []
                _nn = 0
                for _ra,_dec,_fovw,_fovh in zip(self.data['ra'], self.data['dec'],
                                                self.data['fovra'], self.data['fovdec']):
                        _ra0,_dec0,_fovw0,_fovh0 = self.divide_OB_one(_ra, _dec, _fovw, _fovh,
                                                        self.conf['obra'], self.conf['obdec'])
                        _nn+=1
                        for _ra00,_dec00,_fovw00,_fovh00 in zip(_ra0,_dec0,_fovw0,_fovh0):
                                ralist.append(float('%.5f'%_ra00))
                                declist.append(float('%.5f'%_dec00))
                                fovralist.append(float('%.5f'%_fovw00))
                                fovdeclist.append(float('%.5f'%_fovh00))
                                oblist.append(_nn)
                        
                self.data = {'ra':     np.array(ralist),
                             'dec':    np.array(declist),
                             'n':      np.array(oblist),
                             'fovra':  np.array(fovralist),
                             'fovdec': np.array(fovdeclist)}                
        @staticmethod
        def common_ele(list1, list2, nullv=-99):
                """ find common elements of list1 inside list2

                Parameters
                ----------   
                list1 :   sequence
                  an index list, to be checked
                list2 :   sequence
                  an index list

                returns
                ----------   
                ele :     sequence
                  common elements
                """               
                list1 = KBGetTilings.flatten(list1, nullv=nullv)
                list2 = KBGetTilings.flatten(list2, nullv=nullv)
                if not list1 is None and not list2 is None:                        
                        _mask = np.in1d(list1, list2)                                                
                        return _mask

        @staticmethod
        def flatten(_list, nullv=-99):
                from kobe.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(_list):
                        # check length of each elements
                        llist = np.unique([len(ll) for ll in _list])
                        lmax = max([len(ll) for ll in _list])                                
                        if len(llist) > 1:
                                # arrange each elements to same lenght
                                _listm = []
                                for ll in _list:
                                        lla = [nullv]*(lmax - len(ll))
                                        ll = np.append(ll, lla)
                                        _listm.append(ll)
                                _list = _listm
                        else:
                                pass
                        # flatten to 1d
                        _listf = [item for sublist in _list for item in sublist]
                elif is_seq(_list):
                        _listf = _list
                else:
                        return
                return _listf
        
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

        @staticmethod
        def divide_OB_one(rac, decc, fovw, fovh, nobw, nobh):
                """divide one pointing to a list of sub-pointings

                Parameters
                ----------   
                rac :         `float`
                  ra center of a tiling
                decc :        `float`           
                  dec center of a tiling
                fovw :        `float`           
                  fov in ra direction of a tiling          
                fovh :        `float`           
                  fov in dec direction of a tiling    
                nobw :        `int`
                  number of pointings in OB, in ra direction
                nobh :        `int`
                  number of pointings in OB, in dec direction
                
                returns
                ----------   
                ra :         `list`           
                  tiling ra list
                dec :        `list`           
                  tiling dec list
                fovra :      `list`           
                  tiling field of view list in ra direction              
                fovdec :     `list`           
                  tiling field of view list in dec direction
                """
                
                _ndec = np.arange(nobh)-(nobh-1)/2               
                _decdiff = fovh
                _decspread=decc+_ndec*_decdiff

                ralist,declist, fovralist,fovdeclist = [],[],[],[]
                for _dec in _decspread:
                        npoint = 360*np.cos(_dec*np.pi/180)/fovw
                        _radiff = 360/npoint
                        _nra = np.arange(nobw)-(nobw-1)/2
                        _raspread=rac+_nra*_radiff
                
                        for _ra in _raspread:      
                                ralist.append(_ra)
                                declist.append(_dec)
                                fovralist.append(fovw/nobw)
                                fovdeclist.append(fovh/nobh)
                return  ralist,declist, fovralist,fovdeclist 
        
        @staticmethod
        def ipix_in_box(ra,dec,width,height,nside,nest):
                """finding the healpix indices of a given box

                Parameters
                ----------   
                ra :         `float` or `numpy array`
                  center ra of a box
                dec :        `float` or `numpy array`
                  center dec of a box
                width :      `float` or `numpy array`
                  width ra of a box        
                height :     `float` or `numpy array`
                  height of a box
                nside :      `int`
                  healpix nside parameter, must be a power of 2, less than 2**30
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead

                returns
                ----------   
                indexlist :     `list`           
                  list of indices
                """
                                                
                v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = \
                        KBGetTilings.vertices(ra, dec, width, height)
                ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                             [v1_dec, v2_dec, v4_dec, v3_dec])                
                theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
                phi = np.deg2rad(ra_vertices)
                xyz = hp.ang2vec(theta, phi)               
                from kobe.cookbook import is_seq                            
                if is_seq(ra) and is_seq(dec) and is_seq(width) and is_seq(height):
                        ipix_fov_box = []
                        for _xyz in xyz:
                                ipix_fov_box.append(hp.query_polygon(nside, _xyz, nest=nest))
                else:
                        ipix_fov_box = hp.query_polygon(nside, xyz, nest=nest)
                return ipix_fov_box
        
        @staticmethod
        def vertices(ra,dec,fovw,fovh):
                """finding the vertices of a FoV given the central location (ra[deg], dec[deg])
                and the FoV size (fovw [deg], fovh [deg]).

                Parameters
                ----------   
                ra :         `float`  or `numpy array`      
                  center ra of a box
                dec :        `float`  or `numpy array`       
                  center dec of a box
                fovw :       `float`  or `numpy array`         
                  width ra of a box        
                fovh :       `float`  or `numpy array`       
                  height of a box               

                returns
                ----------   
                vertices :     `list`                  
                """
                
                fovw,fovh = fovw/2.,fovh/2.
                vert_ra,vert_dec=[],[]
                ra_rad,dec_rad,fovw_rad,fovh_rad = np.deg2rad(ra), np.deg2rad(dec),\
                        np.deg2rad(fovw), np.deg2rad(fovh)
                for i,j in zip([-fovw_rad, fovw_rad, fovw_rad, -fovw_rad],\
                               [fovh_rad, fovh_rad, -fovh_rad, -fovh_rad]):
                        arg = -i/(np.cos(dec_rad)-j*np.sin(dec_rad))
                        v_ra = np.rad2deg(ra_rad+np.arctan(arg))       
                        v_dec = np.rad2deg(np.arcsin((np.sin(dec_rad)+\
                                j*np.cos(dec_rad))/(1+i**2+j**2)**0.5))
                        vert_ra.append(v_ra)
                        vert_dec.append(v_dec)
                return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], \
                        vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]
