#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/tilings.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import numpy as np
import healpy as hp
import random
from astropy.table import Table
from astroquery.vizier import Vizier
import logging
from kobe import vtilings, vgalaxies, circulate

__all__ = ['pointings', 'tilings', 'pointings']

class tilings(vtilings, circulate):
        """tilings: Generate, read and process tilings

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
        skipfrac:        `float`  
           one vertice is skipped when how much fraction of it is covered        
        wdir    : `str`     
           define a working directory, that used to read and build tiling list        
        clobber   : `bool`  
           if any products already available in self, clobber=True will redo it 

        See Also
        --------
        kobe.galaxies, kobe.pointings    

        Notes
        ----------   
        *tilings* inherients from **vtilings** and **circulate**,
        and inheriented by *pointings*.    
        """               
        defkwargs = {**vtilings.defkwargs, **circulate.defkwargs}
        '''
        default parameter settings
        '''
 
        data = None
        '''
        tiling data
        '''

        hpmapm = None
        '''
        healpix map for mass
        '''
        def __new__(cls, logger=None, **kwargs):
                """
                generate a new tiling object with input logger and parameter settings        
                """
                instance = super(tilings, cls).__new__(cls)
                instance.__init__(logger=logger, **kwargs)        
                return instance
        
        def __init__(self, logger=None, **kwargs):
                """
                initialize tiling object with input logger and parameter settings
        
                Parameters
                ----------
                logger :     `class`
                   logging object

                **kwargs :     
                   check kobe.vtilings.defkwargs, 
                   kobe.circulate.defkwargs, kobe.tilings.defkwargs
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set keys ----- #        
                self.defkwargs = self.setkeys(kwargs)
        
        def generatep(self, limra=[0,360.], limdec=[-89,89], fovra=1.,
                      fovdec=1., shiftra=0., shiftdec=0., returndata=False, **kwargs):
                """generate pointings by tiling sky

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
                clobber   : `bool`  
                if any products already available in self, clobber=True will redo it 

                Examples
                --------                
                >>> from kobe import tilings
                >>> a = tilings()
                >>> a.generate_pointings(limdec=[-20,90])
                >>> a.data
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
                ...       ...     ...     ...     ...
                27783 143.26854    88.0     1.0     1.0
                27784 171.92225    88.0     1.0     1.0
                27785 200.57596    88.0     1.0     1.0
                27786 229.22967    88.0     1.0     1.0
                27787 257.88338    88.0     1.0     1.0
                27788 286.53708    88.0     1.0     1.0
                27789 315.19079    88.0     1.0     1.0
                27790  343.8445    88.0     1.0     1.0
                """
                kwargs = self.setkeys(kwargs)
        
                if not self.data is None and not kwargs['clobber']:
                        self.logger.info ('Warning: tiling data already parsed')
                        return                
        
                # limit shift to proper range
                shiftra = shiftra % fovra
                shiftdec = shiftdec % fovdec
                
                # cut ra,dec
                ramin, ramax = min(limra)+shiftra, max(limra)+shiftra
                decmin, decmax = min(limdec)+shiftdec, max(limdec)+shiftdec
                ramin = max(ramin, 0)
                ramax = min(ramax, 360)
                decmin = max(decmin, -89)
                decmax = min(decmax, 89)
                
                # dec range
                decrange= np.arange(decmin,decmax,fovdec)

                # get network, i.e. ra,dec list
                ralist,declist=[],[]
                fovralist,fovdeclist=[],[]
                for _dec in decrange:
                        npoint = 360*np.cos(_dec*np.pi/180)/fovra
                        for nn in np.arange(0,npoint,1):  
                                _ra = 360./npoint*nn+shiftra
                                if _ra < ramax and _ra > ramin:
                                        ralist.append(float('%.5f'%_ra))
                                        declist.append(float('%.5f'%_dec))
                                        fovralist.append(float('%.5f'%fovra))
                                        fovdeclist.append(float('%.5f'%fovdec))
                                        
                _data = Table([np.arange(len(ralist)), np.array(ralist),
                               np.array(declist), np.array(fovralist),
                               np.array(fovralist),],
                              names=('n', 'ra', 'dec', 'fovra', 'fovdec'))
                if returndata: return _data
                else: self.data = _data

        def inputp(self, ra=None, dec=None, fovra=None, fovdec=None, hms=False):
                
                radecs = self.radecs(ra=ra, dec=dec, hms=hms)                
                if radecs is None: return                
                _data = {}
                for _i, _j in zip([fovra,fovdec],['fovra','fovdec']):
                        assert not _i is None
                        if self.is_seq(_i):
                                assert len(_i) == len(radecs[0])
                        else:
                                _i = [_i]*len(radecs[0])
                        _data[_j] = np.array(_i)
                        
                self.data = Table([np.arange(len(radecs[0])),
                                radecs[0], radecs[1], _data['fovra'], _data['fovdec']],
                                names=('n', 'ra', 'dec', 'fovra', 'fovdec'))                
                                                        
        def readp(self, filename, filetype='npz', split=' ', returndata=False,
                           keys=['n','ra','dec','fovra','fovdec'], **kwargs):
                """read pointings from file

                Parameters
                ----------               
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file
                  options: `txt` (plain text), `npz` (pickle file), `fits` (astropy.fits)
                wdir :             `string`
                  working directory for file
                clobber   : `bool`  
                  if any products already available in self, clobber=True will redo it
                keys :   sequence
                  keys that read from file and construct data 
                split :      `string`
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)
                """                
                kwargs = self.setkeys(kwargs)
                
                if not self.data is None and not kwargs['clobber']:
                        self.logger.info ('Warning: tiling data already parsed')
                        return
                
                _r = self.readlist(split=split, filename=filename,
                                   filetype=filetype, keys=keys, **kwargs)
                if _r is None: return
               
                _data = []
                for _key in keys:
                        _data.append(np.array(_r[_key]))
                if returndata:
                        return Table(_data, names=keys)
                else:
                        self.data = Table(_data, names=keys)
                                
        def removep_coo(self, ra, dec, fovra, fovdec, data=None,
                        skipfrac=0, nside=512, nest=False):
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
                >>> a.remove_pointings_coo([10],[10],[10],[10],skipfrac=0)
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
                if not data is None:
                        datain = data
                else:
                        datain = self.data
                        
                if datain is None:
                        self.logger.info ('Warning: tiling data not found')
                        return
                
                areasingle =  hp.nside2pixarea(nside, degrees=True)
                _data = {'n':np.array([]), 'ra':np.array([]), 'dec':np.array([]),
                         'fovra':np.array([]), 'fovdec':np.array([])}
                idtiles = self.ipix_in_box(datain['ra'], datain['dec'],
                        datain['fovra'], datain['fovdec'], nside, nest)                
                id2skip = self.ipix_in_box(np.array(ra), np.array(dec), np.array(fovra),
                        np.array(fovdec), nside, nest)
                ele = self.common_ele(idtiles, id2skip)
                lmax = max([len(ll) for ll in idtiles])
                
                nn = 0
                for ii in range(len(idtiles)):
                        _ele = ele[ii*lmax: (ii+1)*lmax-1]                        
                        _areasi = np.count_nonzero(_ele)*areasingle                        
                        _frac = _areasi/datain['fovra'][ii]/datain['fovra'][ii]
                        if _frac <= skipfrac:
                                _data['n']=np.append(_data['n'], datain['n'][ii])
                                _data['ra']=np.append(_data['ra'], datain['ra'][ii])
                                _data['dec']=np.append(_data['dec'], datain['dec'][ii])
                                _data['fovra']=np.append(_data['fovra'], datain['fovra'][ii])
                                _data['fovdec']=np.append(_data['fovdec'], datain['fovdec'][ii])
                if not data is None:
                        return Table(_data)
                else:
                        self.data = Table(_data)   
        
        def removep_file(self, split=' ', nest=False, skipfrac=0, data=None,
                         filename=None, filetype=None, nside=512, **kwargs):
                """remove fields by input file
                
                Parameters
                ------------         
                filename :         `string`
                  name of file
                  mush contain fovra and fovdec, so that define vertices
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
                split :      `string`
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)
                
                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> a.remove_pointings_file(filename='tmp.npz', filetype='npz')
                """                               
                _r = self.readlist(split=split, filename=filename,
                                filetype=filetype, keys=keys, **kwargs)         
                if _r is None: return
                
                # remove coolist
                return self.removep_coo(_r['ra'], _r['dec'], _r['fovra'], _r['fovdec'],
                                        data=data, nside=nside, skipfrac=skipfrac)
                
        def savep(self, filename, filetype='npz',split=' ',data=None,
                           keys=['n','ra','dec','fovra','fovdec'], **kwargs):
                """save tilings to a file
                
                Parameters
                ------------         
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file
                  options: `txt` (plain text), `npz` (pickle file), `fits` (astropy.fits)
                wdir :             `string`
                  working directory for file                        
                split :      `string`
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)

                Examples
                --------                
                >>> from kobe.pipeline.KBGetTilings import KBGetTilings
                >>> a = KBGetTilings()
                >>> a.generate(limdec=[-20,90])               
                >>> a.save_pointings('tmp.npz', filetype='npz')
                """                
                self.writelist(datain=data, filename=filename, filetype=filetype,
                               split=split, keys=keys, **kwargs)

        def preport(self, split=' '):       
                """make a report for data

                Parameters
                ----------   
                split     : `str`
                """        
                _report = self.dic2txt(split=split)
                if not _report in self.texts: self.texts += _report
        
        def group_OB(self, obra, obdec):
                """group sub-tilings to a joint tiling

                Parameters
                ----------  
                obra :         `float`           
                  tiling field of view in ra direction to group pointings
                obdec :        `float`           
                  tiling field of view in dec direction to group pointings
                """
                assert type(obra) in [float,int] and type(obdec) in [float,int]
                
                if self.data is None:
                        self.logger.info ('Warning: pointings not parsed')
                        return

                ral = np.arange(min(self.data['ra']), max(self.data['ra']), obra)
                decl = np.arange(min(self.data['dec']), max(self.data['dec']), obdec)
                X, Y = np.meshgrid(ral, decl)
                                
                _n=0
                for z in [i for i in zip(X.flat,Y.flat)]:                        
                        _id1 = (self.data['ra']>=z[0])                        
                        _id2 = (self.data['ra']< z[0]+obra)
                        _id3 = (self.data['dec']>=z[1])
                        _id4 = (self.data['dec']< z[1]+obdec)
                        
                        self.data['n'][_id1 & _id2 & _id3 & _id4] = _n                        
                        if len(self.data[_id1 & _id2 & _id3 & _id4])>0:_n+=1
                        
        def divide_OB(self, obra, obdec):
                """divide each tiling to a serious of sub-tilings

                Parameters
                ----------  
                obra :         `int`
                  number of pointings in OB, in ra direction
                obdec :         `int`
                  number of pointings in OB, in dec direction

                Examples
                -------- 
                >>> from kobe import tilings
                >>> a = tilings()
                >>> a.generate_pointings(limdec=[-20,90]) 
                >>> a.divide_OB(obra=3,obdec=3)
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
                assert type(obra) is int and type(obdec) is int
                
                if self.data is None:
                        self.logger.info ('Warning: pointings not parsed')
                        return
                
                ralist, declist, fovralist, fovdeclist, oblist = [], [], [], [], []                
                for _nn, _ra,_dec,_fovw,_fovh in zip(self.data['n'], self.data['ra'],
                        self.data['dec'], self.data['fovra'], self.data['fovdec']):
                        _ra0,_dec0,_fovw0,_fovh0 = self.divide_OB_one(_ra,_dec,_fovw,_fovh,obra,obdec)
                        for _ra00,_dec00,_fovw00,_fovh00 in zip(_ra0,_dec0,_fovw0,_fovh0):
                                ralist.append(float('%.5f'%_ra00))
                                declist.append(float('%.5f'%_dec00))
                                fovralist.append(float('%.5f'%_fovw00))
                                fovdeclist.append(float('%.5f'%_fovh00))
                                oblist.append(_nn)
                        
                self.data = Table({'n':      np.array(oblist),
                                   'ra':     np.array(ralist),
                                   'dec':    np.array(declist),                                   
                                   'fovra':  np.array(fovralist),
                                   'fovdec': np.array(fovdeclist)})
                
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
                _decdiff = fovh/nobh
                _decspread=decc+_ndec*_decdiff
        
                ralist,declist, fovralist,fovdeclist = [],[],[],[]
                for _dec in _decspread:
                        npoint = 360*np.cos(_dec*np.pi/180)/fovw
                        _radiff = 360/npoint/nobw
                        _nra = np.arange(nobw)-(nobw-1)/2
                        _raspread=rac+_nra*_radiff
                        for _ra in _raspread:      
                                ralist.append(_ra)
                                declist.append(_dec)
                                fovralist.append(fovw/nobw)
                                fovdeclist.append(fovh/nobh)                
                return  ralist,declist, fovralist,fovdeclist 

        def rankp(self, mode=1, nside=64, nest=False, threshold=1,sort=1):
                """rank pointings: rank from west to east

                Parameters
                ----------                   
                mode     : `int`
                  1. strict approach: follow above strategy strictly
                  2. adjacent approach: start from the westest point or 
                       highest probability point and arange the next pointing 
                       either adjacent or closest to the previous one
                """                
                assert 'data' in self.__dict__                     
                assert not self.data is None, 'Warning: tiling data not parsed'
                _data = self.data        
        
                # define start point                
                _idstart = np.argmin(_data['ra'])                
        
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
                
                self.data = _data[ipx]
                self.data['n'] = np.arange(len(self.data['n']), dtype=int)
 


class galaxies(vgalaxies, circulate):
        """galaxies: Generate (by querying astropy.Vizier), 
        read and process galaxies

        Parameters
        ----------
        limra :         `range`
           galaxies: ra range, default: [0, 360]
        limdec :        `range`
           galaxies: dec range, default: [-90, 90]
        limdist :       `range`
           galaxies: distance range, default: [0, 100]           
        limmag :        `range`
           galaxies: absolute magnitude range, default: [-18, -99]
        catalog :       `string`
           Vizier galaxy catalog
           options: `GWGC`, `GLADE`
           default: GLADE
        filter :        `string`
           magnitude filter. 
           `B` and `K` are available for `GLADE` 
           `B` available for `GWGV`
           default: B
        size :      `float`
           size of the querying galaxies, -1 for the full query. 
           default: -1                      
        wdir    : `str`     
           define a working directory, that used to read and build tiling list       
        clobber   : `bool`  
           if any products already available in self, clobber=True will redo it

        See Also
        --------
        kobe.tilings, kobe.telescope
        """                           
        defkwargs = {**vgalaxies.defkwargs, **circulate.defkwargs,}
        '''
        default parameter settings
        '''
 
        data = None
        '''
        galaxy data
        '''
        def __new__(cls, logger=None, **kwargs):
                """
                generate a new galaxy object with input logger and parameter settings        
                """
                instance = super(galaxies, cls).__new__(cls)
                instance.__init__(logger=logger, **kwargs)        
                return instance
        
        def __init__(self, logger=None, **kwargs):
                """
                initialize galaxy object with input logger and parameter settings
        
                Parameters
                ----------
                logger :     `class`
                   logging object

                **kwargs :     
                   check kobe.vgalaxies.defkwargs, 
                   kobe.circulate.defkwargs, kobe.galaxy.defkwargs
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger
                        
                # ----- set keys ----- #        
                self.defkwargs = self.setkeys(kwargs)
                
        def generatep(self, catalog='GLADE', filter='B', size=-1,
                      limra=[0,360.], limdec=[-89,89], limdist=[0,1000],
                      limmag=[-18,-99], returndata=None, **kwargs):
                """generate galaxies by querying Vizier
                
                Parameters
                ----------
                limra :         `range`
                  galaxies: ra range, default: [0, 360]
                limdec :        `range`           
                  galaxies: dec range, default: [-90, 90]
                limdist :       `range`           
                  galaxies: distance range, default: [0, 1000]
                limmag :        `range`           
                  galaxies: absolute magnitude range, default: [-18, -99]
                catalog :       `string`                  
                  Vizier galaxy catalog, options: `GWGC`, `GLADE`, default: GLADE
                filter :        `string`           
                  magnitude filter. 
                  `B` and `K` are available for `GLADE` 
                  while only `B` available for `GLADE`
                  default: B
                size :      `float`           
                  size of the querying galaxies, -1 for the full query. default: -1
                clobber   : `bool`  
                  if any products already available in self, clobber=True will redo it

                Examples
                --------                
                >>> from kobe import galaxies
                >>> a = galaxies()
                >>> a.generate_pointings(limdec=[-20,90], limdist=[0,40])
                INFO:kobe.KBtelescope:2492 galaxies selected from VII/281
                >>> a.data
                <Table length=2492>
                n                    name                      ra        dec       distance     mag
                int64                 str62                   float64    float64     float64    float64
                ----- -------------------------------------- ---------- --------- ------------- --------
                0 28655:NGC3034:NGC3034:09555243+6940469  148.96846 69.679703 4.70228466231 -19.2115
                1 42407:NGC4594:NGC4594:12395949-1137230 189.997894 -11.62307 3.65995814653 -19.2974
                2            --:---:---:12564369+2140575 194.182068 21.682659 3.87868462045 -18.6564
                3 41220:NGC4472:NGC4472:12294679+0800014 187.444992   8.00041 14.6955181559 -20.9739
                4 47404:NGC5194:NGC5194:13295269+4711429 202.469574 47.195259 2.49371678269 -18.3042
                5 10266:NGC1068:NGC1068:02424077-0000478   40.66988  -0.01329 18.4142043388 -21.9158
                6 42831:NGC4649:NGC4649:12434000+1133093 190.916702  11.55261 16.4529675514 -21.3112
                7 33550:NGC3521:NGC3521:11054859-0002092 166.452469   -0.0359 11.2915350655 -20.8638
                8 41361:NGC4486:NGC4486:12304942+1223279 187.705933   12.3911 18.5306296927 -21.4225
                9 29265:NGC3115:NGC3115:10051397-0743068 151.308243  -7.71858 6.51337285663   -19.04
                10              34695:NGC3627:NGC3627:---    170.062   12.9916 12.9133823516 -21.1252
                11   14508:IC0356:IC0356:04074690+6948447   61.94545 69.812431 16.4486586923 -20.9407
                12              69327:NGC7331:NGC7331:---    339.267  34.41562 18.1813685082 -21.6281
                13              27077:NGC2903:NGC2903:---    143.042  21.50141 10.9303004689 -21.0632
                ...                                    ...        ...       ...           ...      ...
                2477                  34426:NGC3607:---:---    169.227  18.05301 24.9589462382 -21.0361
                2478                 38905:UGC07207:---:---     183.08  37.01366  21.227046734 -18.3744
                2479                  49354:NGC5354:---:---    208.362  40.30397 38.9986946803 -20.7053
                2480                  60315:NGC6368:---:---    261.798  11.54362 32.8410432488 -20.7821
                2481               135878:PGC135878:---:---    344.386    -2.501 38.2733015989 -18.3245
                2482               166072:PGC166072:---:---    50.7864   62.7893 8.30849177211 -18.3676
                2483             2801052:PGC2801052:---:---    315.851  57.28721 13.5503403375 -18.7498
                2484                2802656:NGC5904:---:---    229.639   2.08277 2.39941051318 -19.5605
                2485              38897:NGC4173:NGC4173:---    183.089  29.20702 21.0457642607 -19.2158
                2486     --:---:SDSSJ022739.95-010913.2:---    36.9165  -1.15368 21.2227303453  -18.494
                2487     --:---:SDSSJ123616.69+260006.9:---     189.07  26.00194 25.2571985761 -18.3719
                2488     --:---:SDSSJ123626.88+255738.2:---    189.112  25.96063 20.7004974885 -18.9699
                2489     --:---:SDSSJ122250.38+155056.9:---     185.71  15.84916 24.8768244613  -18.159
                2490     --:---:SDSSJ124211.24+074016.0:---    190.547    7.6711 29.7521819633 -18.0476
                2491                     --:---:NGC5496:---    212.908  -1.15744 23.4249708242 -18.3384
                """
                
                # ----- set keys ----- #        
                kwargs = self.setkeys(kwargs)    

                if not self.data is None and not kwargs['clobber']:
                        self.logger.info ('Warning: galaxy data already parsed')
                        return
                
                # specify columns
                if catalog == 'GLADE':
                        if not filter in ['B', 'K']:
                                self.logger.info ('Error: wrong filters for GLADE')
                                return
                        catid = 'VII/281'
                        columns = ['RAJ2000', 'DEJ2000', '%sMAG'%filter, \
                                   'Dist', 'PGC', 'GWGC', 'HyperLEDA', '2MASS']
                elif catalog == 'GWGC':
                        if filter != 'B':
                                self.logger.info ('Error: wrong filters for GWGC')
                                return
                        catid = 'VII/267'
                        columns =  ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                                    'Dist', 'Name']
                else:
                        self.logger.info ('Error: wrong galaxy catalog')
                        return

                # download catalog with vizier                
                v = Vizier(columns=columns, column_filters={\
                        columns[0]:'%s..%s'%(limra[0],limra[1]),\
                        columns[1]:'%s..%s'%(limdec[0],limdec[1]),\
                        columns[2]:'%s..%s'%(limmag[0],limmag[1]),\
                        columns[3]:'%s..%s'%(limdist[0],limdist[1])
                })
                v.ROW_LIMIT = size
                catalogs = v.get_catalogs(catid)[0]
                self.logger.info ("%i galaxies selected from %s"%(len(catalogs),catid))

                # return infos    
                if catalog == 'GLADE':
                        _name = []
                        for ii in range(len(catalogs)):
                                _name.append('%s:%s:%s:%s'%(catalogs[columns[4]][ii], \
                                                            catalogs[columns[5]][ii], \
                                                            catalogs[columns[6]][ii], \
                                                            catalogs['_%s'%columns[7]][ii]))
                else:
                        _name = catalogs[columns[4]]
                        
                _data = Table({'n':       np.arange(len(catalogs[columns[0]])),
                               'name':    np.array(_name),
                               'ra':      np.array(catalogs[columns[0]]),
                               'dec':     np.array(catalogs[columns[1]]),
                               'mag':     np.array(catalogs[columns[2]]),
                               'dist':    np.array(catalogs[columns[3]])})                        
                if returndata:
                        return _data
                else:
                        self.data = _data

        def readp(self, filename, filetype='npz', split=' ', returndata=None,
                           keys=['n','name','ra','dec','mag','dist'], **kwargs):
                """read pointings from file

                Parameters
                ----------               
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file
                  options: `txt` (plain text), `npz` (pickle file), `fits` (astropy.fits)
                wdir :             `string`
                  working directory for file
                clobber   : `bool`  
                if any products already available in self, clobber=True will redo it 
                split :      `string`
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)
                """                
                kwargs = self.setkeys(kwargs)
                
                if not self.data is None and not kwargs['clobber']:
                        self.logger.info ('Warning: galaxy data already parsed')
                        return
                
                _r = self.readlist(split=split, filename=filename,
                                   filetype=filetype, keys=keys, **kwargs)
                if _r is None: return
               
                _data = []
                for _key in keys:
                        _data.append(np.array(_r[_key]))

                if returndata:
                        return Table(_data, names=keys)   
                else:
                        self.data = Table(_data, names=keys)
                        
        def removep_coo(self, ra, dec, fovra, fovdec, data=None,
                        nside=512, nest=False):
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

                Examples
                --------                
                >>> from kobe import galaxies
                >>> a = galaxies()
                >>> a.generate_pointings(limdec=[-20,90], limdist=[0,40])
                INFO:kobe.KBtelescope:2492 galaxies selected from VII/281
                >>> a.remove_fileds_coo([10],[10],[10],[10])
                >>> a.data
                <Table length=2490>
                n                    name                      ra        dec      mag         dist
                int64                 str62                   float64    float64  float64     float64
                ----- -------------------------------------- ---------- --------- -------- -------------
                0 28655:NGC3034:NGC3034:09555243+6940469  148.96846 69.679703 -19.2115 4.70228466231
                1 42407:NGC4594:NGC4594:12395949-1137230 189.997894 -11.62307 -19.2974 3.65995814653
                2            --:---:---:12564369+2140575 194.182068 21.682659 -18.6564 3.87868462045
                3 41220:NGC4472:NGC4472:12294679+0800014 187.444992   8.00041 -20.9739 14.6955181559
                4 47404:NGC5194:NGC5194:13295269+4711429 202.469574 47.195259 -18.3042 2.49371678269
                5 10266:NGC1068:NGC1068:02424077-0000478   40.66988  -0.01329 -21.9158 18.4142043388
                6 42831:NGC4649:NGC4649:12434000+1133093 190.916702  11.55261 -21.3112 16.4529675514
                ...                                    ...        ...       ...      ...           ...
                2485              38897:NGC4173:NGC4173:---    183.089  29.20702 -19.2158 21.0457642607
                2486     --:---:SDSSJ022739.95-010913.2:---    36.9165  -1.15368  -18.494 21.2227303453
                2487     --:---:SDSSJ123616.69+260006.9:---     189.07  26.00194 -18.3719 25.2571985761
                2488     --:---:SDSSJ123626.88+255738.2:---    189.112  25.96063 -18.9699 20.7004974885
                2489     --:---:SDSSJ122250.38+155056.9:---     185.71  15.84916  -18.159 24.8768244613
                2490     --:---:SDSSJ124211.24+074016.0:---    190.547    7.6711 -18.0476 29.7521819633
                2491                     --:---:NGC5496:---    212.908  -1.15744 -18.3384 23.4249708242
                """
                if not data is None:
                        datain = data
                else:
                        datain = self.data
                        
                if datain is None:
                        self.logger.info ('Warning: galaxy data not found')
                        return                
                
                id2skip = self.ipix_in_box(np.array(ra), np.array(dec), np.array(fovra),
                        np.array(fovdec), nside, nest)                                           
                _theta, _phi = np.pi/2.-np.radians(datain['dec']),np.radians(datain['ra'])
                idgalaxies = hp.ang2pix(nside, _theta, _phi, nest=nest)
                ele = self.common_ele(idgalaxies, id2skip)
                if not data is None:
                        return datain[[not i for i in ele]]                
                else:                                                
                        self.data = datain[[not i for i in ele]]                
        
        def removep_file(self, split=' ', nest=False, data=None,
                         filename=None, filetype=None, nside=512, **kwargs):
                """remove fields by input file
                
                Parameters
                ------------         
                filename :         `string`
                  name of file
                  mush contain fovra and fovdec, so that define vertices
                filetype :         `string`
                  type of file, options: `txt` (plain text), `npz` (pickle file)
                wdir :             `string`
                  working directory for file       
                nside :      `int`
                  healpix nside parameter, must be a power of 2, less than 2**30
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead
                split :      `string`
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)
                
                Examples
                --------                
                >>> from kobe import galaxies
                >>> a = galaxies()
                >>> a.generate_pointings(limdec=[-20,90], limdist=[0,40])              
                >>> a.remove_pointings_file(filename='tmp.npz')
                """                
                _r = self.readlist(split=' ', filename=filename,
                                   filetype=filetype, keys=keys, **kwargs)         
                if _r is None: return
                
                # remove coolist
                return self.removep_coo(_r['ra'], _r['dec'], _r['fovra'],
                                _r['fovdec'], data=data, nside=nside, nest=nest)

        def savep(self, filename, filetype='npz',split=' ', data=None,
                           keys=['n','ra','dec','name','dist','mag'], **kwargs):
                """save tilings to a file
                
                Parameters
                ------------         
                filename :         `string`
                  name of file
                filetype :         `string`
                  type of file
                  options: `txt` (plain text), `npz` (pickle file), `fits` (astropy.fits)
                wdir :             `string`
                  working directory for file                        
                split :      `string`
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)

                Examples
                --------                
                >>> from kobe import galaxies
                >>> a = galaxies()
                >>> a.generate(limdec=[-20,90])               
                >>> a.save_pointings('tmp.npz', filetype='npz')
                """                  
                self.writelist(filename=filename, filetype=filetype, 
                               datain=data, keys=keys, **kwargs)                        

        def group_OB(self, obra, obdec):
                """group sub-galaxies, assign them with same index

                Parameters
                ----------  
                obra :         `float`           
                  tiling field of view in ra direction to group pointings
                obdec :        `float`           
                  tiling field of view in dec direction to group pointings
                """
                assert type(obra) in [float,int] and type(obdec) in [float,int]                
                
                if self.data is None:
                        self.logger.info ('Warning: pointings not parsed')
                        return

                ral = np.arange(min(self.data['ra']), max(self.data['ra']), obra)
                decl = np.arange(min(self.data['dec']), max(self.data['dec']), obdec)
                X, Y = np.meshgrid(ral, decl)
                                
                _n=0
                for z in [i for i in zip(X.flat,Y.flat)]:                        
                        _id1 = (self.data['ra']>=z[0])                        
                        _id2 = (self.data['ra']< z[0]+obra)
                        _id3 = (self.data['dec']>=z[1])
                        _id4 = (self.data['dec']< z[1]+obdec)
                        
                        self.data['n'][_id1 & _id2 & _id3 & _id4] = _n                        
                        if len(self.data[_id1 & _id2 & _id3 & _id4])>0:_n+=1
                        
        def preport(self, split=' '):       
                """make a report for data

                Parameters
                ----------   
                split     : `str`
                """        
                _report = self.dic2txt(split=split)
                if not _report in self.texts: self.texts += _report                        
        
        def galaxies2hpmap(self, tracer='l', nside=512, nest=False):
                """build a healpix map with KBGetGalaxies.data
        
                Parameters
                ------------
                tracer :      `str`
                  what was used to build the galaxy healpy map
                  options: [c]ounts, [l]uminosities
                nside :      `int`
                  healpix nside parameter, must be a power of 2, less than 2**30
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead

                returns
                ----------               
                hpmap :      array-like shape (Npix,)
                  healpix fits map for trigger localization probability (1d) 
                """                                 
                if self.data is None:
                        self.logger.info ('Warning: pointings not parsed')
                        return
                
                hpx = np.zeros(hp.nside2npix(nside))
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=nest)
                if tracer == 'l':
                        hpx[ipix] += 10**((-1)*(self.data['mag']/2.5))
                elif tracer == 'c':
                        hpx[ipix] += 1
                else:
                        self.logger.info ('Error: wrong tracer option')
                        return
                self.hpmapm = hpx/sum(hpx)

        def rankp(self, mode=1, nside=64, nest=False, threshold=1,sort=1):
                """rank pointings: rank from west to east

                Parameters
                ----------                   
                mode     : `int`
                  1. strict approach: follow above strategy strictly
                  2. adjacent approach: start from the westest point or highest probability point and arange the next pointing either adjacent or closest to the previous one
                """                
                assert 'data' in self.__dict__                     
                assert not self.data is None, 'Warning: tiling data not parsed'
                _data = self.data        
        
                # define start point                
                _idstart = np.argmin(_data['ra'])                
        
                if mode == 1:
                        ipx = np.argsort(_data['ra'])
            
                elif mode == 2:            
                        ipx = self.adjacent_move(_data['ra'],_data['dec'],
                                                 _idstart,threshold=1,sort=1)
                        
                else:
                        self.logger.info ('wrong mode')
                        return
                
                self.data = _data[ipx]
                self.data['n'] = np.arange(len(self.data['n']), dtype=int)
                
class pointings:
        """telescope: define the name and a pointing list for one telescope
        inherient from tilings and galaxies              
        """
        def __new__(cls, strategy, logger=None, **kwargs):                
                if cls is not pointings:
                        return super().__new__(cls)
                class_map = {'T': tilings, 'G': galaxies}
                cls = class_map[strategy]
                return cls.__new__(cls, logger=logger, **kwargs)
