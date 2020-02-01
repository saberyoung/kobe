#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/PstGetGalaxies.py
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
from astroquery.vizier import Vizier
import logging

__all__ = ('PstGetGalaxies')

class PstGetGalaxies():
        """PstGetGalaxies: Generate, read and process galaxies

        Parameters
        ----------
        name :       `string`
          telescope name, if none, will be constructed with coordinates
        ra :         `list`           
          galaxies ra list, default: None
        dec :        `list`           
          galaxies dec list, default: None
        distance :        `list`           
          galaxies distance list, default: None
        mag :        `list`           
          galaxies magnitude list, default: None
        gname :      `list`           
          galaxies name list, default: None      
        index :      `list`           
          galaxies index list, default: None
        defconf :    `dict`           
          default configure, if any PstGetGalaxies parameter was included in defconf dictionary, 
          then its default value would be applied
        logger :     `class`
          logging object

        See Also
        --------
        PstGetTilings, PstGetGalaxies

        Examples
        --------
        see also https://github.com/saberyoung/pstools/blob/master/notebook/test_galaxies.ipynb
        
        >>> from pst.pipeline.Galaxies import PstGetGalaxies
        >>> a = Galaxies()
        >>> a.data
        """ 

        # Static version info
        version = 1.0
   
        def __init__(self, name=None, ra=None, dec=None,
                     distance=None, gname=None, mag=None,
                     index=None, defconf=None, logger=None):
                """
                initialize
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set galaxy properties ----- #
                self.data   =   {'ra':ra, 'dec':dec, 'n':index,
                                 'dist':distance, 'name':gname,
                                 'mag':mag}

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
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> dict = {'catalog': 'GWGC', 'limdist': [0,40]}
                >>> a.run_config(dict)
                """
                self.conf = {
                        'catalog':    'GLADE',
                        'filter':     'B',
                        'size':       -1,
                        'limra':      [0,360.],        # `range`   [0,360.]
                        'limdec':     [-89,89],        # `range`   [-89.,89.]
                        'limdist':    [0,1000],        # `range`   [0,1000]
                        'limmag':     [-18,-99],       # `range`   [-10,-99]                        
                        'nside':      512,             # `int`     default healpix resolution
                        'nest':       False,     # `bool`  healpix map defualt ordering: nest or ring
                        'wdir':       '/tmp/',         # `str`     working directory
                        'filetype':   'npz',           # `str`     options: `txt`, `npz`
                        'filename':   'pst_galaxies',  # `str`     tilings file name
                        'obstime':    None,            #
                        'lat':        24.625,          # la parma
                        'lon':        70.403,          # la parma
                        'alt':        2635             # la parma
                }
        
                if defconf is None: return        
                for _k in self.conf.keys():
                        if _k in defconf.keys():
                                self.conf[_k] = defconf[_k]
                        else:
                                self.logger.info ('### Warning: use default value for %s'%_k)
                
        def generate(self,limra=None,limdec=None,
                     catalog=None,filtro=None,size=None,
                     limdist=None,limmag=None):
                """generate galaxies by querying Vizier
                
                Parameters
                ----------
                limra :         `range`
                  galaxies: ra range, default: [0, 360]
                limdec :        `range`           
                  galaxies: dec range, default: [-90, 90]
                limdist :       `range`           
                  galaxies: distance range, default: [0, 1000]
                  Notice: depending on the network, it will take time to download Vizier galaxies, 
                    we recommend you to download your galaxy dataset, and use the store option to remain the list, 
                    so that to be faster for the next time
                limmag :        `range`           
                  galaxies: absolute magnitude range, default: [-18, -99]
                catalog :       `string`                  
                  Vizier galaxy catalog, options: `GWGC`, `GLADE`, default: GLADE
                filtro :        `string`           
                  magnitude filter. `B` and `K` are available for `GLADE` while only `B` available for `GLADE`
                  default: B
                size :      `float`           
                  size of the querying galaxies, -1 for the full query. default: -1

                Examples
                --------                
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90], limdist=[0,40])
                >>> a.data
                {'n': array([   0,    1,    2, ..., 2489, 2490, 2491]), 'name': array(['28655:NGC3034:NGC3034:09555243+6940469', '42407:NGC4594:NGC4594:12395949-1137230', '--:---:---:12564369+2140575', ..., '--:---:SDSSJ122250.38+155056.9:---', '--:---:SDSSJ124211.24+074016.0:---', '--:---:NGC5496:---'], dtype='<U62'), 'ra': array([148.96846 , 189.997894, 194.182068, ..., 185.71    , 190.547   ,  212.908   ]), 'dec': array([ 69.679703, -11.62307 ,  21.682659, ...,  15.84916 ,   7.6711  , -1.15744 ]), 'mag': array([-19.2115, -19.2974, -18.6564, ..., -18.159 , -18.0476, -18.3384]), 'dist': array([ 4.70228466,  3.65995815,  3.87868462, ..., 24.87682446, 29.75218196, 23.42497082])}
                >>> a.astrotab()
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
                
                _hp = self.checkdata()
                if _hp: return
                
                if limra is None:    limra   =  self.conf['limra']
                if limdec is None:   limdec  =  self.conf['limdec']
                if limdist is None:  limdist =  self.conf['limdist']
                if limmag is None:   limmag  =  self.conf['limmag']
                if catalog is None:  catalog =  self.conf['catalog']
                if filtro is None:   filtro  =  self.conf['filter']
                if size is None:     size    =  self.conf['size']                

                # specify columns
                if catalog == 'GLADE':
                        if not filtro in ['B', 'K']:
                                self.logger.info ('### Error: wrong filters for GLADE')
                                return
                        catid, columns = 'VII/281', \
                                ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                                 'Dist', 'PGC', 'GWGC', 'HyperLEDA', '2MASS']
                elif catalog == 'GWGC':
                        if not filtro in ['B']:
                                self.logger.info ('### Error: wrong filters for GWGC')
                                return
                        catid, columns = 'VII/267', \
                                ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                                 'Dist', 'Name']
                else:
                        self.logger.info ('### Error: wrong galaxy catalogs')
                        return

                # download catalog with vizier                
                v = Vizier(columns=columns, \
                           column_filters={columns[0]:'%s..%s'%(str(limra[0]),str(limra[1])),\
                                           columns[1]:'%s..%s'%(str(limdec[0]),str(limdec[1])),\
                                           columns[2]:'%s..%s'%(str(limmag[0]),str(limmag[1])),\
                                           columns[3]:'%s..%s'%(str(limdist[0]),str(limdist[1]))
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

                self.data = {'n':       np.arange(len(catalogs[columns[0]])),
                             'name':    np.array(_name),
                             'ra':      np.array(catalogs[columns[0]]),
                             'dec':     np.array(catalogs[columns[1]]),
                             'mag':     np.array(catalogs[columns[2]]),
                             'dist':    np.array(catalogs[columns[3]])}                

        def checkdata(self):
                """check if PstGetGalaxies is parsed or not
        
                returns
                ----------
                res :        `bool`                  
                """
                
                if not self.data['ra'] is None and \
                   not self.data['dec'] is None and \
                   not self.data['dist'] is None and \
                   not self.data['mag'] is None:
                        self.logger.info ('galaxies has already been parsed')
                        return True
                else:
                        return False

        def hpmap(self, nside=None, nest=None):
                """build a healpix map with PstGetGalaxies.data
        
                Parameters
                ------------
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
                 
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if nside is None: nside = self.conf['nside']
                if nest is None: nest = self.conf['nest']
                
                hpx = np.zeros(hp.nside2npix(nside))
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=nest)
                hpx[ipix] += 10**((-1)*(self.data['mag']/2.5))                
                hpx = hpx/sum(hpx)
                return hpx
        
        def read(self, filename=None, filetype=None, wdir=None):
                """read galaxies from file

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
        
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.read(filename='pst_galaxies', filetype='npz', wdir='./')
                """
                
                _res = self.readfile(filename=filename, filetype=filetype, wdir=wdir)
                if not _res is None:                        
                        self.data = {'ra':_res['ra'], 'dec':_res['dec'],'n':_res['n'],
                                'dist':_res['dist'], 'mag':_res['mag'], 'name':_res['name']}
                
        def readfile(self, filename=None, filetype=None, wdir=None):
                """read galaxies from file

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
                  Galaxies: `n`, `ra`, `dec`, `dist`, `mag`, `name`
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
                                for _k in ['ra','dec','dist','mag','n','name']:
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
                                for _k in ['ra','dec','dist','mag','n','name']:
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
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdist=[60,100])               
                >>> a.astrotab()
                <Table length=23991>
                n                       name                         ra        dec        distance     mag
                int64                    str82                      float64    float64      float64    float64
                ----- -------------------------------------------- ---------- ---------- ------------- --------
                0              37670:ESO217-017:ESO217-017:---    179.586  -50.97639  63.717247595 -20.4313
                1                  --:---:---:20181914-7051317 304.579773 -70.858833 62.2100694769 -22.9563
                2                    17625:NGC1961:NGC1961:---    85.5199   69.37866 63.8005161226 -23.0541
                3 64041:ESO185-054:ESO185-054:20032702-5556498  300.86261 -55.947193 68.8713202645 -22.3002
                4                  --:---:---:00574891+3021083  14.453818  30.352314 78.3041870859 -22.4849
                5       43423:NGC4709:NGC4709:12500394-4122554 192.516449 -41.382072 66.6945103643 -22.6005
                6       15406:NGC1600:NGC1600:04313985-0505099  67.916077  -5.086102 67.2166182679 -21.8334
                7         49025:IC4329:IC4329:13490531-3017452 207.272141 -30.295889 63.8224292971 -22.4189
                8       43924:NGC4782:NGC4782:12543569-1234069 193.648743  -12.56859  69.376260513 -22.2661
                9         62407:IC4765:IC4765:18471814-6319521 281.825623  -63.33115 70.2414517955  -21.838
                10       10302:NGC1060:NGC1060:02431504+3225300  40.812687  32.425011 79.1940957915 -22.2765
                11 57649:ESO137-008:ESO137-008:16154609-6055071 243.942078 -60.918663 63.0029911687 -21.9668
                12 57612:ESO137-006:ESO137-006:16150386-6054261 243.766098 -60.907261 79.8947610861 -22.1526
                13       19622:NGC2258:NGC2258:06474618+7428546 101.942436  74.481842 65.6418095541 -21.7609
                ...                                          ...        ...        ...           ...      ...
                23976                        --:---:ESO217-006:---    177.444   -52.1222 87.2376818255 -18.6435
                23977                        --:---:ESO576-047:---     201.04    -17.902 98.4997238141 -18.3072
                23978                        --:---:ESO383-019:---     203.45  -37.73178 97.6881556667 -20.9692
                23979                        --:---:ESO383-059:---     205.01  -32.89598  96.216373912 -18.6662
                23980                        --:---:ESO383-084:---    207.186  -35.80501 98.0074344629 -18.7263
                23981                        --:---:ESO383-093:---    207.684  -36.75278 99.6797792736  -20.203
                23982                          --:---:UGC10489:---    249.335   62.74141  88.497194846 -18.5246
                23983                          --:---:UGC10557:---    252.061   13.90453 69.7802767525 -18.4587
                23984                        --:---:ESO102-006:---    262.235  -66.15067 90.7612652378 -19.6995
                23985                        --:---:ESO399-008:---    299.471  -32.46556 82.8182312366 -19.7206
                23986                          --:---:UGC11869:---    329.926   44.35817 85.4620306264 -18.2489
                23987                          --:---:UGC11899:---    331.187   44.94891 83.2020728228 -18.7407
                23988                          --:---:UGC12062:---    337.729   39.44567  80.035795991 -18.5164
                23989                          --:---:UGC12761:---     356.14    6.04077 93.7045229755 -18.4388
                23990                        --:---:ESO215-011:---    164.564  -49.42383 82.4167960688 -19.6501
                """

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                
                return Table([self.data['n'],
                              self.data['name'],
                              self.data['ra'],
                              self.data['dec'],
                              self.data['dist'],
                              self.data['mag']],
                             names=('n', 'name', 'ra', 'dec', 'distance', 'mag'))
        
        def save(self, filetype=None, filename=None, wdir=None):
                """save galaxies to a file
                
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
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90])               
                >>> a.save()
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
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
                                 dist=self.data['dist'],mag=self.data['mag'],
                                 name=self.data['name'],n=self.data['n'])
                                
                elif self.conf['filetype'] == 'txt':
                        cachefile = '%s/%s.txt' % (self.conf['wdir'], self.conf['filename'])
                        if os.path.exists(cachefile): os.remove(cachefile)

                        # store to cachefile
                        ww = open(cachefile,'w')                       
                        ww.write('# n ra dec distance mag name \n')
                        for _ra,_dec,_dist,_mag,_n,_name in zip(self.data['ra'], self.data['dec'],
                                                self.data['dist'],self.data['mag'],
                                                self.data['n'],self.data['name']):
                                ww.write('%i %.5f %.5f %.2f %.2f %s \n'%(_n,_ra,_dec,_dist,_mag,_name))
                        ww.close()
                else:
                        self.logger.info ('### Warning: filetype %s unknown'%self.conf['filetype'])

        def calc_prob_loc(self, triggerobj, nest=None):
                """calculate the localization probability of each galaxy in their center pixel
                
                Parameters
                ------------         
                triggerobj :       `class`
                  PstParseTriggers object
                nest :       `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead                   
                
                Examples
                --------                
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90], limdist=[0,40])               
                >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
                >>> b = PstParseTriggers()   
                >>> b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
                >>> c = a.calc_prob_loc(b)
                >>> c
                <Table length=2492>
                n            prob
                int64        float64
                ----- ----------------------
                0 1.5637598103792824e-05
                1 4.1564742774336236e-05
                2 7.4099969227877845e-06
                3  1.392453548154686e-05
                4 3.1442444023962773e-07
                5 3.7145227221084495e-05
                6 2.1804973957569932e-05
                7  2.326217672248188e-07
                8  6.785513402446836e-06
                9  9.436367017384702e-08
                10 3.3040075395909226e-07
                11  1.395602640846228e-08
                12  1.958021450430296e-06
                13 2.2640017389577904e-07
                ...                    ...
                2477  2.232319447165003e-07
                2478  5.814795979191668e-07
                2479    4.8723962399549e-07
                2480 1.3827155101346376e-06
                2481 2.8209811142971617e-06
                2482  3.430237663418057e-08
                2483 3.3050864016173866e-09
                2484  6.652237415193616e-05
                2485  3.168978960448268e-07
                2486  8.649814370731974e-06
                2487   8.80995743012443e-07
                2488   8.80995743012443e-07
                2489 1.0677603786094531e-06
                2490  1.904402472728465e-05
                2491  0.0001926440173553161
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if not nest is None: self.conf['nest'] = nest
                
                from pst.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                elif is_seq(triggerobj.data['hpmap']):
                        hpx = triggerobj.data['hpmap']
                else: return
                
                probs, ns = [], []
                nside = hp.get_nside(hpx)
                pixarea = hp.nside2pixarea(nside, degrees=True)
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=self.conf['nest'])
                probs = hpx[ipix]/pixarea                        
                return Table([np.arange(len(probs)), probs], names=('n', 'prob'))  
        
        def calc_prob_dis(self, triggerobj, nest=None):
                """calculate each galaxy the probability from distance distribution.
                
                Parameters
                ------------         
                triggerobj :       `class`
                  PstParseTriggers object
                nest :             `bool`
                  healpix ordering options: 
                  if True, healpix map use `nest` ordering, otherwise, use `ring` instead
                
                Examples
                --------                
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90], limdist=[0,40])               
                >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
                >>> b = PstParseTriggers()   
                >>> b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
                >>> c = a.calc_prob_dis(b)
                >>> c
                <Table length=2492>
                n            prob
                int64        float64
                ----- ----------------------
                0  4.689546195431277e-07
                1 1.1885425467637767e-06
                2  7.236225284574737e-08
                3 1.8311021499842933e-07
                4  2.299045888552183e-08
                5  3.970894399243602e-06
                6 1.7237244486275295e-07
                7 1.5800110539396612e-06
                8  9.533719932024894e-08
                9  1.694104788141493e-05
                10 3.6634249776590484e-08
                11 5.3503570784002356e-08
                12                    0.0
                13 1.6129515983493087e-07
                ...                    ...
                2477 1.2071687148609212e-08
                2478  2.785388678971533e-08
                2479  6.522168308072937e-08
                2480  4.655106737253541e-05
                2481  5.309869829433903e-07
                2482  6.040797494995478e-07
                2483  8.661872354220442e-09
                2484  4.199977807128782e-06
                2485 3.6786813759222278e-09
                2486  2.787555103923971e-06
                2487  1.508035341020697e-08
                2488 1.3361634353665521e-08
                2489 1.9511270857343716e-08
                2490 2.7607735557675464e-07
                2491 3.0980960021459507e-06
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if not nest is None: self.conf['nest'] = nest
                
                from pst.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                else:
                        return

                from scipy.stats import norm        
                nside = hp.get_nside(hpx)
                pixarea = hp.nside2pixarea(nside, degrees=True)
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=self.conf['nest']) 
                dmu, dsigma, dnorm = hpd1[ipix], hpd2[ipix], hpd3[ipix]
                probs = dnorm * norm(dmu, dsigma).pdf(self.data['dist'])/pixarea
                return Table([self.data['n'], probs], names=('n', 'prob'))

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
                if lat is None: lat = self.conf['lat']
                if lon is None: lon = self.conf['lon']
                if alt is None: alt = self.conf['alt']     
                observatory = astropy.coordinates.EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg, dec=self.data['dec']*u.deg)               

                # Alt/az reference frame at observatory, now
                frame = astropy.coordinates.AltAz(obstime=obstime, location=observatory)

                # Transform grid to alt/az coordinates at observatory, now
                altaz = radecs.transform_to(frame)
                return altaz, obstime, observatory
                
        def calc_airmass(self, obstime=None, lat=None, lon=None, alt=None):
                """calculate airmass for galaxies

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
                  aismass of each galaxy

                Examples
                -------- 
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90], limdist=[0,40]) 
                >>> b=a.calc_airmass()
                >>> b
                <Table length=2492>
                n         airmass
                int64       float64
                ----- -------------------
                0    2.45822474769971
                1  -1.582820782124962
                2  -2.418465393396415
                3 -2.1552326830446353
                4  -5.836324948920255
                5  1.1808382514756013
                6  -2.100571797667892
                7  -4.409531346475804
                8 -2.3322078594904623
                9 -22.456682403975346
                10  -5.530047411614138
                11  1.4202164017530534
                12   3.035017007059533
                13   3.514106192534303
                ...                 ...
                2478  -8.050072007646094
                2479  -3.244213057431077
                2480 -1.3238910756227278
                2481   5.675834562202663
                2482   1.285876078526185
                2483   4.697281278442855
                2484 -1.1460986096538353
                2485  -4.839581546710565
                2486  1.2271304406718495
                2487 -3.2252172351830937
                2488 -3.2156676633639614
                2489 -2.6912267582548473
                2490 -1.9772720228546403
                2491  -1.245843683922017
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                return Table([self.data['n'], altaz.secz], names=('n', 'airmass'))  
        
        def calc_sun(self, obstime=None, lat=None, lon=None, alt=None):
                """calculate sun height for galaxies

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
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90], limdist=[0, 40]) 
                >>> b=a.calc_sun(obstime='2020-01-01 00:00:00')
                >>> b
                <Latitude -27.0241102 deg>
                """

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the sun, now?
                sun_altaz = astropy.coordinates.get_sun(obstime).transform_to(altaz)
                return sun_altaz.alt

        def calc_solar(self, sobj, obstime=None, lat=None, lon=None, alt=None):
                """calculate distance of galaxies to solar object

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
                >>> from pst.pipeline.PstGetGalaxies import PstGetGalaxies
                >>> a = PstGetGalaxies()
                >>> a.generate(limdec=[-20,90], limdist=[0, 40]) 
                >>> b=a.calc_solar('moon')
                >>> b
                <Table length=2492>
                n          dist
                deg
                int64      float64
                ----- ------------------
                0  53.33330588792869
                1  88.24625042763911
                2  79.55175108487516
                3  78.30869490502155
                4  77.32031077621345
                5  68.09252455916292
                6  80.20545198639483
                7 62.159401885764616
                8  76.97208445032689
                9  53.00344701561276
                10  60.55196604922355
                11 54.169964482303634
                12 104.36704371919019
                13  33.49649104707532
                ...                ...
                2478  66.04212038717846
                2479   83.8518979831337
                2480 138.01191140250333
                2481 120.79226149839266
                2482 54.790569376554906
                2483  96.55882206244557
                2484 119.08533624379777
                2485   67.7534238448892
                2486  71.98254120864843
                2487  73.72821267289586
                2488  73.77667723826598
                2489  73.96900697494519
                2490  81.28711502508183
                2491 105.26879312115406
                >>> b=a.calc_solar('sun')
                >>> b
                <Table length=2492>
                n          dist
                deg
                int64      float64
                ----- ------------------
                0 51.959195379172066
                1  76.49962734958095
                2  68.86971807113785
                3  66.63391392017586
                4  70.21262580349818
                5  79.54835437287645
                6  68.73822577642416
                7 50.405703769869724
                8  65.49131587811843
                9 42.359346813019016
                10  48.84373745219334
                11  59.52268646001445
                12 113.02831343004696
                13 21.967900070925094
                ...                ...
                2478 56.910345907873406
                2479  75.87621789437932
                2480  130.8055874507337
                2481 132.58578867970064
                2482 61.845696164577355
                2483 100.66624864878354
                2484 108.00146970397446
                2485  57.57861000515114
                2486   83.4599566267763
                2487 63.331472768688236
                2488  63.37667486721704
                2489  62.64938693196699
                2490  69.64205481246496
                2491  93.62311346982547
                """
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if not sobj in ['earth', 'sun', 'earth-moon-barycenter', 'moon',
                                'mercury','venus','mars','jupiter', 'saturn',
                                'uranus','neptune']:
                        self.logger.info ('### Error: sobj wrong')
                        return
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg, dec=self.data['dec']*u.deg)
                
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the solar object, now
                moonloc = astropy.coordinates.get_body(sobj, obstime, observatory).transform_to(altaz)
                return Table([self.data['n'], radecs.separation(moonloc)], names=('n', 'dist'))

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
