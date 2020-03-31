#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : kobe/KButils.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

# Standard library
import os
import sys

# third-party module
import astropy.time
from astropy.table import vstack, Column, Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import healpy as hp

__all__ = ['utils']

class utils:
    """
    provide various generic useful functions kobe objects   

    Parameters
    ----------    
    wdir      : str   
       working directory.
       default: ./.
    clobber   : bool  
       if any product has already been parsed, 
       set clobber as True will remove and redo the calculation.
       default: False
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

    Notes
    ----------   
    *utils* is the basic **KOBE** class
    """
        
    version = 1.0
    '''
    kobe class version
    '''
    
    defkwargs = {'wdir'           : './',
                 'clobber'        : False,
                 # plotter
                 'rot_theta'      : 0,
                 'rot_phi'        : 0,
                 'rot_psi'        : 0,                 
                 'num'            : 1,                 
    }
    '''
    default parameter settings
    '''
    
    exkwargs = {'readmap'    : ['dtype', 'nest', 'partial', 'hdu', 'verbose', 'memmap'],
                'writemap'   : ['nest', 'coord', 'overwrite'],
                'healpixview': ['coord', 'unit', 'xsize', 'title', 'xsize', 'ysize',
                                'nest', 'min', 'max', 'flip', 'remove_dip',
                                'remove_mono', 'gal_cut', 'format', 'format2',
                                'cbar', 'notext', 'norm', 'hold', 'margins', 'sub'],
                'matplotlib' : ['figsize', 'dpi','facecolor', 'edgecolor', 'frameon', 'clear'],
                'graticule'  : ['dpar', 'dmer', 'coord', 'local'],
                'projplot'   : ['ms', 'label', 'marker', 'color', 'coord'],
                'projtext'   : ['label', 'color', 'coord'],
                'savefig'    : ['dpi', 'transparent', 'pad_inches', 'format'],
                'coord'      : ['coord'],
                'pixtool'    : ['nest'],
    }
    '''
    set parameters for external package calls: healpy, matplotlib, etc
    '''
    
    def ckpython(self):
        '''check python version
        '''
        if sys.version_info>(3,0,0): return 3
        else: return 2
            
    def ckdir(self, **kwargs):
        ''' check working directory,
        see if it's exists, and if readable/writable        

        Parameters
        ----------    
        wdir      : str   
           working directory        
        '''
        kwargs = self.setkeys(kwargs)
        wdir = kwargs['wdir']
        assert type(wdir) is str, 'directory should be a string'
        
        try:
            from pathlib import Path
            output = Path(wdir)
            if output.exists():
                assert output.is_dir()
            else:
                output.mkdir()
        except:            
            if not os.path.exists(wdir):
                os.mkdir(wdir)
                
        if os.access(wdir, os.R_OK | os.W_OK):
            return wdir
        else:
            self.logger.info ('Error: insufficient priority to read and write in %s'%wdir)            

    def setkeys(self, kwargs):
        '''set kobe key arguments:
        if one kwarg provide, set it, otherwise use default

        Parameters
        ----------    
        wdir      : str   
           working directory   
        '''
        for _key in self.defkwargs:
            kwargs.setdefault(_key, self.defkwargs[_key])

        for _key in kwargs:
            if _key in self.defkwargs: continue
            _ifex = False
            for _ex in self.exkwargs:
                if _key in self.exkwargs[_ex]: _ifex = True
            assert _ifex, 'Error: key %s not recognized' % _key            
        return kwargs
    
    def getkeys(self, kwargs, _class=None):
        '''get specific kwargs for kobe and external function calls.
        e.g. _class='matplotlib' will return all matplotlib parameters
        '''
        kwargs = self.setkeys(kwargs)
        if _class is None:
            # select kobe parameters
            kwargs = {key: value for key, value in
                      kwargs.items() if key in self.defkwargs}            
        else:
            # select external parameters            
            assert _class in self.exkwargs
            kwargs = {key: value for key, value in
                      kwargs.items() if key in self.exkwargs[_class]}
        return kwargs
                
    @staticmethod
    def is_seq(o):
        """Check if the object is a sequence

        Parameters
        ----------
        o : any object
           The object to check
    
        Returns
        -------
        is_seq : bool, scalar
           True if *o* is a sequence, False otherwise
        """
        return hasattr(o, "__len__")

    def is_seq_of_seq(self, o):
        """Check if the object is a sequence of sequences

        Parameters
        ----------
        o : any object
           The object to check
    
        Returns
        -------
        is_seq_of_seq : bool
           True if *o* is a sequence of sequences, False otherwise.
        """
        if not self.is_seq(o):
            return False
        for s in o:
            if not self.is_seq(s):
                return False
        return True
            
    def ipix_in_box(self,ra,dec,width,height,nside,nest):
        """finding the healpix indices of a given box   

        Parameters
        ----------
        ra,dec,width,height : sequence, or single number
           the center and size of a vertice box
        nside :      int
            healpix nside parameter, must be a power of 2, less than 2**30
        nest :       bool
            healpix ordering options: 
            if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        Returns
        -------
        ipix_in_box : list
           a sequence of healpix indices    
        """
                                                
        v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = \
            self.vertices(ra, dec, width, height)
        ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                     [v1_dec, v2_dec, v4_dec, v3_dec])                
        theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
        phi = np.deg2rad(ra_vertices)
        xyz = hp.ang2vec(theta, phi)                           
        if self.is_seq(ra) and self.is_seq(dec) and \
           self.is_seq(width) and self.is_seq(height):
            ipix_fov_box = []
            for _xyz in xyz:
                ipix_fov_box.append(hp.query_polygon(nside, _xyz, nest=nest))
        else:
            ipix_fov_box = hp.query_polygon(nside, xyz, nest=nest)
        return ipix_fov_box
    
    @staticmethod
    def vertices(ra,dec,fovw,fovh):
        """finding the vertices of a FoV by giving:
        the central location, i.e. ra, dec in deg) 
        and the FoV size: fovw, fovh in deg
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

    def common_ele(self, list1, list2, nullv=-99):
        """ find the common elements of 2 lists
        
        Parameters
        ----------   
        list1 :   sequence
                  an index list
        list2 :   sequence
                  an index list

        returns
        ----------   
        mask :     sequence
             a bool list, that indeces of common elements is True, otherwise False
        """               
        list1 = self.flatten(list1, nullv=nullv)
        list2 = self.flatten(list2, nullv=nullv)
        assert not list1 is None and not list2 is None
        _mask = np.in1d(list1, list2)        
        return _mask
    
    def flatten(self, _list, nullv=-99):
        """ flatten a sequence if it's a sequence of sequence
        
        Parameters
        ----------   
        _list :   sequence
                  a list
        nullv :   float
                  if the lengths of sequences are not the same,
                  fill in the shorter ones with nullv, so that they have the same length

        returns
        ----------   
        mask :     sequence
             a bool list, that indeces of common elements is True, otherwise False
        """ 
        if self.is_seq_of_seq(_list):
            # check length of each elements
            llist = np.unique([len(ll) for ll in _list])
            lmax = max([len(ll) for ll in _list])                                
            if len(llist) > 1:
                # arrange each elements to same length
                _listm = []
                for ll in _list:
                    lla = [nullv]*(lmax - len(ll))
                    ll = np.append(ll, lla)
                    _listm.append(ll)
                _list = _listm                        
            # flatten to 1d
            _listf = [item for sublist in _list for item in sublist]
        elif self.is_seq(_list): _listf = _list
        else: return
        return _listf
            
    def readlist(self, filename=None, filetype=None,
                 split=' ', keys=None, **kwargs):
        """ read data from a file
        
        Parameters
        ----------   
        filename :         string
                  name of file, basename is needed (without dirname)
        filetype :         string
                  type of file
                  options: `txt` (plain text), `npz` (pickle file), `fits` (astropy.fits)
        split :      string
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)
        keys :   sequence
                  keys that read from file and construct data
        
        Other Parameters
        ----------   
        wdir :         string
                  working directory

        returns
        ----------   
        data :     dictionary
             a dictionary
        """ 
        if keys is None: return

        wdir = self.ckdir(**kwargs)
        assert not wdir is None
        filename = '%s/%s' % (wdir, filename)
            
        data = {}  
        if filetype == 'fits':
            if not '.fits' in filename: filename += '.fits'
            assert os.path.exists(filename)
            return Table.read(filename, format='ascii')
        
        if filetype == 'txt':
            if not '.txt' in filename: filename += '.txt'
            assert os.path.exists(filename)
            _keys, _data = {}, {}            
            for nll, ll in enumerate(open(filename).readlines()):
                if ll[0] == '#': # read header
                    ll = ll.replace('#','').replace('\n','')            
                    _nk = 0
                    for _k in ll.split(split):
                        if len(_k) != 0:
                            _keys[_nk] = _k
                            _data[_k] = []
                            _nk += 1
                else: # read data
                    _nk = 0
                    for _k in ll.split(split):                    
                        try:  float(_k)                            
                        except:  continue
                        if _nk in _keys:
                            _data[_keys[_nk]].append(float(_k))
                            _nk += 1                      
            for _k in keys:
                if _k in _data.keys(): data[_k] = _data[_k]                        
        
        if filetype == 'npz':
            if not '.npz' in filename: filename += '.npz'
            assert os.path.exists(filename)               
            try: _npdata = np.load(filename)
            except: return            
            for _k in keys:               
                try:
                    data[_k] = _npdata['data'][_k]
                except:
                    pass                
        return data

    def writelist(self, datain=None, filename=None, filetype=None,
                  split=' ', keys=None, **kwargs):
        """ write data to a file
        
        Parameters
        ----------   
        filename :         string
                  name of file, basename is needed (without dirname)
        filetype :         string
                  type of file
                  options: `txt` (plain text), `npz` (pickle file), `fits` (astropy.fits)
        split :      string
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)
        keys :   sequence
                  keys that read from file and construct data        
        wdir :         string
                  working directory
        clobber   : bool  
                  if any products already available in self, clobber=True will redo it 
        """
        if keys is None: return
        kwargs = self.setkeys(kwargs)        
        wdir = self.ckdir(**kwargs)        

        assert not wdir is None
        assert not filename is None
        filename = '%s/%s' % (wdir, filename)
        
        if datain is None:            
            assert 'data' in self.__dict__
            assert not self.data is None 
            datain = self.data
            
        if datain is None:
            self.logger.info ('Error: no data')
            return
            
        if filetype == 'fits':
            if not '.fits' in filename: filename += '.fits'
            if os.path.exists(filename) and not kwargs['clobber']:
                return
            datain.write(filename, overwrite=True)
                        
        if filetype == 'npz':
            if not '.npz' in filename: filename += '.npz'
            if os.path.exists(filename) and not kwargs['clobber']:
                return
            np.savez(filename, data=datain)
                        
        if filetype == 'txt':            
            if not '.txt' in filename: filename += '.txt'
            if os.path.exists(filename) and not kwargs['clobber']:
                return
            ww = open(filename,'w')                
            ww.write('# %s \n'%(split.join(keys)))
            _klist = list(datain.keys())            
            for ii in range(len(datain[_klist[0]])):
                _text=''
                for _key in keys:
                    if _key in _klist:
                        _text += '%s%s'%(datain[_key][ii],split)
                    else:
                        self.logger.info ('Error: key %s not found'%_key)                
                ww.write('%s \n'%_text)
            ww.close()

    def dic2txt(self, split=' '):
        """ read dictionary like object, generate standard texts    
        
        Parameters
        ----------   
        split :      string
                  if filetype is txt, the string defined to split each lines
                  default: ' ' (space)           
        """
        assert self.data
        assert not self.data is None       
        keys = list(self.data.keys())
        _txt = '%s \n'%(split.join(keys))
        for ii in range(len(self.data[keys[0]])):            
            for _key in keys:                
                _txt += '%s%s'%(self.data[_key][ii],split)          
            _txt += '\n'
        return _txt
    
    def spherical_angle(self,ra1,dec1,ra2,dec2):
        """ calculate spherical angle
        
        Parameters
        ----------   
        ra1, dec2 :      float, float
                  should be a point
        ra, dec2 :      float, float
                  could be a point or a sequence of points
        """        
        assert not self.is_seq(ra1) and not self.is_seq(dec1)        
        if self.is_seq(ra2):
            assert len(ra2) == len(dec2)
            ra2, dec2 = np.array(ra2), np.array(dec2)
        else:
            assert not self.is_seq(dec2)

        c1 = SkyCoord(ra1, dec1, unit='deg')
        c2 = SkyCoord(ra2, dec2, unit='deg')
        dist = c1.separation(c2)                
        return dist

    def adjacent_move(self,ra,dec,ral,decl,threshold=5,sort=1):
        '''
        sort = 1: clockwise
        sort = 2: reverse clockwise
        '''                   
        _dist = self.spherical_angle(ra,dec,ral,decl)                        
        assert min(_dist) != 0, 'remove duplicated pointings'
                            
        _sortid = np.argsort(_dist)
        al, il = [], []
        for _nd in _sortid:               
            _dd = _dist[_nd]                
            if _dd.deg > threshold + min(_dist).deg:
                break
            if ral[_nd] > ra and decl[_nd] > dec:
                al.append(0)
                il.append(_nd)
            if ral[_nd] > ra and decl[_nd] == dec:
                al.append(1)
                il.append(_nd)
            if ral[_nd] > ra and decl[_nd] < dec:
                al.append(2)
                il.append(_nd)
            if ral[_nd] == ra and decl[_nd] < dec:
                al.append(3)
                il.append(_nd)
            if ral[_nd] < ra and decl[_nd] < dec:
                al.append(4)
                il.append(_nd)
            if ral[_nd] < ra and decl[_nd] == dec:
                al.append(5)
                il.append(_nd)
            if ral[_nd] < ra and decl[_nd] > dec:
                al.append(6)
                il.append(_nd)
            if ral[_nd] == ra and decl[_nd] > dec:
                al.append(7)
                il.append(_nd)
                    
        assert sort in [1,2,3] and len(al) > 0
        if sort == 1:
            _idrr = il[np.argmin(al)]
        elif sort == 2:
            _idrr = il[np.argmax(al)]
        elif sort == 3:
            for _a in [1, 5, 7, 3]:
                if _a in al:
                    _idrr = _a
                    break
            for _a in [0, 4, 6, 2]:
                if _a in al:
                    _idrr = _a                                           
        return _idrr

    def arrange_filter(self, data=None, add=False, filters='r'):
        """set telescope name

        Parameters
        ----------
        filter:
                telescope filter
                if multiple, split with comma        
        """
        if data is None:
            assert self.pointings.data
            data = self.pointings.data
                            
        ff=[]
        for _nf, _filt in enumerate(filters.split(',')):
            if _nf == 0:                                
                _data = data                                
            else:                                
                _data = vstack([_data, data], join_type='outer')
            ff += [_filt]*len(data)                        
        ffc = Column(ff, name='filter')

        if add:
            self.pointings.data = _data
            self.pointings.data.add_column(ffc)
        else:
            return ffc
        
    def arrange_time(self, data=None, add=False, ovhtime=120, exptime=45):
        """set telescope observe time

        Parameters
        ----------        
        exptime:        
                estimated time for telescope exposure per frame + field slewing
                unit in second               
        ovhtime:       
                estimated time for telescope overhat
                unit in second
        """
        if data is None:
            assert self.pointings.data
            data = self.pointings.data
                    
        # arrange times
        tt = Column(np.zeros(len(data['n'])), name='time')
                
        _time = 0
        for _nn in np.unique(data['n']):
            
            _idx = np.where(data['n']==_nn)            
            _time += ovhtime            

            _tt = []
            for _ii in range(len(data[_idx])):
                
                _time += exptime
                _tt.append(_time)

            tt[_idx] = _tt
            
        if add:
            assert self.pointings.data
            self.pointings.data.add_column(tt)
        else:            
            return tt

        
    def radecs(self, ra=None, dec=None, theta=None, phi=None, hms=False):
        '''read ra dec
        '''
        
        if hms: ra,dec = self.hms2deg(ra,dec)        
        if not ra is None and not dec is None:            
            if self.is_seq(ra):
                assert len(ra) == len(dec)                             
                ra, dec = np.array(ra), np.array(dec)
                
        elif not theta is None and not phi is None:
            if self.is_seq(theta):
                assert len(theta) == len(phi)
                theta, phi = np.array(thata), np.array(phi)                                
            dec, ra = -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
                        
        else:                        
            return
        
        if not self.is_seq(ra):
            ra, dec = np.array([ra]), np.array([dec])
        return ra, dec    

    def hms2deg(self, ra, dec):
        
        _p = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))                
        return _p.ra.deg, _p.dec.deg
        
    def deg2hms(self, ra, dec):
        
        _p = SkyCoord(ra=ra, dec=dec, unit='deg')                                  
        _rahms = '%.2i:%.2i:%.3f'%(_p.ra.hms[0],abs(_p.ra.hms[1]),abs(_p.ra.hms[2]))
        _dechms = '%.2i:%.2i:%.3f'%(_p.dec.dms[0],abs(_p.dec.dms[1]),abs(_p.dec.dms[2]))
        return _rahms, _dechms

    @staticmethod
    def gcn_server(_server,_port):
        '''get host and port so that pygcn could use for monitoring triggers.  
      
        Parameters
        ----------  
        _server: 
            local, eApps, Atlantic_2, Atlantic_3, Linode
        _port: 
            public, LVC, AMON
        '''        
        assert _server in ['local', 'eApps', 'Atlantic_2', 'Atlantic_3', 'Linode']            
        assert _port in ['public', 'LVC', 'AMON']            
        
        if _server == 'local':        
            return '127.0.0.1',8099

        if _server == 'eApps':            
            if _port == 'public':
                return '68.169.57.253', 8099
            if _port == 'LVC':
                return '68.169.57.253', 8092
            if _port == 'AMON':
                return '68.169.57.253', 8096

        if _server == 'Atlantic_2':            
            if _port == 'public':
                return '209.208.78.170', 8099
            if _port == 'LVC':
                return '209.208.78.170', 8092
            if _port == 'AMON':
                return '209.208.78.170', 8096
            
        if _server == 'Atlantic_3':           
            _server,_port = '45.58.43.186',8099

        if _server == 'Linode':            
            if _port == 'public':
                return '50.116.49.68', 8099
            if _port != 'public':
                return '50.116.49.68', 8096

        return _server,_port
