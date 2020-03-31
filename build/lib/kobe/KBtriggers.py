#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/triggers.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
import numpy as np
import healpy as hp
import astropy.time
import astropy.coordinates        
import astropy.units as u
from kobe import visualization, cookbook

__all__ = ('triggers')

class triggers(visualization):
    """
    triggers: Read and process trigger healpix map from a given source
    class inherit from kobe.KBplotter.visualization

    Parameters
    ----------
    wdir      : `str`   
       working directory       
    savefits  : `str`   
       save healpix fits file (if parsed trigger via voevent) or not
    nside     : `int`   
       defualt healpix map resolution, applied if no healpix map obtained and will generate one
    coord     : `str`   
       defualt healpix map coordinate system
       options: `G`, `C`, `E`
    nest      : `bool`  
       default healpix map ordering
       options: True (for `nest` ordering) or False (`ring`)    
    cls       : `list` 
       confidence levels for contours
       e.g. cls=[.5, .9] would construct 50% and 90% CL contours
    style     : `str`
       when making report, specify the stype
       current options: `email`, `slack`, `sms`
       can be added as you want in kobe.KBtriggers.triggers.make_trigger_report()
    clobber   : `bool`  
       if any products already available in self, clobber=True will redo it
    skymap_key: `str`
       once user intend to build healpix map via voevent xml file, 
       which key is respond for skymap url
       default: `skymap_fits`
    obst_key  : `str`
       which key is used to build the trigger time in voevent XML file
       default: `ISOTime`
    radec_keys: `list`
       ra, dec, uncertainty keys that would be used for the constrction of healpix map
       default: ['C1','C2','Error2Radius']
    keys_to_checkxml  : `list`
       which parameters from voevent XML would be reported
    keys_to_checkfits : `list`
       which parameters from healpix fits header XML would be reported   

    See Also
    --------
    KBplotter, KBcirculate

    Examples
    --------
    see https://github.com/saberyoung/kobe/blob/master/notebook/test_triggers.ipynb
    """    

    version = 1.0
    '''
    class kobe.KBtrigger version
    '''
    
    defkwargs = {'wdir'        : '/tmp/',
                 'savefits'    : None,
                 'nside'       : 512,
                 'coord'       : 'C',
                 'nest'        : False,
                 'cls'         : [.5,.9],
                 'style'       : 'sms',
                 'clobber'     : True,
                 'skymap_key'  : 'skymap_fits',
                 'obst_key'    : 'ISOTime',
                 'radec_keys'  : ['C1','C2','Error2Radius'],
                 'keys_to_checkxml'  : ('GraceID', 'AlertType', 'Group', 'FAR',
                                        'Terrestrial', 'HasNS', 'HasRemnant',
                                        'BNS', 'BBH', 'NSBH', 'Instruments',
                                        'EventPage'),
                 'keys_to_checkfits' : ('DISTMEAN', 'DISTSTD', 'DATE-OBS',
                                        'MJD-OBS', 'OBJECT', 'INSTRUME',
                                        'CREATOR')
    }
    '''
    default parameter settings
    '''

    hpmap = None    
    '''
    array-like shape (Npix,)  healpix fits map for 2d sphere probability
    '''

    hpd1 = None
    '''
    array-like shape (Npix,)  healpix fits map for distance mean in each direction (for CBC GW only)
    '''

    hpd2 = None
    '''
    array-like shape (Npix,)  healpix fits map for distance variance in each direction (for CBC GW only)
    '''
    
    hpd3 = None
    '''
    array-like shape (Npix,)  healpix fits map for distance normalization in each direction (for CBC GW only)
    '''

    hphdr = None
    '''
    `dict`  healpix fits header
    '''

    voinf = None
    '''
    `dict`  voevent informations
    '''
                           
    def __init__(self, logger=None, **kwargs):
        """
        initialize trigger object with input logger and parameter settings
        
        Parameters
        ----------
        logger :     `class`
            logging object

        **kwargs :     
            check kobe.triggers.defkwargs, kobe.visualization.defkwargs
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
        self.kwargs = {**kwargs, **visualization.defkwargs}
        
    def url(self, url, **kwargs):
        """parser hpmay via skymap url

        Parameters
        ----------
        url :        `string`         
          url of healpix fits

        Other Parameters
        ----------
        `savefits`, `wdir, `nest`

        Examples
        --------       
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.hpmap
        array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
               2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
        >>> a.hpd1
        array([144.84488626, 154.92359484, 140.96139722, ..., 272.89617891,
               221.54107632, 282.52651738])
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        # parse fits_url
        if type(url) is str:
            self.fits(url,nest=kwargs['nest'])
            if not self.hpmap is None: return
        else:
            self.logger.info ('Error: wrong url format')
            return
        
        # savefits
        if not kwargs['savefits'] is None:
            flag = self.download_skymap(kwargs['wdir'], kwargs['savefits'],
                                        url, kwargs['clobber'])            
            if flag == 1:
                self.logger.info ('download skymap via wget')
            elif flag == 2:
                self.logger.info ('download skymap via requests')
            elif flag == 0:
                self.logger.info ('Warning: fits exists already')
            else:
                self.logger.info (flag)
            
    def fits(self, skymap, **kwargs):
        """parser hpmap via skymap file 

        Parameters
        ----------
        skymap :        `string`         
          file name or url of healpix fits

        Other Parameters
        ----------
        `nest`
        
        Examples
        --------       
        >>> from kobe import triggers
        >>> a = triggers(nest=True)
        >>> a.fits('S190510g_bayestar.fits.gz')
        >>> a.hpmap
        array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
               2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpmap is None and not kwargs['clobber']:
            self.logger.info ('Warning: hpmap already parsed')
            return
        
        self.hpmap, self.hpd1, self.hpd2, self.hpd3, self.hphdr = \
            self.read_skymap(skymap, kwargs['nest'])
        
    def xml(self, xml, **kwargs):       
        """parser hpmap via voenevt XML file

        Parameters
        ----------
        xml :        `string`
          file name or url of voenevt XML file

        Other Parameters
        ----------
        `savefits`, `wdir`, `nside`, `coord`

        Examples
        --------       
        >>> from kobe import triggers
        >>> a = triggers()
        # build healpix fits by downloading one via the url in voevent XML file
        >>> a.xml('LVC#S190510g-2-Initial.xml')
        array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
               2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
        # build healpix fits by construction one via the ra,dec,loc parameters in voevent XML file
        >>> a.xml('notebook/icecube_ehe.xml')
        >>> a.hpmap
        array([ 7.63510930e-07,  7.54514778e-07,  7.43099315e-07, ...,
              -8.13487380e-08, -7.78317088e-08, -8.04872986e-08])
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpmap is None and not kwargs['clobber']:
            self.logger.info ('Warning: hpmap already parsed')
            return
        
        # read xml file, obtain root obj        
        root = self.parse_xml(xml)

        # parse root
        self.root(root, **kwargs)
        
    def root(self, root, **kwargs):        
        """parser hpmap via a root element

        Parameters
        ----------
        root :        `class`
          root element of the XML document, parsed by lxml.etree

        Other Parameters
        ----------
        `savefits`, `wdir`, `nside`, `coord`

        Examples
        --------              
        >>> from kobe import triggers
        >>> a = triggers(clobber=True, savefits=True)        
        # use pygcn to moniror ligo server
        >>> import gcn
        >>> def handler(payload, root): 
        >>>     a.root(root)
        >>>     b = a.calc_area()
        >>>     print (b)
        >>> gcn.listen(handler=handler)
        """        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpmap is None and not kwargs['clobber']:
            self.logger.info ('Warning: hpmap already parsed')
            return

        try:            
            self.voinf = {elem.attrib['name']: elem.attrib['value'] for elem in root.iterfind('.//Param')}   
        except:
            return

        # 1. find skymap url
        url = self.parse_root_skymap(root, kwargs['skymap_key'])        
        if not url is None:
            self.url(url, **kwargs)
            return
        
        # 2. find ra,dec,error box        
        coo = self.parse_root_coo(root, kwargs['nside'], kwargs['radec_keys'], kwargs['obst_key'])
        if not coo is None:
            self.coo(coo, **kwargs)
            
    def coo(self, coo, **kwargs):        
        """parse hpmap with ra, dec, error box
        
        Parameters
        ----------
        coo :        `list`
          [ra, dec, err]: ra, dec is the center of trigger 
          while err is the estimated error box region

        Other Parameters
        ----------
        `savefits`, `wdir`, `nside`, `coord`, `nest`
       
        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        # generate a skymap with center at ra=30 deg, dec=-20 deg, with 10 deg sky area
        >>> a.coo([30, -20, 10])
        >>> a.hpmap
        array([8.92728175e-11, 8.92185286e-11, 8.91040539e-11, ...,
               6.12155463e-11, 6.08690152e-11, 6.10351951e-11])
        """        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpmap is None and not kwargs['clobber']:
            self.logger.info ('Warning: hpmap already parsed')
            return

        _ra,_dec,_loc = coo
        _radius = np.radians(_loc) # from deg to radians
        _theta, _phi = np.pi/2.-np.radians(_dec),np.radians(_ra)        
        _pmap = np.zeros(hp.nside2npix(kwargs['nside']))
        _index = hp.ang2pix(kwargs['nside'], _theta, _phi, nest=kwargs['nest'])
        _pmap[_index]+=1
        _pmap=hp.smoothing(_pmap,fwhm=_radius,sigma=None)                
        self.hpmap = _pmap/sum(_pmap)

        # savefits
        if not kwargs['savefits'] is None:
            self.make_hpmap_coo(**kwargs)            
                
    def savefits(self, **kwargs):
        """save self.hpmap to '%s/%s' % (kwargs['wdir'], kwargs['savefits'])

        Parameters
        ----------        
        `savefits`, `wdir`, `coord`, `nest`         
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        skymap = '%s/%s' % (kwargs['wdir'], kwargs['savefits'])
        try:
            hp.write_map(skymap, self.hpmap, nest=kwargs['nest'], coord=kwargs['coord'],
                         extra_header=self.hphdr, overwrite=True)            
            self.logger.info ('make healpix map %s' % skymap)            
        except:
            self.logger.info ('Warning: failed to make healpix map')

    def parse_root_skymap(self, root, skymap_key):       
        """parse root object, try find healpix fits map url
        used by kobe.triggers.root
        """        

        # skymap url
        xmlinfos = {elem.attrib['name']:
                    elem.attrib['value']
                    for elem in root.iterfind('.//Param')}        
        if skymap_key in xmlinfos:
            self.logger.info ('obtain hpmap via %s' % skymap_key)
            skymap_url = xmlinfos[skymap_key]
        else:
            self.logger.info ('Warning: no %s found in voevent'%skymap_key)
            skymap_url = None
        return skymap_url
            
    def parse_root_coo(self, root, nside, radec_keys, obst_key):        
        """parse root object, try find ra,dec and localization of trigger
        used by kobe.triggers.root
        """
        
        _coolist = []
        for _kk in radec_keys:
            _val = root.findtext('.//%s'%_kk)
            if not _val is None:
                try:
                    _coolist.append(float(_val))                    
                except:
                    self.logger.info ('Warning: %s not float, skip'%_kk)
            else:
                self.logger.info ('Warning: no %s found in voevent'%_kk)
                
        if len(_coolist) == 3:
            _ra,_dec,_loc = _coolist
        else:
            return None
        
        """
        make header            
        """        
        # obs time
        timeobs = root.findtext('.//%s'%obst_key)

        # obs mjd        
        mjdobs = astropy.time.Time(timeobs, format='isot', scale='utc').mjd

        # read classifier
        try:
            tid = root.attrib['ivorn']
        except:
            tid = None

        # header
        hdr = [('CREATOR','KOBE'),
               ('IVORN',tid),
               ('RA',_ra),
               ('DEC',_dec),
               ('ERROR-BOX',_loc),
               ('NSIDE',nside),
               ('MJD-OBS',mjdobs),
               ('DATE-OBS',timeobs)]
        
        self.hphdr = dict(hdr)
        return _coolist
    
    def make_trigger_report(self, **kwargs):       
        """make a summary report for input source. 
        report is built from healpix fits header and XML file.
        specify keys in keys_to_checkxml, keys_to_checkfits

        Parameters
        ----------   
        `style`

        returns
        ----------
        report :        `string`
          summary report

        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.make_report()
        'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 OBJECT:G331903 INSTRUME:H1,L1,V1 CREATOR:BAYESTAR '
        >>> a.defkwargs['keys_to_checkfits']
        ('DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS', 'OBJECT', 'INSTRUME', 'CREATOR')        
        >>> a.make_report(keys_to_checkfits=['DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS'])
        'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 '
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        _dict =  {}
        for _kk,_vv in zip(['keys_to_checkxml', 'keys_to_checkfits'],
                           [self.voinf, self.hphdr]):
            for el in kwargs[_kk]:
                if el in _vv:                    
                    _dict[el] = _vv[el]
                    
        alist = self.calc_area(**kwargs)
        if not alist is None:
            for el in alist:                
                _dict[el] = alist[el]
                
        _report = ''
        for el in _dict:
            try:
                _vv = '%.2f' % _dict[el]
            except:
                _vv = '%s' % _dict[el]                
            if kwargs['style'] == 'sms':                
                _report += '%s:%s '%(el,_vv)
            elif kwargs['style'] == 'slack':
                _report += '`%s`:%s\n'%(el,_vv)
            elif kwargs['style'] == 'email':
                _report += '#\t%s:\t%s\n'%(el,_vv)                
        return _report

    def calc_contours(self, **kwargs):        
        """calculate indices located in specific confidence level region of trigger map

        Parameters
        ----------   
        `cls`
        
        returns
        ----------
        indexlist :   `dictionary`
          dictionary of healpix indices corresponding input C.L.

        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_contours()
        {0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.9: array([2414072, 2418168, 2416119, ...,  783552,  734398,  771264])}
        >>> a.calc_contours(cls=[.1,.5,.99])
        {0.1: array([2414072, 2418168, 2416119, ..., 1570953, 1573001, 1573000]), 0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.99: array([2414072, 2418168, 2416119, ..., 1309038, 1309052, 1309051])}
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        
        indexlist = {}
        sortid = self.hpmap.argsort()[::-1]
        sorted_hpx = self.hpmap[sortid]
        cumulative_hpx = np.cumsum(sorted_hpx)
        for _cl in kwargs['cls']:
            if len(sorted_hpx[cumulative_hpx<_cl]) > 0:
                _limit = sorted_hpx[cumulative_hpx<_cl][-1]
                indexlist[_cl] = np.array([i for i in sortid if self.hpmap[i] >= _limit])
        return indexlist
    
    def calc_area(self, **kwargs):       
        """calculate sky localization region area (unit in sq. deg) for different confidence level region of trigger map

        Parameters
        ----------   
        `cls`

        returns
        ----------
        arealist :   `dictionary`
          dictionary of area corresponding input C.L.

        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_area()
        {0.5: 575.5456172035578, 0.9: 3463.7648737864856}
        >>> a.calc_area([.1,.5,.99])
        {0.1: 36.351906008208665, 0.5: 575.5456172035578, 0.99: 11508.39446313552}
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])   
        _ilist = self.calc_contours(**kwargs)
        if _ilist is None: return        
        
        nside = hp.get_nside(self.hpmap)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2
        _alist = {}
        for _cl in kwargs['cls']:
            _area = areapix*len(_ilist[_cl])
            _alist[_cl] = _area
        return _alist    

    def calc_prob_loc(self, ra=None, dec=None, theta=None, phi=None,
                      fovra=None, fovdec=None, **kwargs):
        """calculate localization probability for specific regions
        if fov is None, return probability of input coordinates
        otherwise, calculate sum of one or a serios of vertices

        Parameters
        ----------   
        ra :   None, `int`, `float`-or sequence
          ra, unit in deg
        dec :  None, `int`, `float`-or sequence
          dec, unit in deg
        theta :  None, `int`, `float`-or sequence
          healpix theta, unit in radians
        phi :   None, `int`, `float`-or sequence
          healpix phi, unit in radians
        fovra :  None, `int`, `float`-or sequence
          set the width for a tiling
        fovdec :  None, `int`, `float`-or sequence
          set the height for a tiling

        returns
        ----------
        arealist :   `dictionary`
          dictionary of area corresponding input C.L.

        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        # calculate prob for point at ra=1, dec=1, and ra=10, dec=-20
        >>> a.calc_prob_loc(ra=[1,10], dec=[1,-20])
        array([7.36631063e-07, 8.50340927e-06])
        # calculate total probs for vertice 1: ra=1, dec=1, fovra=1, fovdec=1 and vertice 2: ra=10, dec=-20, fovra=10, fovdec=10
        >>> a.calc_prob_loc(ra=[1,10], dec=[1,-20],fovra=[1,10],fovdec=[1,10])
        [5.8577539241927334e-05, 0.042806591674792185]
        """
         
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        nside = hp.get_nside(self.hpmap)

        _list = True
        if not ra is None and not dec is None:            
            if cookbook.is_seq(ra) and cookbook.is_seq(dec):
                ra, dec = np.array(ra), np.array(dec)
            else:
                _list = False
            theta, phi = np.pi/2.-np.radians(dec),np.radians(ra)
                
        elif not theta is None and not phi is None:
            if cookbook.is_seq(theta) and cookbook.is_seq(phi):
                theta, phi = np.array(thata), np.array(phi)
            else:
                _list = False
            ra, dec = np.degrees(phi),np.degrees(np.pi/2.-theta)
            
        else:
            self.logger.info ('Error: wrong coo input')
            return

        if not fovra is None and not fovdec is None:
            if cookbook.is_seq(fovra) and cookbook.is_seq(fovdec):
                if _list:
                    if len(fovra) != len(ra) or len(fovdec) != len(dec):
                        self.logger.info('Error: ra, dec and fovra, fovdec length not the same')
                        return
                else:
                    self.logger.info('Error: ra, dec has no length')
                    return
            else:
                if _list:
                    fovra = np.zeros(len(ra)) + fovra
                    fovdec = np.zeros(len(dec)) + fovdec
                else:
                    ra,dec,fovra,fovdec = [ra], [dec], [fovra], [fovdec]
            probs = []
            for _ra,_dec,_fovw,_fovh in zip(ra, dec, fovra, fovdec):
                ipix_poly=(self.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,kwargs['nest']))
                _probs = self.hpmap[ipix_poly].sum()
                probs.append(_probs)                
        else:
            probs = hp.get_interp_val(self.hpmap, theta, phi,
                                      nest=kwargs['nest'], lonlat=False)            
        return probs

    def calc_prob_dis(self, ra=None, dec=None, theta=None, phi=None,
                      details=False, limdist=1000, **kwargs):        
        """calculate trigger probability for specific directions, based on distance informations
        for CBC GW, the distance mean, variance and normalization at each direction is provided
        by assuming the limiting distance of one telescope on target, e.g. 200 Mpc on KNe,
        thus, triggers.calc_prob_dis will calculate the Gaussian probability from 0 to 200 Mpc,
        obtain the largest probability within such region
                        
        Parameters
        ------------         
        ra :   None, `int`, `float`-or sequence
          ra, unit in deg
        dec :  None, `int`, `float`-or sequence
          dec, unit in deg
        theta :  None, `int`, `float`-or sequence
          healpix theta, unit in radians
        phi :   None, `int`, `float`-or sequence
          healpix phi, unit in radians       
        limdist :  `float`
          limiting distance of telescope on targeting sources, unit in Mpc
        details :  `bool`
          if True, will return distance mean, var, norm at such direction
            
        Examples
        --------                
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_prob_dis(ra=[1,10], dec=[1,-20])
        [3.6100349724294025e-06, 3.089517808764708e-06]
        >>> a.calc_prob_dis(ra=[1,10], dec=[1,-20], details=True)
        (array([233.86346031, 273.32868364]), array([121.25377326, 112.45704833]), array([1.44320589e-05, 1.14500239e-05]))
        """
                                
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.hpd1 is None or self.hpd2 is None or self.hpd3 is None:            
            self.logger.info ('Warning: no distance information')
            return
        nside = hp.get_nside(self.hpmap)
        pixarea = hp.nside2pixarea(nside, degrees=True)

        if not ra is None and not dec is None:            
            if cookbook.is_seq(ra) and cookbook.is_seq(dec):
                ra, dec = np.array(ra), np.array(dec)            
            theta, phi = np.pi/2.-np.radians(dec),np.radians(ra)
                
        elif not theta is None and not phi is None:
            if cookbook.is_seq(theta) and cookbook.is_seq(phi):
                theta, phi = np.array(thata), np.array(phi)           
            ra, dec = np.degrees(phi),np.degrees(np.pi/2.-theta)
            
        else:
            self.logger.info ('Error: wrong coo input')
            return
        
        from scipy.stats import norm        
        r = np.linspace(0, limdist)        
        ipix = hp.ang2pix(nside,theta,phi,nest=kwargs['nest'])                
        dmu, dsigma, dnorm = self.hpd1[ipix], self.hpd2[ipix], self.hpd3[ipix]
        if details:
            return dmu, dsigma, dnorm
        probl = [dnorm * norm(dmu, dsigma).pdf(rr)/pixarea for rr in r]
        probs = [max(ii) for ii in list(map(list, zip(*probl)))]               
        return probs    
    
    def parse_xml(self, xmlfile):
        """parse xmlfile via lxml
        used by kobe.triggers.xml
        """
         
        try:
            from lxml import etree
            self.logger.info ("running with lxml.etree")
        except ImportError:
            try:
                # Python 2.5
                import xml.etree.cElementTree as etree
                self.logger.info ("running with cElementTree on Python 2.5+")
            except ImportError:
                try:
                    # Python 2.5
                    import xml.etree.ElementTree as etree
                    self.logger.info ("running with ElementTree on Python 2.5+")
                except ImportError:
                    try:
                        # normal cElementTree install
                        import cElementTree as etree
                        self.logger.info ("running with cElementTree")
                    except ImportError:
                        try:
                            # normal ElementTree install
                            import elementtree.ElementTree as etree
                            self.logger.info ("running with ElementTree")
                        except ImportError:
                            self.logger.info ("Failed to import ElementTree from any known place")
                            return None
        tree = etree.parse(xmlfile)
        root = tree.getroot()        
        return root
        
    @staticmethod
    def read_skymap(_fits,nest):
        """read healpix fits or fits url, obtain map and header
        """
               
        try:
            # read 3D trigger healpix map (for GW)            
            tmap, header = hp.read_map(_fits,field=[0, 1, 2, 3], nest=nest,h=True)
            return tmap[0], tmap[1], tmap[2], tmap[3], dict(header)
        except:
            # if failed, read 2D map            
            try:
                tmap, header = hp.read_map(_fits, nest=nest, h=True)
                return tmap, None, None, None, dict(header)
            except:
                return None, None, None, None, None            
                
    @staticmethod
    def download_skymap(wdir, tmpfile, url, clobber):            
        """look up URL of sky map, download sky map, and parse FITS file    
        """
        
        skymap = '%s/%s' % (wdir, tmpfile)
        if os.path.exists(skymap) and not clobber:
            return 0
        try: 
            import wget
            flag = 1
        except ImportError as e:
            try:
                import requests,tempfile,shutil
                flag = 2
            except ImportError as e:
                return e
        if flag == 1:            
            wget.download(url, out=skymap)            
        else:                      
            # Send HTTP request for sky map   
            response = requests.get(url, stream=True)   

            # Raise an exception unless the download succeeded (HTTP 200 OK)
            response.raise_for_status()

            # Create a temporary file to store the downloaded FITS file
            with tempfile.NamedTemporaryFile() as tmpfile:
                # Save the FITS file to the temporary file
                shutil.copyfileobj(response.raw, tmpfile)
                tmpfile.flush()

                # Uncomment to save FITS payload to file       
                shutil.copy(tmpfile.name, skymap)                  
        return flag

    @staticmethod
    def ipix_in_box(ra,dec,width,height,nside,nest):
        """finding the healpix indices of a given box        
        """
                                                
        v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = \
            triggers.vertices(ra, dec, width, height)
        ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                     [v1_dec, v2_dec, v4_dec, v3_dec])                
        theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
        phi = np.deg2rad(ra_vertices)
        xyz = hp.ang2vec(theta, phi)                           
        if cookbook.is_seq(ra) and cookbook.is_seq(dec) and \
           cookbook.is_seq(width) and cookbook.is_seq(height):
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
