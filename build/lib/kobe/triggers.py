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
from kobe import cookbook

__all__ = ('triggers')

class triggers:
    """triggers: Read and process trigger healpix map from a given source

    Parameters
    ----------
    hpmap :      array-like shape (Npix,) or (4, Npix)
      healpix fits map, including trigger localization probability (1d)
      and/or distance informations (3d).
    fithdr :     `dict`           
      healpix fits header, from fits file
    xmlinf :     `dict`           
      voevent informations, from xml file
    defconf :    `dict`           
      default configure, if any KBParseTriggers parameter was included in defconf dictionary, 
      then its default value would be applied
    logger :     `class`          
      logging object

    See Also
    --------
    KBGetTilings, KBGetGalaxies

    Examples
    --------
    see also https://github.com/saberyoung/kobe/blob/master/notebook/test_triggers.ipynb

    >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
    >>> a = KBParseTriggers()
    parse sources to obtain self.data
    Source options:           
        skymap      `string`         filename of healpix fits
        xmlfile     `string`         the filename of XML
        root        `class`          root element of the XML document parsed by lxml.etree
        coolist     `list`           [ra,dec,err], `float`, `float`, `float`
    >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
    >>> a.data_list()
    {'hpmap': array([[1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
        2.46590458e-10, 9.75865543e-10, 6.87730424e-10],
       [1.44844886e+02, 1.54923595e+02, 1.40961397e+02, ...,
        2.72896179e+02, 2.21541076e+02, 2.82526517e+02],
       [1.41256230e+02, 1.39624119e+02, 1.41143310e+02, ...,
        1.68491517e+02, 1.68595934e+02, 1.55918798e+02],
       [2.53084235e-05, 2.36136978e-05, 2.61189424e-05, ...,
        9.76686037e-06, 1.30767396e-05, 9.62544372e-06]]), 'xmlinf': None, 'fithdr': {'XTENSION': 'BINTABLE', 'BITPIX': 8, 'NAXIS': 2, 'NAXIS1': 32, 'NAXIS2': 3145728, 'PCOUNT': 0, 'GCOUNT': 1, 'TFIELDS': 4, 'TTYPE1': 'PROB', 'TFORM1': 'D', 'TUNIT1': 'pix-1', 'TTYPE2': 'DISTMU', 'TFORM2': 'D', 'TUNIT2': 'Mpc', 'TTYPE3': 'DISTSIGMA', 'TFORM3': 'D', 'TUNIT3': 'Mpc', 'TTYPE4': 'DISTNORM', 'TFORM4': 'D', 'TUNIT4': 'Mpc-2', 'MOC': True, 'PIXTYPE': 'HEALPIX', 'ORDERING': 'NESTED', 'COORDSYS': 'C', 'NSIDE': 512, 'INDXSCHM': 'IMPLICIT', 'OBJECT': 'G331903', 'REFERENC': 'https://gracedb.ligo.org/events/G331903', 'INSTRUME': 'H1,L1,V1', 'DATE-OBS': '2019-05-10T02:59:39.292500', 'MJD-OBS': 58613.12476032978, 'DATE': '2019-05-10T03:00:47.000000', 'CREATOR': 'BAYESTAR', 'ORIGIN': 'LIGO/Virgo', 'RUNTIME': 18.0, 'DISTMEAN': 268.8566049372629, 'DISTSTD': 108.0709050006497, 'LOGBCI': 0.6949211109947058, 'LOGBSN': 7.032293281836687, 'VCSVERS': 'ligo.skymap 0.1.6', 'VCSREV': '79504ec9fb1890fa91665bd69d7aa66cdaf11184', 'DATE-BLD': '2019-03-26T18:11:21', 'HISTORY': 'gwcelery worker -l info -n gwcelery-openmp-worker -Q openmp -c 1'}}
    >>> a.make_report()
    'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 OBJECT:G331903 INSTRUME:H1,L1,V1 CREATOR:BAYESTAR '
    >>> a.calc_area()   
    {0.5: 575.5456172035578, 0.9: 3463.7648737864856} 
    """    

    # Static version info
    version = 1.0

    # default key args
    defkwargs = {'wdir'        : '/tmp/',      # `str`   working directory
                 'savefits'    : None,         # `str`   save file or not
                 'nside'       : 512,          # `int`   healpix map default resolution
                 'coord'       : 'C',          # `str`   healpix map defualt coordinate system: G, C, E
                 'nest'        : False,        # `bool`  healpix map defualt ordering: nest or ring
                 'cls'         : [.5,.9],      # `list`  confidence levels for contours
                 'style'       : 'sms',        # `str`   report type
                 'clobber'     : True,         # `bool`  if fits parsed, redo or not
                 'skymap_key'  : 'skymap_fits',# key for skymap url, in XML
                 'obst_key'    : 'ISOTime',    # key for timeobs, in XML     
                 'radec_keys'  : ['C1','C2','Error2Radius'],  # key for coordinate and uncertainty, in XML
                 'keys_to_checkxml'  : ('GraceID', 'AlertType', 'Group', 'FAR',
                                        'Terrestrial', 'HasNS', 'HasRemnant',
                                        'BNS', 'BBH', 'NSBH', 'Instruments',
                                        'EventPage'),
                 'keys_to_checkfits' : ('DISTMEAN', 'DISTSTD', 'DATE-OBS',
                                        'MJD-OBS', 'OBJECT', 'INSTRUME',
                                        'CREATOR')
    }
    
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
        if not 'hpmap' in self.__dict__: self.hpmap = None
        if not 'hphdr' in self.__dict__: self.hphdr = None
        if not 'voinf' in self.__dict__: self.voinf = None        

    def _check(self):        
        """check if KBParseTriggers is parsed from one source or not  
        
        returns
        ----------
        res :        `bool`
        """
                
        if not self.hpmap is None: return True
        else: return False
        
    def url(self, url, **kwargs):
        """parser skymap url, and build KBParseTriggers.data  

        Parameters
        ----------
        url :        `string`         
          url of healpix fits
        savefits :   `string` or None    
          save healpix fits or not,
          if None, not save; 
          if a string, savefits would be used to construct the downloaded file.
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.
        nest :       `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        Examples
        --------       
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        """

        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        # parse fits_url
        if type(url) is str:
            self.fits(url,nest=kwargs['nest'])
        else:
            self.logger.info ('Error: wrong url format')
            return        
        if not self._check(): return
        
        # savefits
        if not kwargs['savefits'] is None:
            flag = self.download_skymap(kwargs['wdir'], kwargs['savefits'], url, kwargs['clobber'])            
            if flag == 1:
                self.logger.info ('download skymap via wget')
            elif flag == 2:
                self.logger.info ('download skymap via requests')
            elif flag == 0:
                self.logger.info ('Warning: fits exists already')
            else:
                self.logger.info (flag)
            
    def fits(self, skymap, **kwargs):
        """parser skymap file, and build KBParseTriggers.data  

        Parameters
        ----------
        skymap :        `string`         
          file name or url of healpix fits
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead
        
        Examples
        --------       
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()
        >>> a.url('S190510g_bayestar.fits.gz', nest=True)
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self._check() and not kwargs['clobber']:
            self.logger.info ('Warning: hpmap already parsed')
            return
        
        hpmap, hphdr = self.read_skymap(skymap, kwargs['nest'])
        self.hpmap = hpmap            
        self.hphdr = hphdr

    def xml(self, xml, **kwargs):       
        """parser voenevt XML file, and build KBParseTriggers.data  

        Parameters
        ----------
        xml :        `string`
          file name or url of voenevt XML file
        savefits :   `string` or None    
          save healpix fits or not,
          if None, not save; 
          if a string, savefits would be used to construct the downloaded file.
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.
        nside :      `int`
          healpix nside parameter, must be a power of 2, less than 2**30
        coord :      sequence of character, optional
          Either one of ‘G’, ‘E’ or ‘C’ to describe the coordinate system of the map, 
          or a sequence of 2 of these to rotate the map from the first to the 
          second coordinate system.

        Examples
        --------       
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()
        >>> a.xml('LVC#S190510g-2-Initial.xml')
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self._check() and not kwargs['clobber']:
            self.logger.info ('Warning: hpmap already parsed')
            return
        
        # read xml file, obtain root obj        
        root = self.parse_xml(xml)

        # parse root
        self.root(root, **kwargs)
        
    def root(self, root, **kwargs):        
        """parser root element, and build KBParseTriggers.data  

        Parameters
        ----------
        root :        `class`
          root element of the XML document, parsed by lxml.etree
        savefits :   `string` or None    
          save healpix fits or not,
          if None, not save; 
          if a string, savefits would be used to construct the downloaded file.
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.
        nside :      `int`
          healpix nside parameter, must be a power of 2, less than 2**30
        coord :      sequence of character, optional
          Either one of ‘G’, ‘E’ or ‘C’ to describe the coordinate system of the map, 
          or a sequence of 2 of these to rotate the map from the first to the 
          second coordinate system.

        Examples
        --------       
        check usage of pygcn (https://github.com/lpsinger/pygcn)
        
        >>> import gcn
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()
        >>> def handler(payload, root): 
        >>>     a.root(root)
        >>>     b = a.calc_area()
        >>>     print (b)
        >>> gcn.listen(handler=handler)
        """        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self._check() and not kwargs['clobber']:
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
        """use ra, dec, error box to build KBParseTriggers.data  
        
        Parameters
        ----------
        coo :        `list`
          [ra, dec, err]: ra, dec is the center of trigger 
          while err is the estimated error box region
        savefits :   `string` or None    
          save healpix fits or not,
          if None, not save; 
          if a string, savefits would be used to construct the downloaded file.
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.
        nside :      `int`
          healpix nside parameter, must be a power of 2, less than 2**30
        coord :      sequence of character, optional
          Either one of ‘G’, ‘E’ or ‘C’ to describe the coordinate system of the map, 
          or a sequence of 2 of these to rotate the map from the first to the 
          second coordinate system.
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        Examples
        --------                               
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()       
        >>> a.coo([30, -20, 10])
        """        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self._check() and not kwargs['clobber']:
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
            self.make_hpmap(**kwargs)            
                
    def make_hpmap(self, **kwargs):
        """generate and store healpix map, from KBParseTriggers.data['hpmap']

        Parameters
        ----------        
        tempfile :   `string` 
           healpix fits file name
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.        
        coord :      sequence of character, optional
          Either one of ‘G’, ‘E’ or ‘C’ to describe the coordinate system of the map, 
          or a sequence of 2 of these to rotate the map from the first to the 
          second coordinate system.
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead
        returns
        ----------
        res :        `bool`
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

        Parameters
        ----------   
        root :        `class`
          root element of the XML document, parsed by lxml.etree

        returns
        ----------
        res :        `bool`
          if True, url found (in case of GW preminary, initial, update);
          otherwise, url missing (for GW retraction, GRB, AMON, etc)
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

        Parameters
        ----------   
        root :        `class`
          root element of the XML document, parsed by lxml.etree
        nside :      `int`
          healpix nside parameter, must be a power of 2, less than 2**30

        returns
        ----------
        res :        `bool`
          if True, ra,dec,error box found
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
    
    def make_report(self, **kwargs):       
        """make a summary report on input source. 
        report is built from healpix fits header and XML file.
        specify keys in keys_to_checkxml, keys_to_checkfits

        Parameters
        ----------   
        style :        `string`
          Options: `sms`, `email`, `slack`

        returns
        ----------
        report :        `string`
          summary report

        Examples
        --------                               
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.make_report()
        'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 OBJECT:G331903 INSTRUME:H1,L1,V1 CREATOR:BAYESTAR '
        >>> a.keys_to_checkfits
        ('DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS', 'OBJECT', 'INSTRUME', 'CREATOR')
        >>> a.keys_to_checkfits=('DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS')
        >>> a.make_report()
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
        cls :         `list`
          list of confidence level, default: [.5, .9]

        returns
        ----------
        indexlist :   `dictionary`
          dictionary of healpix indices corresponding input C.L.

        Examples
        --------                               
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_contours()
        {0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.9: array([2414072, 2418168, 2416119, ...,  783552,  734398,  771264])}
        >>> a.calc_contours(cls=[.1,.5,.99])
        {0.1: array([2414072, 2418168, 2416119, ..., 1570953, 1573001, 1573000]), 0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.99: array([2414072, 2418168, 2416119, ..., 1309038, 1309052, 1309051])}
        """
                
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not self._check():
            self.logger.info ('Warning: hpmap not parsed')
            return
                
        if cookbook.is_seq_of_seq(self.hpmap):
            (hpx, hpd1, hpd2, hpd3) = self.hpmap
        elif cookbook.is_seq(self.hpmap):
            hpx = self.hpmap
        else:
            return
        
        indexlist = {}
        sortid = hpx.argsort()[::-1]
        sorted_hpx = hpx[sortid]
        cumulative_hpx = np.cumsum(sorted_hpx)
        for _cl in kwargs['cls']:
            if len(sorted_hpx[cumulative_hpx<_cl]) > 0:
                _limit = sorted_hpx[cumulative_hpx<_cl][-1]
                indexlist[_cl] = np.array([i for i in sortid if hpx[i] >= _limit])
        return indexlist
    
    def calc_area(self, **kwargs):       
        """calculate sky localization region area (unit in sq. deg) for different confidence level region of trigger map

        Parameters
        ----------   
        cls :         `list`
          list of confidence level, default: [.5, .9]

        returns
        ----------
        arealist :   `dictionary`
          dictionary of area corresponding input C.L.

        Examples
        --------                               
        >>> from kobe.pipeline.KBParseTriggers import KBParseTriggers
        >>> a = KBParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_area()
        {0.5: 575.5456172035578, 0.9: 3463.7648737864856}
        >>> a.calc_area([.1,.5,.99])
        {0.1: 36.351906008208665, 0.5: 575.5456172035578, 0.99: 11508.39446313552}
        """
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])   
        _ilist = self.calc_contours(**kwargs)
        if _ilist is None: return
        
        if cookbook.is_seq_of_seq(self.hpmap):
            (hpx, hpd1, hpd2, hpd3) = self.hpmap            
        elif cookbook.is_seq(self.hpmap):
            hpx = self.hpmap
        else:
            return
        
        nside = hp.get_nside(hpx)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2
        _alist = {}
        for _cl in kwargs['cls']:
            _area = areapix*len(_ilist[_cl])
            _alist[_cl] = _area
        return _alist    
    
    def parse_xml(self, xmlfile):
        """parse xmlfile via lxml

        Parameters
        ----------   
        xmlfile :      `string`
          file name of XML        

        returns
        ----------
        root :        `class`          
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

        Parameters
        ----------   
        _fits :    `string`
           healpix fits or fits url
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        returns
        ----------
        tmap :     array-like shape (Npix,) or (4, Npix)
          Either an array representing a map, or a sequence of 4 arrays 
            representing map, distance mean, variance, normalization           
        header : `dictionary`
          header of healpix fits
        """
               
        try:
            # read 3D trigger healpix map (for GW)            
            tmap, header = hp.read_map(_fits,field=[0, 1, 2, 3], nest=nest,h=True)
        except:
            # if failed, read 2D map            
            try:
                tmap, header = hp.read_map(_fits, nest=nest, h=True)                
            except:
                return None, None
        return tmap, dict(header)
                
    @staticmethod
    def download_skymap(wdir, tmpfile, url, clobber):            
        """look up URL of sky map, download sky map, and parse FITS file

        Parameters
        ----------   
        tmpfile :   `string`    
          file name of healpix fits
        wdir :       `string`
          working directory for downloaded healpix file if savefits is not None.
        url :        `string`
          url of healpix fits         
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
