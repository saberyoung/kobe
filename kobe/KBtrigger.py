#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : kobe/KBtrigger.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import os
import logging
import numpy as np
import healpy as hp
import random
import astropy.time
import astropy.coordinates        
import astropy.units as u
from astropy.table import Column, Table
from kobe import vtrigger, circulate, candidates

__all__ = ['trigger']

class trigger(vtrigger, circulate, candidates):
    """Read and process trigger healpix map from a given source.
    inherit from kobe.vtrigger and kobe.circulate;
    inherited by kobe.schedule.

    Parameters
    ----------
    url :        `string`         
       url of healpix fits
    hpfile :        `string`         
          file name or url of healpix fits
    vofile :        `string`
          file name or url of voenevt XML file
    root :        `class`
          root element of the XML document, parsed by lxml.etree
    coo :        `list`
          [ra, dec, err]: ra, dec is the center of trigger,
          while err is the estimated error box region
    wdir      : `str`   
       working directory 
    savefits  : `str`
       if set to None, will not save fits.
       Otherwise, the input string would be set as the filename 
       for the saved healpix fits file, if trigger map was achived remotely.    
    clobber   :  `bool`
       if any product has already been parsed, 
       set clobber as True will remove and redo the calculation.              
    """    
    
    defkwargs = {'url'         : None,
                 'hpfile'      : None,
                 'vofile'      : None,
                 'root'        : None,
                 'coo'         : None,
                 'savefits'    : None,                            
    }
    '''
    default parameter settings
    '''
    defkwargs = {**defkwargs,
                 **vtrigger.defkwargs,
                 **circulate.defkwargs}

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
        """

        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # ----- set keys ----- #
        self.defkwargs = self.setkeys(kwargs)

    def pygcn_serve(self, server='local', port='public'):
        """serve root elements via pygcn

        Parameters
        ----------        
        root :        `class`
          root element of the XML document, parsed by lxml.etree
        server:   `str`
           pygcn server, defined in utils.gcn_server()
           options: local, eApps, Atlantic_2, Atlantic_3, Linode
        port:    `str`
           pygcn port, defined in utils.gcn_server()
           options: public, LVC, AMON

        Examples
        --------       
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.pygcn_serve(root='bayestar.xml',server='local',port='public')
        """        
        try:            
            from gcn.cmdline import serve_main
        except:
            self.logger.info ('Error: no pygcn installed')
            return
        
        kwargs = self.setkeys(kwargs)
        
        if kwargs['server']!='local':
            self.logger.info ('Warning: server set to send alert should be local')
            return
        if kwargs['root'] is None:
            self.logger.info ('Warning: root file which intend to send is not set')
            return
        
        server, port = self.gcn_server(kwargs['server'], kwargs['port'])
        if server is None: return
        
        sys.exit(serve_main(args=[kwargs['root'],'--host', '%s:%s'%(server,port)]))        
        
    def pygcn_monitor(self, func, server='local', port='public'):
        """serve root elements via pygcn

        Parameters
        ----------        
        func :        `function`
          function that was called once a gcn was received
        server:   `str`
           pygcn server, defined in utils.gcn_server()
           options: local, eApps, Atlantic_2, Atlantic_3, Linode
        port:    `str`
           pygcn port, defined in utils.gcn_server()
           options: public, LVC, AMON

        Examples
        --------       
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.pygcn_monitor(server='local',port='public')
        """ 
        try:
            import gcn
        except:
            self.logger.info ('Error: no pygcn installed')
            return
        
        kwargs = self.setkeys(kwargs)
        if kwargs is None: return
        
        server, port = self.gcn_server(kwargs['server'], kwargs['port'])
        if server is None: return
        
        gcn.listen(handler = func, host=server, port=port)            
    
    def parse_trigger(self, **kwargs):
        """parser hpmay via root, fits, url, xml, and coo approach in order.

        Parameters
        ----------
        url :        `string`         
          url of healpix fits
        hpfile :        `string`         
          file name or url of healpix fits
        vofile :        `string`
          file name or url of voenevt XML file
        root :        `class`
          root element of the XML document, parsed by lxml.etree
        coo :        `list`
          [ra, dec, err]: ra, dec is the center of trigger 
          while err is the estimated error box region

        Other Parameters
        ----------
        `savefits`, `wdir, `nest`, `clobber`

        Examples
        --------       
        >>> from kobe import trigger
        >>> a = trigger(url='https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.parse_trigger()
        >>> a.hpmap
        array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
               2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
        >>> a.hpd1
        array([144.84488626, 154.92359484, 140.96139722, ..., 272.89617891,
               221.54107632, 282.52651738])
        """        
        for _func in [self.root, self.fits, self.url, self.xml, self.coo]:                 
            if not self.hpmap is None:
                continue
            _func(**kwargs)            
            if not self.hpmap is None:
                self.logger.info ('parse trigger from %s'%_func)
            
    def url(self, *args, **kwargs):
        """parser hpmay via skymap url

        Parameters
        ----------
        url :        `string`         
          url of healpix fits

        Other Parameters
        ----------
        `savefits`, `wdir, `nest`, `clobber`

        Examples
        --------       
        >>> from kobe import trigger
        >>> a = trigger()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        # or 
        >>> a.url(url='https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.hpmap
        array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
               2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
        >>> a.hpd1
        array([144.84488626, 154.92359484, 140.96139722, ..., 272.89617891,
               221.54107632, 282.52651738])
        """       
        kwargs = self.setkeys(kwargs)
        
        if not self.hpmap is None:
            if kwargs['clobber']:
                self.hpmap = None
            else:
                self.logger.info ('Warning: hpmap already parsed')
                return
            
        # parse fits_url
        if len(args) > 0:
            url = args[0]
        else:
            url = kwargs['url']
        if not type(url) is str:return
        self.hpmap, self.hpd1, self.hpd2, self.hpd3, self.hphdr = self.read_skymap(url, **kwargs)
                    
        # savefits
        if not kwargs['savefits'] is None:
            self.download_skymap(url, kwargs['wdir'], kwargs['savefits'], kwargs['clobber'])
            
    def fits(self, *args, **kwargs):
        """parser hpmap via skymap file 
        Parameters
        ----------
        hpfile :        `string`         
          file name or url of healpix fits

        Other Parameters
        ----------
        `nest`, `clobber`
        
        Examples
        --------       
        >>> from kobe import triggers
        >>> a = triggers(nest=True)
        >>> a.fits('S190510g_bayestar.fits.gz')
        >>> a.hpmap
        array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
               2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
        """        
        kwargs = self.setkeys(kwargs)
        
        if not self.hpmap is None:
            if kwargs['clobber']:
                self.hpmap = None
            else:
                self.logger.info ('Warning: hpmap already parsed')
                return
            
        if len(args) > 0:
            skymap1 = args[0]            
        else:
            skymap1 = kwargs['hpfile']
        if not type(skymap1) is str:return

        if self.ckdir(**kwargs):
            skymap2 = '%s/%s' % (kwargs['wdir'], skymap1)
            
        for _skymap in [skymap1, skymap2]:
            if not self.hpmap is None:continue
            self.hpmap, self.hpd1, self.hpd2, self.hpd3, self.hphdr = \
                self.read_skymap(_skymap, **kwargs)            
                    
    def xml(self, *args, skymap_key='skymap_fits', obst_key='ISOTime',
            radec_keys=['C1','C2','Error2Radius'], **kwargs):       
        """parser hpmap via voenevt XML file

        Parameters
        ----------
        vofile :        `string`
          file name or url of voenevt XML file
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

        Other Parameters
        ----------
        `savefits`, `wdir`, `nside`, `coord`, `clobber`

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
        kwargs = self.setkeys(kwargs)
        
        if not self.hpmap is None:
            if kwargs['clobber']:
                self.hpmap = None
            else:
                self.logger.info ('Warning: hpmap already parsed')
                return
            
        if len(args) > 0:
            xml = args[0]     
        else:
            xml = kwargs['vofile']
        if not type(xml) is str:return

        if self.ckdir(**kwargs):
            xml2 = '%s/%s' % (kwargs['wdir'], xml)
                
        # read xml file, obtain root obj
        for _xml in [xml, xml2]:
            if not self.hpmap is None:continue                        
            root = self.parse_xml(_xml)
            # parse root
            if not root is None:
                self.root(root, skymap_key=skymap_key,
                          obst_key=obst_key, radec_keys=radec_keys, **kwargs)                
        
    def root(self, *args, skymap_key='skymap_fits', obst_key='ISOTime',
             radec_keys=['C1','C2','Error2Radius'], nside=512, **kwargs):        
        """parser hpmap via a root element

        Parameters
        ----------
        root :        `class`
          root element of the XML document, parsed by lxml.etree
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
        nside : int
           healpix resolution

        Other Parameters
        ----------
        `savefits`, `wdir`, `coord`,  `clobber`

        Examples
        --------              
        >>> from kobe import triggers
        >>> a = triggers(clobber=True, savefits=True)        
        # use pygcn to moniror ligo server
        >>> import gcn
        >>> def handler(payload, root): 
        >>>     a.root(root)
        >>>     b = a.calc_loc_area()
        >>>     print (b)
        >>> gcn.listen(handler=handler)
        """        
        kwargs = self.setkeys(kwargs)
        
        if not self.hpmap is None:
            if kwargs['clobber']:
                self.hpmap = None
            else:
                self.logger.info ('Warning: hpmap already parsed')
                return
        
        if len(args) > 0:
            root = args[0]
        else:
            root = kwargs['root']
            
        try:            
            self.voinf = {elem.attrib['name']: elem.attrib['value']
                          for elem in root.iterfind('.//Param')}   
        except:
            self.logger.info ('root object wrong format') 
            return

        # 1. find skymap url
        url = self.parse_root_skymap(root, skymap_key)        
        if not url is None:
            self.url(url, **kwargs)
        if not self.hpmap is None:
            return
                
        # 2. find ra,dec,error box        
        coo = self.parse_root_coo(root, nside, radec_keys, obst_key)
        if not coo is None:
            self.coo(coo, **kwargs)                         
        
    def coo(self, *args, nside=512, **kwargs):        
        """parse hpmap with ra, dec, error box
        
        Parameters
        ----------
        coo :        `list`
          [ra, dec, err]: ra, dec is the center of trigger 
          while err is the estimated error box region
        nside : int
           healpix resolution

        Other Parameters
        ----------
        `savefits`, `wdir`, `coord`, `nest`, `clobber`
       
        Examples
        --------                               
        >>> from kobe import trigger
        >>> a = trigger()       
        # generate a skymap with center at ra=30 deg, dec=-20 deg, with 10 deg sky area
        >>> a.coo([30, -20, 10])
        >>> a.hpmap
        array([8.92728175e-11, 8.92185286e-11, 8.91040539e-11, ...,
               6.12155463e-11, 6.08690152e-11, 6.10351951e-11])
        """        
        kwargs = self.setkeys(kwargs)
        
        if not self.hpmap is None:
            if kwargs['clobber']:
                self.hpmap = None
            else:
                self.logger.info ('Warning: hpmap already parsed')
                return

        if len(args) > 0:
            coo = args[0]
        else:
            coo = kwargs['coo']

        if len(coo) != 3: return
        _ra, _dec, _loc = float(coo[0]), float(coo[1]), float(coo[2])
        
        _radius = np.radians(_loc) # from deg to radians
        _theta, _phi = np.pi/2.-np.radians(_dec),np.radians(_ra)        
        _pmap = np.zeros(hp.nside2npix(nside))    
        _index = hp.ang2pix(nside, _theta, _phi, **self.getkeys(kwargs, 'pixtool'))
        _pmap[_index]+=1
        _pmap=hp.smoothing(_pmap,fwhm=_radius,sigma=None)        
        self.hpmap = _pmap/sum(_pmap)

        # savefits
        if not kwargs['savefits'] is None:
            self.savefits(**kwargs)            
                
    def savefits(self, **kwargs):
        """save self.hpmap to '%s/%s' % (kwargs['wdir'], filename)

        Parameters
        ----------        
        `savefits`, `wdir`, `coord`, `nest`         
        """
        kwargs = self.setkeys(kwargs)
        if self.ckdir(**kwargs):
            skymap = '%s/%s' % (kwargs['wdir'], kwargs['savefits'])            
            hp.write_map(skymap, self.hpmap, extra_header=self.hphdr, **kwargs)

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
        
        coo = []
        for _kk in radec_keys:
            _val = root.findtext('.//%s'%_kk)                      
            coo.append(_val)
            
        if len(coo) != 3: return
        
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
               ('RA',float(coo[0])),
               ('DEC',float(coo[1])),
               ('ERROR-BOX',float(coo[2])),
               ('NSIDE',nside),
               ('MJD-OBS',mjdobs),
               ('DATE-OBS',timeobs)]
        
        self.hphdr = dict(hdr)
        return coo
    
    def trigger_report(self, cls=[.5, .9], style='sms', keys_to_checkxml=('GraceID',
                'AlertType', 'Group', 'FAR', 'Terrestrial', 'HasNS', 'HasRemnant',
                'BNS', 'BBH', 'NSBH', 'Instruments', 'EventPage'),
                keys_to_checkfits = ('DISTMEAN', 'DISTSTD', 'DATE-OBS',
                'MJD-OBS', 'OBJECT', 'INSTRUME', 'CREATOR'), **kwargs):       
        """make a summary report for input source. 
        report is built from healpix fits header and XML file.
        specify keys in keys_to_checkxml, keys_to_checkfits

        Parameters
        ----------   
        style     : `str`
            when making report, specify the stype
            current options: `email`, `slack`, `sms`
            can be added as you want in kobe.KBtriggers.triggers.trigger_report()
        keys_to_checkxml  : `list`
           which parameters from voevent XML would be reported
        keys_to_checkfits : `list`
           which parameters from healpix fits header XML would be reported   
        cls       : `list` 
        confidence levels for contours,
        e.g. cls=[.5, .9] would construct 50% and 90% CL contours    

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
        >>> a.make_report(keys_to_checkfits=['DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS'])
        'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 '
        """        
        kwargs = self.setkeys(kwargs)
        
        _dict =  {}
        for _kk,_vv in zip([keys_to_checkxml, keys_to_checkfits],
                           [self.voinf, self.hphdr]):
            if _vv is None: continue
            for el in _kk:
                if el in _vv:                    
                    _dict[el] = _vv[el]
                    
        alist = self.calc_loc_area(**kwargs)
        if not alist is None:
            for el in alist:                
                _dict[el] = alist[el]
                
        _report = ''
        for el in _dict:
            try:
                _vv = '%.2f' % _dict[el]
            except:
                _vv = '%s' % _dict[el]                
            if style == 'sms':   _report += '%s:%s '%(el,_vv)
            elif style == 'slack':  _report += '`%s`:%s\n'%(el,_vv)
            elif style == 'email':  _report += '#\t%s:\t%s\n'%(el,_vv)
            
        if not _report in self.texts:
            self.texts += _report

    def sim_candidates(self, k, nest=False, **kwargs):

        assert type(k) is int
        assert self.hpmap is not None
        kwargs = self.setkeys(kwargs)
        
        _id = random.choices(range(len(self.hpmap)), weights=self.hpmap, k=k)
        theta, phi = hp.pix2ang(hp.npix2nside(len(self.hpmap)), _id, nest=nest)
        ra, dec = np.degrees(phi),np.degrees(np.pi/2.-theta)
        infos = self.calc_dis_inf(ra=ra, dec=dec, nest=nest, **kwargs)
        if infos is None: return
        dmu, dsigma, dnorm = infos
        _dist=[]
        for _dmu, _dsigma in zip(dmu, dsigma):
            _g0, _nn = False, 0
            while not _g0:
                _nn += 1
                _dd = np.random.normal(_dmu, _dsigma, 1)[0]
                if _dd > 0 or _nn > 100: _g0=True                
            _dist.append(_dd)            
        self.candidates = Table([range(len(ra)),ra,dec,_dist],
                                names=('n', 'ra', 'dec', 'dist'))
        
    def calc_loc_contours(self, cls=[.5,.9], **kwargs):        
        """calculate indices located in specific confidence level region of trigger map

        Parameters
        ----------   
        cls       : `list` 
        confidence levels for contours,
        e.g. cls=[.5, .9] would construct 50% and 90% CL contours    
        
        returns
        ----------
        indexlist :   `dictionary`
          dictionary of healpix indices corresponding input C.L.

        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_loc_contours()
        {0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.9: array([2414072, 2418168, 2416119, ...,  783552,  734398,  771264])}
        >>> a.calc_loc_contours(cls=[.1,.5,.99])
        {0.1: array([2414072, 2418168, 2416119, ..., 1570953, 1573001, 1573000]), 0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.99: array([2414072, 2418168, 2416119, ..., 1309038, 1309052, 1309051])}
        """
        
        kwargs = self.setkeys(kwargs)
        assert self.is_seq(cls)
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        
        indexlist = {}
        sortid = self.hpmap.argsort()[::-1]
        sorted_hpx = self.hpmap[sortid]
        cumulative_hpx = np.cumsum(sorted_hpx)
        for _cl in cls:
            if len(sorted_hpx[cumulative_hpx<_cl]) > 0:
                _limit = sorted_hpx[cumulative_hpx<_cl][-1]
                indexlist[_cl] = np.array([i for i in sortid if self.hpmap[i] >= _limit])
        return indexlist
    
    def calc_loc_area(self, cls=[.5,.9], **kwargs):       
        """calculate sky localization region area (unit in sq. deg) for different confidence level region of trigger map

        Parameters
        ----------   
        cls       : `list` 
        confidence levels for contours,
        e.g. cls=[.5, .9] would construct 50% and 90% CL contours    

        returns
        ----------
        arealist :   `dictionary`
          dictionary of area corresponding input C.L.

        Examples
        --------                               
        >>> from kobe import triggers
        >>> a = triggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_loc_area()
        {0.5: 575.5456172035578, 0.9: 3463.7648737864856}
        >>> a.calc_loc_area(cls=[.1,.5,.99])
        {0.1: 36.351906008208665, 0.5: 575.5456172035578, 0.99: 11508.39446313552}
        """
        kwargs = self.setkeys(kwargs)        
        _ilist = self.calc_loc_contours(cls=cls, **kwargs)
        if _ilist is None: return        
        
        nside = hp.get_nside(self.hpmap)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2
        _alist = {}
        for _cl in cls:
            _area = areapix*len(_ilist[_cl])
            _alist[_cl] = _area
        return _alist

    def calc_loc_prob_candidates(self, nest=False, add=False, **kwargs):
        assert 'candidates' in self.__dict__
        assert not self.candidates is None       
        _prob = self.calc_loc_prob(ra=self.candidates['ra'],
                    dec=self.candidates['dec'], nest=nest, **kwargs)
        if add:
            probs = Column(_prob, name='prob')                                
            self.candidates.add_column(probs)
        else:
            return _prob

    def calc_dis_sigma_candidates(self, add=False, nest=False, **kwargs):
        assert 'candidates' in self.__dict__
        assert not self.candidates is None        
        _dsigma = self.calc_dis_sigma(ra=self.candidates['ra'],
                                      dec=self.candidates['dec'],
                                      dist=self.candidates['dist'],
                                      nest=nest, **kwargs)
        if add:
            probs = Column(_dsigma, name='dsigma')                                
            self.candidates.add_column(probs)
        else:
            return _dsigma
    
    def calc_loc_prob_pointings(self, data=None, add=False, nest=False, **kwargs):

        if data is None:
            assert 'pointings' in self.__dict__
            _data = self.pointings.data
        else:
            _data = data
        assert not _data is None
        
        if 'fovra' in _data.keys():
            _prob = self.calc_loc_prob(ra=_data['ra'],
                        dec=_data['dec'], fovra=_data['fovra'],
                        fovdec=_data['fovdec'], nest=nest, **kwargs)        
        else:
            _prob = self.calc_loc_prob(ra=_data['ra'],
                        dec=_data['dec'], fovra=None, fovdec=None,
                        nest=nest, **kwargs)
        if add:
            probs = Column(_prob, name='prob')           
            self.pointings.data.add_column(probs)            
        else:
            return _prob
        
    def cut_pointings_hpmap(self, data=None, cls=None, nest=None, frac=0, **kwargs):                
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
        kwargs = self.setkeys(kwargs)        
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        nside = hp.get_nside(self.hpmap)                     
        areasingle =  hp.nside2pixarea(nside, degrees=True)
        idlist = self.calc_loc_contours(cls=cls, **kwargs)

        if data is None:
            assert 'pointings' in self.__dict__
            _data = self.pointings.data
        else:
            _data = data
        assert not _data is None
        
        _data0 = {} 
        for cc in idlist:
            if 'fovra' in _data:                        
                idtiles = self.ipix_in_box(_data['ra'], _data['dec'],
                                _data['fovra'],_data['fovdec'], nside, nest)                
                ele = self.common_ele(idtiles, idlist[cc])
                lmax = max([len(ll) for ll in idtiles])                
                nn = 0
                for ii in range(len(idtiles)):
                        _ele = ele[ii*lmax: (ii+1)*lmax-1]                        
                        _areasi = np.count_nonzero(_ele)*areasingle                        
                        _frac = _areasi/_data['fovra'][ii]/_data['fovra'][ii]
                        if _frac > frac: _data0[cc] = _data[ii]
            else:
                _theta, _phi = np.pi/2.-np.radians(_data['dec']),np.radians(_data['ra'])
                idgalaxies = hp.ang2pix(nside, _theta, _phi, nest=nest)
                ele = self.common_ele(idgalaxies, idlist[cc])                
                _data0[cc] = _data[[i for i in ele]]                            
        return _data0
            
    def calc_loc_prob(self, ra=None, dec=None, theta=None, phi=None,
                      fovra=None, fovdec=None, nest=False, **kwargs):
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

        Other Parameters
        ----------   
        `nest`

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
        >>> a.calc_loc_prob(ra=[1,10], dec=[1,-20])
        array([7.36631063e-07, 8.50340927e-06])
        # calculate total probs for vertice 1: ra=1, dec=1, fovra=1, fovdec=1 and vertice 2: ra=10, dec=-20, fovra=10, fovdec=10
        >>> a.calc_loc_prob(ra=[1,10], dec=[1,-20],fovra=[1,10],fovdec=[1,10])
        [5.8577539241927334e-05, 0.042806591674792185]
        """
         
        kwargs = self.setkeys(kwargs)        
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        nside = hp.get_nside(self.hpmap)

        radecs = self.radecs(ra=ra, dec=dec, theta=theta, phi=phi)                
        if radecs is None: return        
        
        if not fovra is None and not fovdec is None:            
            if self.is_seq(fovra):
                assert len(fovra) == len(radecs[0])
            else:                    
                fovra += np.zeros(len(radecs[0]))
            if self.is_seq(fovdec):
                assert len(fovdec) == len(radecs[1])
            else:                    
                fovdec += np.zeros(len(radecs[1]))                    
            probs = []
            for _ra,_dec,_fovw,_fovh in zip(radecs[0], radecs[1], fovra, fovdec):
                ipix_poly=(self.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,nest))
                _probs = self.hpmap[ipix_poly].sum()
                probs.append(_probs)                
        else:
            theta, phi = np.pi/2.-np.radians(radecs[1]),np.radians(radecs[0])
            probs = hp.get_interp_val(self.hpmap, theta, phi, nest=nest, lonlat=False)
        return probs

    def calc_loc_cl(self, ra=None, dec=None, theta=None, phi=None,
                    nest=False, **kwargs):       
         
        kwargs = self.setkeys(kwargs)        
        if self.hpmap is None:
            self.logger.info ('Warning: hpmap not parsed')
            return
        nside = hp.get_nside(self.hpmap)

        radecs = self.radecs(ra=ra, dec=dec, theta=theta, phi=phi)                
        if radecs is None: return        
        theta, phi = np.pi/2.-np.radians(radecs[1]),np.radians(radecs[0])
        idlist = hp.ang2pix(nside, theta, phi, nest=nest)
        
        sortid = self.hpmap.argsort()[::-1]        
        sorted_hpx = self.hpmap[sortid]        
        cumulative_hpx = np.cumsum(sorted_hpx)

        _idcl = []
        for _id in idlist:
            _id0 = np.where(sortid==_id)
            _idcl.append(cumulative_hpx[_id0][0])
        return _idcl

    def calc_dis_inf(self, ra=None, dec=None, theta=None,
                     phi=None, nest=False, **kwargs):
        """get trigger distance information, i.e. mean, variance and normalization, of specific direction
        for CBC GW only
                        
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
                    
        Other Parameters
        ------------         
        `nest`

        returns
        ----------
        mean, var, normalization :   `list`
          mean, variance, normalization lists at specified coordinates

        Examples
        --------                
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_dis_inf(ra=[1,10], dec=[1,-20])
        (array([233.86346031, 273.32868364]), array([121.25377326, 112.45704833]), array([1.44320589e-05, 1.14500239e-05]))
        """
                                
        kwargs = self.setkeys(kwargs)        
        if self.hpd1 is None or self.hpd2 is None or self.hpd3 is None:            
            self.logger.info ('Warning: no distance information')
            return
        nside = hp.get_nside(self.hpmap)
        radecs = self.radecs(ra=ra, dec=dec, theta=theta, phi=phi)                
        if radecs is None: return
        theta, phi = np.pi/2.-np.radians(radecs[1]),np.radians(radecs[0])
        ipix = hp.ang2pix(nside,theta,phi,nest=nest)        
        dmu, dsigma, dnorm = self.hpd1[ipix], self.hpd2[ipix], self.hpd3[ipix]        
        return dmu, dsigma, dnorm
    
    def calc_dis_prob(self, ra=None, dec=None, theta=None, phi=None,
                      nest=False, limmag=20, filt='r', dtemp=40, **kwargs):        
        """calculate trigger probability for specific directions, based on distance informations
        for CBC GW, the distance mean, variance and normalization at each direction is provided
        by assuming the limiting distance of one telescope on target, e.g. 200 Mpc on KNe,
        thus, triggers.calc_dis_prob will calculate the Gaussian probability from 0 to 200 Mpc,
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
            
        Other Parameters
        ------------         
        `nest`

        returns
        ----------
        problist :   `list`
          probabilities at specified coordinates

        Examples
        --------                
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_dis_prob(ra=[1,10], dec=[1,-20])
        [3.6100349724294025e-06, 3.089517808764708e-06]
        """
                                
        kwargs = self.setkeys(kwargs)
        if self.hpd1 is None or self.hpd2 is None or self.hpd3 is None:            
            self.logger.info ('Warning: no distance information')
            return
        nside = hp.get_nside(self.hpmap)
        pixarea = hp.nside2pixarea(nside, degrees=True)        
        infos = self.calc_dis_inf(ra=ra, dec=dec, theta=theta,
                                  phi=phi, nest=nest, **kwargs)
        if infos is None: return
        dmu, dsigma, dnorm = infos
        
        from scipy.stats import norm
        limdist = self.limdist(limmag=limmag, filt=filt, dtemp=dtemp)
        if limdist is None:return
        r = np.linspace(0, limdist)        
        probl = [dnorm * norm(dmu, dsigma).pdf(rr)/pixarea for rr in r]
        probs = [max(ii) for ii in list(map(list, zip(*probl)))]               
        return probs    

    def calc_dis_sigma(self, ra=None, dec=None, theta=None,
                       phi=None, dist=None, nest=False, **kwargs):
        """calculate at specific direction of CBC GW, if a distance is interesting or not
                        
        Parameters
        ------------         
        ra :   `int`, `float`-or sequence
          ra, unit in deg
        dec :  `int`, `float`-or sequence
          dec, unit in deg
        theta :  `int`, `float`-or sequence
          healpix theta, unit in radians
        phi :    `int`, `float`-or sequence
          healpix phi, unit in radians       
        dist :  `int`, `float`
          distance, unit in Mpc        
            
        Other Parameters
        ------------         
        `nest`

        returns
        ----------
        problist :   `list`
          probabilities at specified coordinates

        Examples
        --------                
        >>> from kobe import triggers
        >>> a = triggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_dis_sigma(ra=[1,10], dec=[1,-20])
        [3.6100349724294025e-06, 3.089517808764708e-06]
        """
        if dist is None: return        
        kwargs = self.setkeys(kwargs)
        if self.hpd1 is None or self.hpd2 is None or self.hpd3 is None:            
            self.logger.info ('Warning: no distance information')
            return
        nside = hp.get_nside(self.hpmap)        
        infos = self.calc_dis_inf(ra=ra, dec=dec, theta=theta,
                                  phi=phi, nest=nest, **kwargs)
        if infos is None: return
        dmu, dsigma, dnorm = infos        
        if self.is_seq(dist):
            assert len(dist) == len(dmu)
            return abs(dist - dmu)/dsigma
        else:            
            return abs(dist - dmu)/dsigma
    
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
        try:
            tree = etree.parse(xmlfile)
            root = tree.getroot()        
            return root
        except:
            self.logger.info ('Error: failed to parse trigger from xml')
            return
            
    def read_skymap(self, _fits, **kwargs):
        """read healpix fits or fits url, obtain map and header
        """
        kwargs = self.getkeys(kwargs, 'readmap')        
        try:
            # read 3D trigger healpix map (for GW)            
            tmap, header = hp.read_map(_fits,field=[0, 1, 2, 3], h=True, **kwargs)
            return tmap[0], tmap[1], tmap[2], tmap[3], dict(header)
        except:
            # if failed, read 2D map            
            try:
                tmap, header = hp.read_map(_fits, h=True, **kwargs)
                return tmap, None, None, None, dict(header)
            except:
                return None, None, None, None, None            
                    
    def download_skymap(self, url, wdir, fitsname, clobber):            
        """look up URL of sky map, download sky map, and parse FITS file    
        """
        if fitsname is None: return        
        skymap = '%s/%s' % (wdir, fitsname)
        if os.path.exists(skymap) and not clobber:
            self.logger.info ('Warning: %s exists already'%skymap)
        try: 
            import wget
            self.logger.info ('download skymap via wget')
            wget.download(url, out=skymap)            
        except ImportError as e:
            try:
                import requests,tempfile,shutil
                self.logger.info ('download skymap via requests')
                
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
            except ImportError as e:
                self.logger.info (e)                

    def generatep_trigger(self, num, limra=[0,360.], limdec=[-89,89], fovra=1.,
            fovdec=1., dfov=[5.,10.,20.], looplim=3, nest=False, **kwargs):
        """monte carlo approach on shiftra and shiftdec to maxmize 
        trigger probability with a number of pointings
                
        Parameters
        ----------        
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
        shiftra :       `float`                  
                ra of initial tiling will be 0 + shiftra
        shiftdec :      `float`           
                dec of initial tiling will be 0 + shiftdec
        clobber   : `bool`  
            if any products already available in self, clobber=True will redo it             
        """
        assert type(num) is int
        assert 'pointings' in self.__dict__              
        kwargs = self.setkeys(kwargs)
        if not self.pointings.data is None and not kwargs['clobber']:
            self.logger.info ('Warning: tiling data already parsed')
            return
                
        assert not self.hpmap is None, 'parse healpix map first'       
                            
        # monte carlo for tiling
        _log = [0.]
        
        for nn in dfov:
                        
            self.logger.info ('- searching in fovh/%i fovw/%i'%(nn,nn))
            shiftra, shiftdec = fovra/nn, fovdec/nn

            ntrial, answ1 = 0, False
            while not answ1:  # if angle OK
                angle = random.uniform(0,2*np.pi)
                self.logger.info ('\t %i with angle: %.2f'%(ntrial,angle))
                _shiftra = np.sqrt(shiftra**2+shiftdec**2)*np.sin(angle)
                _shiftdec = np.sqrt(shiftra**2+shiftdec**2)*np.cos(angle)
                
                answ2 = False
                while not answ2:  # if OK, go on, if no, change
                                                                                
                    # generate pointings
                    shiftra += _shiftra
                    shiftdec += _shiftdec
                    _data = self.pointings.generatep(limra=limra, limdec=limdec, fovra=fovra,
                                    fovdec=fovdec, shiftra=shiftra, shiftdec=shiftdec,
                                    returndata=True, **kwargs)
                    
                    # cal prob for tiling list
                    t = self.calc_loc_prob_pointings(nest=nest, data=_data)
                    t = sorted(t)[::-1][:num]
                                        
                    # judge converge or not
                    if sum(t)>_log[-1]:  # accept, direction correct
                        _log.append(sum(t))
                        self.logger.info ('\tcovered %.5e probs'%_log[-1])                        
                    else:  # reject, change direction
                        ntrial += 1
                        answ2 = True
                    if ntrial >= looplim: answ1=True
        return _data

    def rankp_trigger(self, mode=1, nside=64, nest=False, threshold=1,sort=1):
        """rank pointings (greedy approach): rank from high trigger probability to low

        Parameters
        ----------          
        mode     : `int`
                  1. strict approach: follow above strategy strictly
                  2. adjacent approach: start from the westest point or highest probability point and arange the next pointing either adjacent or closest to the previous one
        """                
        assert 'pointings' in self.__dict__                     
        assert not self.pointings.data is None, 'Warning: tiling data not parsed'
        _data = self.pointings.data        
        
        # define start point        
        assert not self.hpmap is None, 'parse healpix map first'
        t = self.calc_loc_prob_pointings(nest=nest)                    
        _idstart = np.argmax(t)
        
        if mode == 1:
            ipx = np.argsort(_data['ra'])
            
        elif mode == 2:            
            ipx = self.adjacent_move(_data['ra'],_data['dec'],_idstart,threshold=1,sort=1)

        else:
            self.logger.info ('wrong mode')
            return
            
        self.pointings.data = _data[ipx]
        self.pointings.data['n'] = np.arange(len(self.pointings.data['n']), dtype=int)
