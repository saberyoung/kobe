#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/PstParseTriggers.py
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
        
__all__ = ('PstParseTriggers')

class PstParseTriggers():
    """PstParseTriggers: Read and process trigger healpix map from a given source

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
      default configure, if any PstParseTriggers parameter was included in defconf dictionary, 
      then its default value would be applied
    logger :     `class`          
      logging object

    See Also
    --------
    PstGetTilings, PstGetGalaxies

    Examples
    --------
    see also https://github.com/saberyoung/pstools/blob/master/notebook/test_triggers.ipynb

    >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
    >>> a = PstParseTriggers()
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
        
    def __init__(self, hpmap=None, fithdr=None,
                 xmlinf=None, defconf=None, logger=None):
        """
        initialize
        """

        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # ----- set trigger properties ----- #
        self.data               =   {}
        self.data['hpmap']      =   hpmap
        self.data['xmlinf']     =   xmlinf
        self.data['fithdr']     =   fithdr

        # ----- define keys to check xml ----- #
        '''
        skymap_key  `string`   key for skymap url, in XML   
        radec_keys  `list`     key for coordinate and uncertainty, in XML
        obst_key    `string`   key for timeobs, in XML
        '''
        self.obst_key           =   'ISOTime'
        self.skymap_key         =   'skymap_fits'
        self.radec_keys         =   ['C1','C2','Error2Radius']
        
        # ----- define informations for report ----- #
        self.keys_to_checkxml = ('GraceID', 'AlertType', 'Group', 'FAR',
                                 'Terrestrial', 'HasNS', 'HasRemnant',
                                 'BNS', 'BBH', 'NSBH', 'Instruments',
                                 'EventPage')
        self.keys_to_checkfits = ('DISTMEAN', 'DISTSTD', 'DATE-OBS',
                                  'MJD-OBS', 'OBJECT', 'INSTRUME',
                                  'CREATOR')

        # ----- define default parameters ----- #
        self.run_config(defconf)

    def run_config(self, defconf):
        """read default value for parameters

        Parameters
        ----------
        defconf :      `dict`
          input parameter settings, if None, use default

        Examples
        --------    
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> dict = {'cls': [.5,.9,.99]}
        >>> a.run_config(dict)
        """
        
        self.conf = {
            'wdir':      '/tmp/',   # `str`   working directory
            'savefits':  None,      # `str`   save file or not
            'nside':     512,       # `int`   healpix map default resolution
            'coord':     'C',       # `str`   healpix map defualt coordinate system: G, C, E
            'nest':      False,     # `bool`  healpix map defualt ordering: nest or ring
            'cls':       [.5,.9],   # `list`  confidence levels for contours
            'style':     'sms',     # `str`   report type           
            'lat':       24.625,    # `float`   la palma chile
            'lon':       70.403,    # `float`   la palma chile
            'alt':       2635,      # `float`   la palma chile
            'obstime':   None
        }
        
        if defconf is None: return
        for _k in self.conf.keys():
            if _k in defconf.keys():
                self.conf[_k] = defconf[_k]
            else:
                self.logger.info ('### Warning: use default value for %s'%_k)
    
    def url(self, url, wdir=None, savefits=None, nest=None):
        """parser skymap url, and build PstParseTriggers.data  

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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        
        # parse fits_url
        if not url is None: self.fits(url,nest=nest)

        # savefits
        if not self.conf['savefits'] is None:
            flag = self.download_skymap(self.conf['wdir'], self.conf['savefits'], url)
            if flag == 1:
                self.logger.info ('download skymap via wget')
            elif flag == 2:
                self.logger.info ('download skymap via requests')
            elif flag == 3:
                self.logger.info ('### Warning: failed to download skymap,'+\
                                  'install wget or requests first')
            
    def fits(self, skymap, nest=None):
        """parser skymap file, and build PstParseTriggers.data  

        Parameters
        ----------
        skymap :        `string`         
          file name or url of healpix fits
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead
        
        Examples
        --------       
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> a.url('S190510g_bayestar.fits.gz', nest=True)
        """
        _hp = self.checkhp()
        if _hp: return
        if not nest is None: self.conf['nest'] = nest

        if not skymap is None:            
            hpmap, fitsinfos = self.read_skymap(skymap, self.conf['nest'])
            self.data['hpmap'] = hpmap            
            self.data['fithdr'] = fitsinfos

    def xml(self, xml, wdir=None, savefits=None, nside=None, coord=None):       
        """parser voenevt XML file, and build PstParseTriggers.data  

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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> a.xml('LVC#S190510g-2-Initial.xml')
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        if not nside is None: self.conf['nside'] = nside
        if not coord is None: self.conf['coord'] = coord
        
        # read xml file, obtain root obj        
        root = self.parse_xml(xml)

        # parse root
        self.root(root, wdir=self.conf['wdir'],
                  savefits=self.conf['savefits'],
                  nside=self.conf['nside'],
                  coord=self.conf['coord'] )
            
    def root(self, root, wdir=None, savefits=None, nside=None, coord=None):        
        """parser root element, and build PstParseTriggers.data  

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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()
        >>> def handler(payload, root): 
        >>>     a.root(root)
        >>>     b = a.calc_area()
        >>>     print (b)
        >>> gcn.listen(handler=handler)
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        if not nside is None: self.conf['nside'] = nside
        if not coord is None: self.conf['coord'] = coord
        
        _hp = self.checkhp()
        if _hp: return
        
        if not root is None:            
            self.data['xmlinf'] = {elem.attrib['name']:
                                   elem.attrib['value']
                                   for elem in root.iterfind('.//Param')}   
            
            # 1. find ra,dec,error box: for AMON, GRB, etc   
            coo = self.parse_root_coo(root, self.conf['nside'])
            if not coo is None:
                self.coo(coo, wdir=self.conf['wdir'],
                         savefits=self.conf['savefits'],
                         nside=self.conf['nside'],
                         coord=self.conf['coord'])

            # 2. find skymap url: for GW
            url = self.parse_root_skymap(root)
            if not url is None:
                self.url(url, wdir=self.conf['wdir'],
                         savefits=self.conf['savefits'])
            
    def coo(self, coo, wdir=None, savefits=None,
            nside=None, coord=None, nest=None):        
        """use ra, dec, error box to build PstParseTriggers.data  
        
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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()       
        >>> a.coo([30, -20, 10])
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        if not nside is None: self.conf['nside'] = nside
        if not coord is None: self.conf['coord'] = coord
        if not nest is None: self.conf['nest'] = nest
        
        _hp = self.checkhp()
        if _hp: return

        _ra,_dec,_loc = coo
        _radius = np.radians(_loc) # from deg to radians
        _theta, _phi = np.pi/2.-np.radians(_dec),np.radians(_ra)        
        _pmap = np.zeros(hp.nside2npix(self.conf['nside']))
        _index = hp.ang2pix(self.conf['nside'], _theta, _phi,
                            nest=self.conf['nest'])
        _pmap[_index]+=1
        _pmap=hp.smoothing(_pmap,fwhm=_radius,sigma=None)                
        self.data['hpmap'] = _pmap/sum(_pmap)

        # savefits
        if not self.conf['savefits'] is None:
            flag = self.make_hpmap(self.conf['wdir'],
                                   self.conf['savefits'],
                                   self.conf['nest'],
                                   self.conf['coord'])
            if flag:
                self.logger.info ('make healpix map')            
            else:
                self.logger.info ('### Warning: failed to make healpix map')
                                  
    def checkhp(self):
        """check if PstParseTriggers is parsed from one source or not  
        
        returns
        ----------
        res :        `bool`
        """
                
        if not self.data['hpmap'] is None:
            self.logger.info ('hpmap has already been parsed')
            return True
        else:
            return False
        
    def make_hpmap(self, wdir, tempfile, nest, coord):
        """generate and store healpix map, from PstParseTriggers.data['hpmap']

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

        skymap = '%s/%s' % (wdir, tempfile)
        try:
            hp.write_map(skymap, self.data['hpmap'], nest=nest, coord=coord,
                         extra_header=self.data['fithdr'], overwrite=True)
            return True
        except:
            return False

    def parse_root_skymap(self, root):       
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
        if self.skymap_key in xmlinfos:
            self.logger.info ('obtain hpmap via %s' % self.skymap_key)
            skymap_url = xmlinfos[self.skymap_key]
        else:
            self.logger.info ('### Warning: no %s found in voevent'%self.skymap_key)
            skymap_url = None
        return skymap_url
            
    def parse_root_coo(self, root, nside):        
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
        
        _coolist,_coo = [],0
        for _kk in self.radec_keys:
            _val = root.findtext('.//%s'%_kk)
            if not _val is None:
                try:
                    _coolist.append(float(_val))
                    _coo += 1
                except:
                    self.logger.info ('### Warning: %s not float, skip'%_kk)
            else:
                self.logger.info ('### Warning: no %s found in voevent'%_kk)
                
        if _coo == 3:
            _ra,_dec,_loc = _coolist
        else:
            return None
        
        """
        make header            
        """        
        # obs time
        timeobs = root.findtext('.//%s'%self.obst_key)

        # obs mjd        
        mjdobs = astropy.time.Time(timeobs, format='isot', scale='utc').mjd

        # read classifier
        try:  tid = root.attrib['ivorn']
        except:  tid = None

        # header
        hdr = [('CREATOR','PSTOOLS'),
               ('IVORN',tid),
               ('RA',_ra),
               ('DEC',_dec),
               ('ERROR-BOX',_loc),
               ('NSIDE',nside),
               ('MJD-OBS',mjdobs),
               ('DATE-OBS',timeobs)]
        
        self.data['fithdr'] = dict(hdr)
        return _coolist
    
    def make_report(self, style=None):       
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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.make_report()
        'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 OBJECT:G331903 INSTRUME:H1,L1,V1 CREATOR:BAYESTAR '
        >>> a.keys_to_checkfits
        ('DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS', 'OBJECT', 'INSTRUME', 'CREATOR')
        >>> a.keys_to_checkfits=('DISTMEAN', 'DISTSTD', 'DATE-OBS', 'MJD-OBS')
        >>> a.make_report()
        'DISTMEAN:268.8566049372629 DISTSTD:108.0709050006497 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12476032978 '
        """
         
        if not style is None: self.conf['style'] = style
        
        _dict =  {}
        for el in self.keys_to_checkxml:
            try:
                _dict[el] = self.data['xmlinf'][el]
            except:
                self.logger.info ('### Warning: no %s found in voevent'%el)
        for el in self.keys_to_checkfits:
            try:
                _dict[el] = self.data['fithdr'][el]
            except:
                self.logger.info ('### Warning: no %s found in fits header'%el)

        _report = ''
        for el in _dict:
            if self.conf['style'] == 'sms':
                _report += '%s:%s '%(el,_dict[el])
            elif self.conf['style'] == 'slack':
                _report += '`%s`:%s\n'%(el,_dict[el])
            elif self.conf['style'] == 'email':
                _report += '#\t%s:\t%s\n'%(el,_dict[el])           

        if not self.conf['style'] in ['sms', 'email', 'slack']:
            self.logger.info ('### Warning: style %s not available'%self.conf['style'])
            self.logger.info ('### select from email, slack, sms')
        return _report

    def calc_contours(self, cls=None):        
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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_contours()
        {0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.9: array([2414072, 2418168, 2416119, ...,  783552,  734398,  771264])}
        >>> a.calc_contours(cls=[.1,.5,.99])
        {0.1: array([2414072, 2418168, 2416119, ..., 1570953, 1573001, 1573000]), 0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.99: array([2414072, 2418168, 2416119, ..., 1309038, 1309052, 1309051])}
        """
                
        if not cls is None: self.conf['cls'] = cls
        
        _hp = self.checkhp()
        if not _hp:
            self.logger.info ('### Warning: no hpmap found')
            return

        from pst.cookbook import is_seq, is_seq_of_seq
        if is_seq_of_seq(self.data['hpmap']):
            (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        elif is_seq(self.data['hpmap']):
            hpx = self.data['hpmap']
        else:
            return
        
        indexlist = {}
        sortid = hpx.argsort()[::-1]
        sorted_hpx = hpx[sortid]
        cumulative_hpx = np.cumsum(sorted_hpx)
        for _cl in self.conf['cls']:
            if len(sorted_hpx[cumulative_hpx<_cl]) > 0:
                _limit = sorted_hpx[cumulative_hpx<_cl][-1]
                indexlist[_cl] = np.array([i for i in sortid if hpx[i] >= _limit])
        return indexlist
    
    def calc_area(self, cls=None):       
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
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_area()
        {0.5: 575.5456172035578, 0.9: 3463.7648737864856}
        >>> a.calc_area([.1,.5,.99])
        {0.1: 36.351906008208665, 0.5: 575.5456172035578, 0.99: 11508.39446313552}
        """
        
        if not cls is None: self.conf['cls'] = cls
        
        _ilist = self.calc_contours(cls=self.conf['cls'])
        if _ilist is None: return
        
        from pst.cookbook import is_seq, is_seq_of_seq
        if is_seq_of_seq(self.data['hpmap']):
            (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        elif is_seq(self.data['hpmap']):
            hpx = self.data['hpmap']
        else:
            return

        nside = hp.get_nside(hpx)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2
        _alist = {}
        for _cl in self.conf['cls']:
            _area = areapix*len(_ilist[_cl])
            _alist[_cl] = _area
        return _alist    

    def altaz(self, obstime=None, cls=None, nest=None,
              lat=None, lon=None, alt=None):
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
        cls :         `list`
          list of confidence level, default: [.5, .9]
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        returns
        ----------
        altaz :         `list`
          altaz list

        obstime :      `astropy.time`
          observing time, 

        observaroty :  `astropy.coordinate`
          observaroty location
        """                    
        
        _hp = self.checkhp()
        if not _hp:
            self.logger.info ('### Warning: no hpmap found')
            return
        
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
        from pst.cookbook import is_seq, is_seq_of_seq
        if is_seq_of_seq(self.data['hpmap']):
            (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        elif is_seq(self.data['hpmap']):
            hpx = self.data['hpmap']
        else:
            return
        nside = hp.get_nside(hpx)        
        if not nest is None: self.conf['nest'] = nest
        if not cls is None: self.conf['cls'] = cls
        _ilist = self.calc_contours(cls=self.conf['cls'])
        altazlist = {}
        for _cc in _ilist:
            _ii = _ilist[_cc]
            theta, phi = hp.pix2ang(nside, _ii, nest=self.conf['nest'], lonlat=False)        
            dec, ra = -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)        
            radecs = astropy.coordinates.SkyCoord(ra=ra*u.deg, dec=dec*u.deg)               

            # Alt/az reference frame at observatory, now
            frame = astropy.coordinates.AltAz(obstime=obstime, location=observatory)

            # Transform grid to alt/az coordinates at observatory, now
            altazlist[_cc] = radecs.transform_to(frame)            
        return altazlist, obstime, observatory

    def calc_vis(self, limalt=2.5, limsun=-18, obstime=None,
                 cls=None, nest=None, lat=None, lon=None, alt=None):
        """calculate available area, considering airmass and sun constrain

        Parameters
        ----------   
        limalt :      `float`
          limitation on airmass, e.g. 2.5, will remain only fields with airmass less than 2.5
        limsun :      `float`
          limitation on sun height, e.g. -18, will remain only fields with sun lower than -18
        lon :         `float` or None
          longtitude of telescope, defaul: la palma chile
        lat :         `float` or None
          latitude of telescope, defaul: la palma chile
        alt :         `float` or None
          altitude of telescope, defaul: la palma chile
        obstime :     `float` or None or astropy.time obj
          observing time.
          None for now; `float` number (unit in seconds) represents how much time before (negative) or after (positive) now; `astropy.time` set specific observing time
        cls :         `list`
          list of confidence level, default: [.5, .9]
        nest :          `bool`
          healpix ordering options: 
          if True, healpix map use `nest` ordering, otherwise, use `ring` instead

        returns
        ----------
        arealist :   `dictionary`
          dictionary of area corresponding input C.L.

        Examples
        --------                               
        >>> from pst.pipeline.PstParseTriggers import PstParseTriggers
        >>> a = PstParseTriggers()       
        >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
        >>> a.calc_vis(obstime='2020-01-06 00:34:52')
        {0.5: 198.07330026983539, 0.9: 1359.6557052420903}
        >>> a.calc_vis(obstime='2020-01-06 10:34:52')
        {0.5: 0.0, 0.9: 0.0}
        >>> a.calc_vis(obstime='2020-01-06 10:34:52',limsun=None)
        {0.5: 326.4852279871439, 0.9: 1809.0843382894639}
        """
        
        _hp = self.checkhp()
        if not _hp:
            self.logger.info ('### Warning: no hpmap found')
            return
        
        from pst.cookbook import is_seq, is_seq_of_seq
        if is_seq_of_seq(self.data['hpmap']):
            (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        elif is_seq(self.data['hpmap']):
            hpx = self.data['hpmap']
        else:
            return

        _alist = {}
        
        if limalt is None:
            self.logger.info ('no airmass constrain')
            return       
             
        nside = hp.get_nside(hpx)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2        
        altazl, obstime, observatory = self.altaz(obstime=obstime,
                        cls=cls, nest=nest, lat=lat, lon=lon, alt=alt)               
        for cc in altazl:
            altaz = altazl[cc]
            sun_altaz = astropy.coordinates.get_sun(obstime).transform_to(altaz)            
            if limsun is None:
                self.logger.info ('no sun constrain')
                cond = np.logical_and(altaz.secz <= limalt, altaz.secz >= 1)
            else:               
                cond = np.logical_and(np.logical_and(altaz.secz <= limalt, altaz.secz >= 1),
                                      sun_altaz.alt <= limsun*u.deg)           
            _alist[cc] = len(altaz[cond])*areapix
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

    def data_list(self):
        """  get PstParseTriggers.data

        returns
        ----------
        data :        `dictionary`
        """
        return self.data
            
    def set(self, _k, _v):
        """  set PstParseTriggers.data

        Parameters
        ----------
        _k :        `string`
          key
        _v :        `string` or `float` or `int` or etc
          value
        """
        self.data[_k] = _v
    
    def get(self, _k):
        """  get PstParseTriggers.data

        Parameters
        ----------
        _k :        `string`
          key

        returns
        ----------
        value :      `list`
        """
        return self.data[_k]

    def resetall(self):
        """  reset all of PstParseTriggers.data
        """
        self.__init__()

    def reset(self, _k):
        """  reset part of PstParseTriggers.data
        
        Parameters
        ----------
        _k :        `string`
          key
        """
        self.data[_k] = None
        
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
               
        try: # read 3D trigger healpix map (for GW)            
            tmap, header = hp.read_map(_fits,field=[0, 1, 2, 3],
                                       nest=nest,h=True)
        except: # if failed, read 2D map            
            try:
                tmap, header = hp.read_map(_fits, nest=nest,h=True)                
            except:
                tmap, header = None, None
        return tmap, dict(header)
    
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
    def download_skymap(wdir, tmpfile, url):            
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
        try: 
            import wget            
            if os.path.exists(skymap): os.remove(skymap)
            wget.download(url, out=skymap)
            
            flag = 1
        except:
            try:
                import requests,tempfile,shutil                

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
                    
                    flag = 2
            except:
                flag = 3
        return flag
