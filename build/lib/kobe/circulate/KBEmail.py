#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/circulate/circulate.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
import numpy as np
import healpy as hp
import astropy.time

__all__ = ('PstCirculateInfos')

class PstCirculateInfos():
    """
    * implement the following functions:
    -> __init__(self)    
    """

    # Static version info
    version = 1.0
        
    def __init__(self, hpmap=None, fithdr=None,
                 xmlinf=None, defconf=None, logger = None):
        """
        Read trigger map from a given source
        only one needed, if more than 1, read the latest        

        Parameters:
        -----------
        hpmap       `list`           healpix fits for 3d localization
        fithdr      `dict`           healpix fits map informations
        xmlinf      `dict`           voevent xml file informations
        defconf     `dict`           default configure
        logger      `class`          logger
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

        self.conf = {
            'wdir':      '/tmp/',   # `str`   working directory
            'savefits':  None,      # `str`   save file or not
            'nside':     512,       # `int`   healpix map default resolution
            'coord':     'C',       # `str`   healpix map defualt coordinate system: G, C, E
            'nest':      False,     # `bool`  healpix map defualt ordering: nest or ring
            'cls':       [.5,.9],   # `list`  confidence levels for contours
            'style':     'sms'      # `str`   report type
        }
        
        if defconf is None: return
        for _k in self.conf.keys():
            if _k in defconf.keys():
                self.conf[_k] = defconf[_k]
            else:
                self.logger.info ('### Warning: use default value for %s'%_k)

    '''
    parse sources, to obtain self.data
    Source options:
    -----------
        skymap_url  `string`         url of healpix fits
        skymap      `string`         filename of healpix fits
        xmlfile     `string`         the filename of XML
        root        `class`          root element of the XML document parsed by lxml.etree
        coolist     `list`           [ra,dec,err], `float`, `float`, `float`   
    '''
    
    def texts(self, url, wdir=None, savefits=None):
        """ 
        apply parser on skymap url        
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        
        # parse fits_url
        if not url is None: self.fits(url)

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
            
    def attachments(self, skymap, nest=None):
        """ 
        apply parser on healpix fits or fits url
        """
        
        _hp = self.checkhp()
        if _hp: return
        if not nest is None: self.conf['nest'] = nest

        if not skymap is None:            
            hpmap, fitsinfos = self.read_skymap(skymap, self.conf['nest'])
            self.data['hpmap'] = hpmap            
            self.data['fithdr'] = fitsinfos

    def slack(self, xml, wdir=None, savefits=None, nside=None, coord=None):
        """ 
        apply parser on voenevt XML file
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
            
    def email(self, root, wdir=None, savefits=None, nside=None, coord=None):
        """ 
        apply parser on root element of the XML document, parsed by lxml.etree
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
            
    def sms(self, coo, wdir=None, savefits=None,
            nside=None, coord=None, nest=None):
        """ 
        apply parser on ra, dec, error box
        obtain hpmap
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
                                  
    def esop2(self):
        if not self.data['hpmap'] is None:
            self.logger.info ('hpmap has already been parsed')
            return True
        else:
            return False

    def ssh(self):
        if not self.data['hpmap'] is None:
            self.logger.info ('hpmap has already been parsed')
            return True
        else:
            return False
