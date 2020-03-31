#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : kobe/KBcandidate.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import numpy as np
from astropy.table import Table
import requests
import logging
from kobe import vcandidates, circulate
import time, sys, fastavro, tarfile                
                
__all__ = ['candidates']

class candidates(vcandidates):
        """candidates: read a list of astronomical sources as candidates
        """
        
        candidates = None
        '''
        candidates data
        '''        

        lc = None
        '''
        target lightcurve
        '''
        
        def __init__(self, logger=None, **kwargs):
                """
                initialize candidate object
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set keys ----- #        
                self.defkwargs = self.setkeys(kwargs)

        def readc_ampel(self, alerts):
                
                _data = {'name':[], 'n':[], 'ra':[],
                         'dec':[], 'time':[], 'mag':[]}
                
                run_start = time.time()
                iter_count = 0
                                                          
                # Iterate over alerts
                for alert_content in alerts:
                                                         
                        _tranid = alert_content.tran_id                        
                        _alert = alert_content['candidate']
                                 
                        self.logger.info('accept alert %s'%_tranid)
                        _data['name'].append(_tranid)
                        _data['ra'].append(_alert['ra'])
                        _data['dec'].append(_alert['dec'])
                        _data['mag'].append(_alert['magpsf'])
                        _data['time'].append(_alert['jd'])
                                                
                self.logger.info(
                        "%i alert(s) processed (time required: %is)" % 
                        (iter_count, int(time.time() - run_start))
                )

                self.candidates = Table([np.arange(len(_data['ra'])),
                        _data['ra'], _data['dec'], _data['mag'], _data['time']],
                        names=('n', 'ra', 'dec', 'mag', 'time'))
        
        def readc_avro(self, tar_file, tar_mode="r:gz", rbt=.3, iter_max=5000):
                """
                For each alert: load, filter, ingest.
                """
                _data = {'name':[], 'n':[], 'ra':[],
                         'dec':[], 'time':[], 'mag':[]}
                
                run_start = time.time()
                iter_count = 0
                                                   
                iterable = tarfile.open(tar_file, mode=tar_mode)                
                
                # Iterate over alerts
                for content in iterable:
                                                               
                        reader = fastavro.reader(iterable.extractfile(content))
                        alert_content = next(reader, None)
                                
                        if alert_content is None:
                                self.logger.info("Reached end of tar files")
                                iterable.close()
                                break
                                
                        _tranid = alert_content['objectId']                       
                        _alert = alert_content['candidate']
                                                             
                        if _alert is None:
                                break

                        if _alert['rb'] > rbt:
                                self.logger.info('accept alert %s'%_tranid)
                                _data['name'].append(_tranid)
                                _data['ra'].append(_alert['ra'])
                                _data['dec'].append(_alert['dec'])
                                _data['mag'].append(_alert['magpsf'])
                                _data['time'].append(_alert['jd'])
                        else:
                                self.logger.info('reject alert %s'%_tranid)
                                
                        iter_count += 1
                        if iter_count == iter_max:
                                self._logger.info("Reached max number of iterations")
                                break
                        
                self.logger.info(
                        "%i alert(s) processed (time required: %is)" % 
                        (iter_count, int(time.time() - run_start))
                )

                self.candidates = Table([np.arange(len(_data['ra'])),
                        _data['ra'], _data['dec'], _data['mag'], _data['time']],
                        names=('n', 'ra', 'dec', 'mag', 'time'))
                        
        def readc(self, ra=None, dec=None, theta=None,
                  phi=None, dist=None, mag=None, time=None, hms=False, **kwargs):
                
                radecs = self.radecs(ra=ra, dec=dec, theta=theta, phi=phi, hms=hms)                
                if radecs is None: return                
                _data = {}
                for _i, _j in zip([dist,mag,time],['dist','mag','time']):
                        if _i is None:
                                _data[_j] = np.zeros(len(radecs[0]))
                                continue
                        assert type(_i) is type(radecs[0])
                        if self.is_seq(_i):
                                assert len(_i) == len(radecs[0])
                        else:
                                _i = [_i]
                        _data[_j] = np.array(_i)
                
                self.candidates = Table([np.arange(len(radecs[0])),
                                radecs[0], radecs[1], _data['dist'],
                                _data['mag'], _data['time']],
                                names=('n', 'ra', 'dec', 'dist', 'mag', 'time'))
                        
        def readc_file(self, filename, filetype='npz', split=' ',
                keys=['n','ra','dec','dist','mag'], **kwargs):
                
                _r = self.readlist(split=split, filename=filename,
                                   filetype=filetype, keys=keys, **kwargs)
                if _r is None: return
               
                _data = []
                for _key in keys:
                        _data.append(np.array(_r[_key]))
                self.candidates = Table(_data, names=keys)    

        def generatelc(self, db='https://api.astrocats.space/', tname='AT2017gfo',
                       keys = ['time','magnitude','band','upperlimit',
                               'e_time','e_lower_time','e_upper_time',
                               'e_magnitude','e_lower_magnitude','e_upper_magnitude'],
                       argkeys=['time', 'magnitude', 'band'], timeout=30,
                       keysfloat = ['time', 'magnitude'], timezero=None):
                """query lightcurve via OACAPI, https://github.com/astrocatalogs/OACAPI    

                Parameters
                -----------------
                db: str
                    https://api.astrocats.space/
                    https://api.sne.space/
                    https://api.tde.space/
                    https://api.kilonova.space/
                    https://api.faststars.space/
                """                
                
                _res = self.buildurl(db, tname, timeout, keys, argkeys)
                if _res is None: return
                
                for _name in _res:
                        assert 'photometry' in _res[_name], _res
                        assert len(_res[_name]['photometry']) > 0
                        self.lc = Table(rows=_res[_name]['photometry'], names=keys)
                        
                self.tab2float(keysfloat)
                if timezero is None:                        
                        self.lc['time'] -= min(self.pps(pps=True)['time'])
                else:
                        assert type(timezero) is float
                        self.lc['time'] -= timezero
                        
        def buildurl(self, db, tname, timeout, keys2query, argkeys):
                
                # base url                                                
                url = '%s%s/photometry/' % (db, tname)
                url = url.replace(' ','')
                
                # query keys
                for _nkey, _key in enumerate(keys2query):
                        if _nkey == len(keys2query)-1: url += '%s'%_key
                        else: url += '%s+'%_key

                # constrains
                _init = True
                for _key in keys2query:
                        if _key in argkeys:
                                if _init:
                                        url += '?'
                                        _init = False
                                else:
                                        url += '&'
                                url += '%s'%_key
                        
                # complete url
                self.logger.info (url)
                try:
                        r = requests.get(url, timeout=timeout)
                        r.raise_for_status()
                except requests.RequestException as e:
                        self.logger.info (e)
                else:
                        _res = r.json()                        
                        return _res

        def readlc(self, filename, filetype='npz', split=' ',
                   keys = ['time','magnitude','band','upperlimit',
                           'e_time','e_lower_time','e_upper_time',
                           'e_magnitude','e_lower_magnitude','e_upper_magnitude'],
                   keysfloat = ['time', 'magnitude'], timezero=None, **kwargs):
                """read lc from file
                """                
                kwargs = self.setkeys(kwargs)

                if not self.lc is None and not kwargs['clobber']:
                        self.logger.info ('Warning: lc data already parsed')
                        return
                
                _r = self.readlist(split=split, filename=filename,
                        filetype=filetype, keys=keys, **kwargs)
                if _r is None: return
                
                _data, _rk = [], []
                for _key in keys:                                            
                        try: _data.append(np.array(_r[_key]))
                        except: _rk.append(_key)
                for _key in _rk:
                        keys.remove(_key)                
                self.lc = Table(_data, names=keys)
                self.tab2float(keysfloat)
                if timezero is None:
                        self.lc['time'] -= min(self.pps(pps=True)['time'])
                else:
                        assert type(timezero) is float
                        self.lc['time'] -= timezero
                        
        def savelc(self, filename, filetype='npz',split=' ',
                   keys = ['time','magnitude','band','upperlimit',
                           'e_time','e_lower_time','e_upper_time',
                           'e_magnitude','e_lower_magnitude','e_upper_magnitude'],
                   **kwargs):
                """save lc to a file                
                """
                kwargs = self.setkeys(kwargs)
                
                if self.lc is None:
                        self.logger.info ('Warning: lc not parsed')
                        return
                
                self.writelist(datain=self.lc, filename=filename,
                        filetype=filetype, split=split, keys=keys, **kwargs)

        def tab2float(self, keys):

                for key in keys:                        
                        if key in self.lc.colnames:                                
                                _idx = []
                                for nn, ii in enumerate(self.lc[key]):
                                        try:
                                                float(ii)
                                                _idx.append(nn)
                                        except:
                                                pass
                                self.lc = self.lc[np.array(_idx)]
                                self.lc[key] = self.lc[key].astype(float)
                                
        def pps(self, pps=True, band=None):
                                   
                if pps:
                        _pps = ''
                else:
                        _pps = 'True'
                if band is None:
                        _id = np.where(self.lc['upperlimit']==_pps)
                else:
                        _id = np.logical_and(self.lc['upperlimit']==_pps,
                                             self.lc['band'] == band)                        
                return self.lc[_id]        
                
        def errors(self, pps=True, band=None):
                
                _xerru, _xerrl, _yerru, _yerrl = [], [], [], []                                
                _data = self.pps(pps=pps, band=band)                
                for _ndd, _dd in enumerate(_data):
                        if _dd['e_time']!='':
                                try: 
                                        _xerru.append(float(_dd['e_time']))
                                        _xerrl.append(float(_dd['e_time']))
                                except:
                                        _xerru.append(0)
                                        _xerrl.append(0)
                        elif _dd['e_lower_time']!='' and _dd['e_upper_time']!='':
                                try:
                                        _xerru.append(float(_dd['e_upper_time']))
                                        _xerrl.append(float(_dd['e_lower_time']))
                                except:
                                        _xerru.append(0)
                                        _xerrl.append(0)
                        else:
                                _xerru.append(0)
                                _xerrl.append(0)
                                                
                        if _dd['e_magnitude']!='':
                                try:
                                        _yerru.append(float(_dd['e_magnitude']))
                                        _yerrl.append(float(_dd['e_magnitude']))
                                except:
                                        _yerru.append(0)
                                        _yerrl.append(0)
                        elif _dd['e_upper_magnitude']!='' and \
                             _dd['e_lower_magnitude']!='':
                                try:
                                        _yerru.append(float(_dd['e_upper_magnitude']))
                                        _yerrl.append(float(_dd['e_lower_magnitude']))
                                except:
                                        _yerru.append(0)
                                        _yerrl.append(0)
                        else:
                                _yerru.append(0)
                                _yerrl.append(0)
                _xerru, _xerrl, _yerru, _yerrl = np.array(_xerru), np.array(_xerrl),\
                        np.array(_yerru), np.array(_yerrl)
                _xerru[np.where(_xerru > 1.)] = 0
                _xerrl[np.where(_xerrl > 1.)] = 0
                _yerru[np.where(_yerru > 1.)] = 0
                _yerrl[np.where(_yerrl > 1.)] = 0                                
                return  [_xerru, _xerrl], [_yerru, _yerrl] 

        def limdist(self, limmag=20, filt='r', dtemp=40):
                
                if self.lc is None:
                        self.logger.info ('lc data not parsed')
                        return
                assert 'time' in self.lc.keys()
                assert 'magnitude' in self.lc.keys()
                assert 'band' in self.lc.keys()        
                                
                _pps = self.pps(pps=True, band=filt)                
                mag0 = min(_pps['magnitude'])
                return 10 ** ((limmag-mag0)/5.) * dtemp
