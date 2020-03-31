#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/tilings.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import numpy as np
import healpy as hp
import logging
from kobe import pointings, circulate, visualization

__all__ = ['telescope']

class telescope(circulate, visualization):
        """telescope: define the name and a pointing list for one telescope
        inherient from tilings and galaxies
        """
        
        telname = None
        '''
        define telescope name
        '''

        pointings = None
        '''
        define telescope pointings
        '''        
        
        def __init__(self, logger=None, **kwargs):
                """
                initialize tiling object with input logger and parameter settings        
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set keys ----- #        
                self.defkwargs = self.setkeys(kwargs)

        def check_tel(self):
                assert not self.telname is None, 'no telescope name parsed'
                assert not self.pointings is None, 'no telescope pointing defined'
                assert not self.pointings.data is None, 'no telescope pointing parsed'
                
        def set_pointings(self, strategy='T'):                
                self.pointings = pointings(strategy)
                
        def set_telname(self, telname):
                """set telescope name

                Parameters
                ----------
                name :      `string` 
                  input a telescope name
                namemax:    `int` 
                  if not telescope name given, construct with an interger number
                  number selected from 0 tp namemax
                """
                
                if not telname is None:
                        self.telname = telname
                else:                        
                        _nlist = []
                        if hasattr(self, 'telescopes'):
                                for ii in self.telescopes.keys():
                                        if type(ii) is int:
                                                _nlist.append(ii)
                                if len(_nlist) == 0:
                                        self.telname = 0
                                else:
                                        self.telname = max(_nlist)+1
                        else:
                                self.telname = 0                        
