#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/kobe.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2020
# Last Modified Date: 01.03.2020
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
from kobe import telescope, triggers, tilings, galaxies

__all__ = ('main')

class kobe(telescope, triggers, tilings, galaxies):

    defkwargs = {**telescope.defkwargs,
                 **triggers.defkwargs,
                 **tilings.defkwargs,
                 **galaxies.defkwargs}
    
    def __init__(self, logger=None, **kwargs):        
        super(kobe, self).__init__(logger=logger, **kwargs)        
