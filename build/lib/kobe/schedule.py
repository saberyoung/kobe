#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/schedule.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2020
# Last Modified Date: 01.03.2020
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
from kobe import telescope, visualization

__all__ = ('schedule')

class schedule(visualization):
    
    def __init__(self, logger=None, **kwargs):
        
        super(schedule, self).__init__(logger=logger, **kwargs)        
        if not hasattr(self, 'tellist'):
            self.tellist = {}
    
    def add_telescope(self, schedule='T', **kwargs):
        
        _tel = telescope(**kwargs)
        _tel.set_name(**kwargs)
        if _tel.telname is None or _tel.location is None:
            return
        self.tellist[_tel.telname] = _tel.location

    def del_telescope(self, telname):
        if telname in self.tellist:
            del self.tellist[telname]
        elif telname == 'all':
            self.tellist.clear()
        else:
            self.logger.info ('Warning: %s not in tellist' % telname)

    def get_telescope(self,telname, **kwargs):
        if telname in self.tellist:
            kwargs['location'] = self.tellist[telname]
            kwargs['telname'] = telname
            return telescope(**kwargs)
        else:
            self.logger.info ('Warning: %s not in tellist' % telname)
