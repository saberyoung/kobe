# coding:utf-8

"""Various generic useful functions 
"""
def is_seq(o):
    """Check if the object is a sequence.

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

def is_seq_of_seq(o):
    """Check if the object is a sequence of sequences. No check is done on
    the sizes of the sequences.

    Parameters
    ----------
    o : any object
      The object to check
    
    Returns
    -------
    is_seq_of_seq : bool
      True if *o* is a sequence of sequences, False otherwise.
    """
    if not is_seq(o):
        return False
    for s in o:
        if not is_seq(s):
            return False
    return True

def is_like2d(o):
    """Check if *o* is conformable to a 2d array.

    Parameters
    ----------
    o : any object
      The object to check
   
    Returns
    -------
    is_like2d : bool, scalar
      True if *o* is conformable to a 2d array, False otherwise.
    """
    if not is_seq(o):
        return False
    size = None
    for s in o:
        if not is_seq(s):
            return False
        if size is None:
            size = len(s)
        if len(s) != size:
            return False
    return True

def len_array_or_arrays(o):
    """Returns the length of a single array or list of arrays
    
    Parameters
    ----------
    o : either array or sequence of arrays

    Returns
    -------
    length : length of array
    """
    if is_seq_of_seq(o):
        length = len(o[0])
    else:
        length = len(o)
    return length
    
class readlist:
    
    def __init__(self, logger=None, filename=None, split=' ',
                 keys=['ra','dec','fovra','fovdec','n']):
        
        # ----- define logger ----- #
        if logger is None:
            try:
                import logging                            
                logging.basicConfig(level = logging.INFO)
                self.logger = logging.getLogger(__name__)
            except:
                self.logger = None
        else:
            self.logger = logger
            
        self.filename = filename
        self.split = split
        self.keys = keys
        self.data = {}
        
    def _exists(self):
        import os
        if os.path.exists(self.filename):
            return True
        else:
            if not self.logger is None:
                self.logger.info ('Error: %s not exists'%self.filename)
            return False

    def fits(self):
        from astropy.table import Table
        self.data = Table.read(self.filename, format='ascii')
        
    def txt(self):
        
        if not self._exists(): return                
        _keys, _data = {}, {}        
        for nll, ll in enumerate(open(self.filename).readlines()):
            if ll[0] == '#':
                # read header
                ll = ll.replace('#','').replace('\n','')            
                _nk = 0
                for _k in ll.split(self.split):
                    if len(_k) != 0:
                        _keys[_nk] = _k
                        _data[_k] = []
                        _nk += 1
            else:
                # read data
                _nk = 0
                for _k in ll.split(self.split):                    
                    try:  float(_k)                            
                    except:  continue
                    if _nk in _keys:
                        _data[_keys[_nk]].append(float(_k))
                        _nk += 1
                        
        for _k in self.keys:
            if _k in _data.keys():
                self.data[_k] = _data[_k]
            else:
                if not self.logger is None:
                    self.logger.info ('Warning: key %s not found'%_k)
                
    def npz(self):
        import numpy as np
        try:
            _npdata = np.load(self.filename)
        except:
            if not self.logger is None:
                self.logger.info ('Warning: npz file %s broken'%self.filename)
            
        for _k in self.keys:
            if _k in _npdata.keys():
                self.data[_k] = _npdata[_k]
            else:
                if not self.logger is None:
                    self.logger.info ('Warning: missing keyword %s'%_k)                
