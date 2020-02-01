# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
PStools is a package intended to contain core functionality and some
common tools needed for performing pointings scheduler for telescopes with
Python.
"""
from .__version__ import version
from .pipeline import *
from .view import *
from .circulate import *
from .interface import *

#__all__ = triggers.__all__ + tilings.__all__
__version__ = version
