# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
KOBE is a package intended to contain core functionality and some
common tools needed for performing pointings scheduler for telescopes with
Python, in order to search for Kilonovae in MMA erea.
"""
from .KButils import utils
from .KBplotter import visualization, vtrigger, \
    vtilings, vgalaxies, vcandidates
from .KBcirculate import circulate
from .KBcandidate import candidates
from .KBtrigger import trigger
from .KBpointings import pointings, tilings, galaxies
from .KBtelescope import telescope
from .KBobservatory import observatory
from .KBschedule import frame, schedule

from .__version__ import version, description
from . import KButils, KBplotter, KBcirculate, \
    KBcandidate, KBtrigger, KBpointings, KBtelescope,\
    KBobservatory, KBschedule
__all__ = KButils.__all__ + KBplotter.__all__ +\
    KBcirculate.__all__ + KBcandidate.__all__ +\
    KBtrigger.__all__ + KBpointings.__all__ +\
    KBtelescope.__all__ + KBobservatory.__all__ +\
    KBschedule.__all__
__version__ = version
__description__ = description
