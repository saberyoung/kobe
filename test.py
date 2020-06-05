import sys
import pylab as pl
import numpy as np
pl.ion()

from kobe import schedule       
a = schedule()
a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
a.set_pointings(strategy='T')
a.pointings.generatep(limdec=[-20,90],fovra=5, fovdec=5)
a.calc_loc_cl_pointings(add=True)

'''
from kobe import schedule
a=schedule()

a.set_pointings(strategy='T')
a.pointings.generatep(fovra=3,fovdec=3,limdec=[-10,10], limra=[0,10])
a.rankp_visibility()


a.fits('/Users/yash0613/Downloads/LALInference.v1.fits.gz')
#a.locshow(cls=[.5,.9])
#a.zoomin([-.3,0],[-.5,-.2])

a.set_pointings(strategy='T')
a.pointings.inputp(ra='00:59:40.20', dec='-24:44:56.00', fovra=1, fovdec=1, hms=True)
b=a.calc_loc_prob_pointings(add=False)

'''


'''
from kobe import schedule
a=schedule()
a.parse_time('20200101')

a.add_observatory('lapalma')
a.set_pointings(strategy='T')
a.pointings.generatep(fovra=3,fovdec=3,limdec=[-90,0])
a.add_telescope('VST')

a.add_observatory('asiago')
a.set_pointings(strategy='T')
a.pointings.generatep(fovra=3,fovdec=3,limdec=[0,90])
a.add_telescope('schmidt')
a.plot_sky()
'''
