import sys
import pylab as pl
import numpy as np
pl.ion()


from kobe import galaxies
a=galaxies()    
a.generatep(limdist=[0,40], limra=[20, 100], limdec=[0,30])
a.rankp(mode=2)
                
'''
from kobe import pointings
a=pointings('T')
a.generatep(fovra=3,fovdec=3,limra=[20, 100], limdec=[0,30])
a.rankp(mode=2)
a.locshow(color='r')
a.routeshow(color='k')
a.zoomin([-1.2,0],[-.5,.6])
'''




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
