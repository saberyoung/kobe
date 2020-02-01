from pst.pipeline.tilings import PstGetTilings
from pst.pipeline.triggers import PstParseTriggers
from pst.pipeline.galaxies import PstGetGalaxies
from pst.view.show import PstPlotter, interactive_show
import pylab as pl

#a=PstGetTilings()
b=PstParseTriggers()
#c=PstGetGalaxies()
d=PstPlotter()

#def plotfunc(pp):
#    pl.clf()
#    d=PstPlotter(defconf=pp)
#    d.locshow(b,vmax=10**-5, ptype='m')
#    d.contshow(b,nside=16)
#    d.frame()
    
#a.generate(limra=[0,360],limdec=[-20,90],fovra=3.,fovdec=3.)
b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
#c.read()


pl.ion()
           
#pparams = {'theta':0,'phi':0,'psi':0}
#optparams = ['theta','phi','psi']
#interactive_show(plotfunc,pparams,optparams)

            
#theta=0
d.locshow(b,vmax=10**-5, ptype='m',figsize=(6,6))

input('..')
