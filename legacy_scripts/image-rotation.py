"""
L1k - 1585989326 imsize: 5799 pix
L4kn- 1586305859 imsize: 5962 pix
U1k - 1586162571 imsize: 3684 pix
U4k - 1586070058 imsize: 3696 pix
L4k - 1586074843 imsize: 5804 pix
"""

import atpy, os, sys, glob
import numpy as N
import pylab as pl
from math import *
from astropy.io import fits
#pl.ion()

f     = "/home/samuel/SARAO/imagin_varification/2020-correlator-test/April/images/catalogues/L1K-x-U4K-x-L4Kv2-x-U1K-x-SUMSS.tbl"
nums  = ['L1K 1585989326 - L1K', 'U4K 1586070058 - L1K', 'L4K-2 1586305859 - L1K', 'U1K 1586162571 - L1K', 'SUMSS - L1K']
files = ['L1K', 'U4K', 'L4K', 'U1K']
outdr = os.path.dirname(f)

sizes = N.asarray([5799, 3696, 5962, 3685, 5799])  # repeat l1k since l1k is ref since sumss has fat beam
xcen, ycen = sizes/2, sizes/2
cell = N.asarray([0.0003448295, 0.0005410507, 0.0003354017, 0.0005428091, 0.0003448295])
cell *= 3600

t  = atpy.Table(f)
#t.describe()  # will show column names

ras, decs = [], []
for i in range(4):
  ras.append(t['RA_'+str(i+1)])
  decs.append(t['DEC_'+str(i+1)])
ras.append(t['_RAJ2000'])
decs.append(t['_DEJ2000'])

xpos, ypos = [], []
for i in range(4):
  xpos.append(t['Xposn_'+str(i+1)])
  ypos.append(t['Yposn_'+str(i+1)])
xpos.append(t['Xposn_1'])
ypos.append(t['Yposn_1'])

# approx, is valid only at equator
# theta is PA of src2 from src1
# phi is PA of src from image centre
# no rotation => random; rotation => linear
#f1 = pl.figure(figsize=(9,9))
#f2 = pl.figure(figsize=(9,9))
#f3 = pl.figure(figsize=(9,9))
#f4 = pl.figure(figsize=(9,9))
ra1, dec1 = ras[-1], decs[-1]
ra1, dec1 = ras[0], decs[0]

for i in range(1,5):
  ra2, dec2 = ras[i], decs[i]
  theta = N.arctan2(dec2-dec1, (ra1-ra2)*cos(N.mean(dec1)/180*pi))*180/pi
  phi   = N.arctan2(ypos[i]-ycen[i], xpos[i]-xcen[i])*180/pi
  phi   = N.where(phi<-90,phi+360,phi)
  dist  = N.sqrt((xpos[i]-xcen[i])*(xpos[i]-xcen[i]) + (ypos[i]-ycen[i])*(ypos[i]-ycen[i]))

  #pl.figure()
  pl.subplot(2,2,i)
  pl.scatter(phi, theta, s=dist*1.0/300)
  pl.title(nums[i])
  if i in [1,3]: pl.ylabel('Phi (deg)')
  if i in [3,4]: pl.xlabel('Theta (deg)')

pl.suptitle('PA of src from centre (phi) vs angle of offset (theta)')
pl.savefig(outdr+'/corrtest_rotation_phi_theta.png')
pl.show()


