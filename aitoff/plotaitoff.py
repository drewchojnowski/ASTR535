from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# asu.fit = Catalog of Milky Way Globular Clusters
x=fits.open('asu.fit')

ra=x[1].data['raj2000']* u.degree
dec=x[1].data['dej2000']* u.degree
glon=x[1].data['glon']* u.degree
glat=x[1].data['glat']* u.degree

c=SkyCoord(ra=ra,dec=dec,frame='icrs')
ra_rad=c.ra.wrap_at(180*u.deg).radian
dec_rad=c.dec.radian

plt.figure(figsize=(14,9))
plt.subplot(111, projection="aitoff")
plt.title("Aitoff projection of Milky Way globular clusters (RA,DEC)", y=1.08)
plt.grid(True)
plt.plot(ra_rad, dec_rad, 'ro', markersize=7, alpha=0.6)
plt.subplots_adjust(top=0.95, bottom=0.0)
plt.show()

c=SkyCoord(ra=glon,dec=glat,frame='icrs')
ra_rad=c.ra.wrap_at(180*u.deg).radian
dec_rad=c.dec.radian

plt.figure(figsize=(14,9))
plt.subplot(111, projection="aitoff")
plt.title("Aitoff projection of Milky Way globular clusters (GLON,GLAT)", y=1.08)
plt.grid(True)
plt.plot(ra_rad, dec_rad, 'ro', markersize=7, alpha=0.6)
plt.subplots_adjust(top=0.95, bottom=0.0)
plt.show()
