from astropy.io import fits

file = "snap.Object_Ids.fits"
img = fits.open(file)[0].data
print img.shape
