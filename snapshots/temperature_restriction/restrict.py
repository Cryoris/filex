"""
  For all pixels in the density cube, only keep it, if the temperature of 
  the according pixel is within specified temperature interval.
"""

import sys
import numpy as np
from astropy.io import fits

if len(sys.argv) < 2:
  print "Usage: python restrict.py <target_id>"
  print "Example: python restrict.py 012_z003p017"
  exit(1)
else:
  target = sys.argv[1] 

# Define needed directories
snaps_dir = "/net/astrogate/export/astrodata/jgacon/filex/snapshots/"
target_dir = snaps_dir + "snapshot_" + target + "/fits3D/"
density_fits = target_dir + "density_cube.fits"
temperature_fits = target_dir + "temperature_cube.fits"
cubex_dir = snaps_dir + "snapshot_" + target + "/cubex/"

# Extract FITS data to 3D numpy array
density_cube = fits.open(density_fits)[0].data
temperature_cube = fits.open(temperature_fits)[0].data

# Compute pixels where temperature condition is met 
T_min = 1e3 
T_max = 5e4
T_mask = np.where((T_min < temperature_cube) & (temperature_cube < T_max))

# Apply temperature restriction to density cube
density_cube[T_mask] = 0

# Store output
fits.PrimaryHDU(density_cube).writeto(cubex_dir + "temperature_restricted_density.fits")
