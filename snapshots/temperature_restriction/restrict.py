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
cubex_dir = snaps_dir + "snapshot_" + target + "/cubex/"

# Target files are of form:
# snap_002_z009p993_256_256_8_box_density
density_fits = target_dir + "snap_" + target + "_256_256_8_box_density.fits"
temperature_fits = target_dir + "snap_" + target + "_256_256_8_box_temp.fits"

# Extract FITS data to 3D numpy array
density_cube = fits.open(density_fits)[0].data
temperature_cube = fits.open(temperature_fits)[0].data

# Compute pixels where temperature condition is *not* met 
T_min = 1e3 
T_max = 5e4
T_mask = np.where((T_max < temperature_cube) | (temperature_cube < T_min))

# Apply temperature restriction
density_cube[T_mask] = 0
temperature_cube[T_mask] = 0

# Store output
fits.PrimaryHDU(density_cube).writeto(cubex_dir + "temperature_restricted_density.fits")
fits.PrimaryHDU(temperature_cube).writeto(cubex_dir + "temperature_restricted.fits")
