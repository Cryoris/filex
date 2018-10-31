#!/usr/bin/python
from astropy.io import fits
import numpy as np
import sys 

print "Run this with -i option, fits data will be stored in the variable called `data`"
if len(sys.argv) < 2:
  print "Please provide the fits file."
else:
  data = fits.open(sys.argv[1])[0].data
  print "Max. val:", data.max()
  print "Min. val:", data.min()
  print "Num. different values:", np.unique(data).size
