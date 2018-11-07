"""
  Convert FITS file to csv, where the value of the FITS file is > 0.
  Furthermore classify the values into bins.
"""

from astropy.io import fits
import numpy as np
import sys

if len(sys.argv) < 2:
  print "Must give an input as cmd-line argument!"
  print "Optionally after that you can specify an output file."
  exit

if len(sys.argv) < 3:
  out = "data.csv"
else:
  out = sys.argv[2]

inp = sys.argv[1]
data = fits.open(inp)[0].data

thres = 0
idx = np.where(data > thres)
X = np.empty((idx[0].size, 4))
# i = [0,1,2]
# coord = [x,y,z]
# Store coordinates
for i, coord in enumerate(idx):
  X[:,i] = coord

# Store labels
labels = data[idx].astype(int) # cluster labels
X[:,3] = labels

# Compute bins
bins = [1e10, 1e5, 1e4, 1e3, 200, 100, 30] # descending!

# for every label count how often it appears
counts = np.empty(labels.max() + 1)
for l in labels:
  counts[l] += 1

# collapse count to bin values 
bin_counts = np.zeros_like(counts)
for val, n in enumerate(bins):
  bin_counts[np.where(counts < n)] = val + 1

# new labels, clusters are now separated in bins
bin_labels = bin_counts[labels]

Y = np.hstack((X, bin_labels.reshape((-1,1))))
np.savetxt(out, Y, delimiter=",", header="x,y,z,vals,bins")
