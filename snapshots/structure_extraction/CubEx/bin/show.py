import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits

#input = "snap.fits"
thres = 0
input = "snap.Objects_Id.fits"
data = fits.open(input)[0].data

idx = np.where(data > thres)
c = np.log(data[idx])
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(idx[0], idx[1], zs=idx[2], s=4, c=c)
fig.savefig("objects.png")
