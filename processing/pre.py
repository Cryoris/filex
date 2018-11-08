import sys
import numpy as np
from cubeprops import CubeProps

# params from cmd line
if len(sys.argv) < 5:
    print "Usage: python pre.py <snap_id> <filament_threshold> " \
          "<halo_threshold> <min_nvoxels>\n" \
          "Use <snap_id> = 0 for all snapshots (except 001, it's buggy)\n" \
          "Examples: * python pre.py 0 8 50 100\n" \
          "          * python pre.py 012_z003p017 8 50 100"
    exit(1)

snap_id = sys.argv[1]
f_thres = int(sys.argv[2])
h_thres = int(sys.argv[3])
nvox = int(sys.argv[4])

opt = "f%d_h%d_v%d" % (f_thres, h_thres, nvox)

print "\nChosen options\n" \
      "-- snap: " + ("all snaps" if snap_id == "0" else snap_id) + "\n" \
      "-- filament threshold: " + str(f_thres) + "\n" \
      "-- halo threshold: " + str(h_thres) + "\n" \
      "-- min. num. voxels: " + str(nvox)

# load data
if snap_id == "0":
    all_files = ["001_z015p132",
                 "002_z009p993",
                 "003_z008p988",
                 "004_z008p075",
                 "005_z007p050",
                 "006_z005p971",
                 "007_z005p487",
                 "008_z005p037",
                 "009_z004p485",
                 "010_z003p984",
                 "011_z003p528",
                 "012_z003p017",
                 "013_z002p478",
                 "014_z002p237",
                 "015_z002p012",
                 "016_z001p737",
                 "017_z001p487",
                 "018_z001p259",
                 "019_z001p004",
                 "020_z000p865",
                 "021_z000p736",
                 "022_z000p615",
                 "023_z000p503",
                 "024_z000p366",
                 "025_z000p271",
                 "026_z000p183",
                 "027_z000p101",
                 "028_z000p000"]

    redshift = [15.132,
                9.993,
                8.988,
                8.075,
                7.050,
                5.971,
                5.487,
                5.037,
                4.485,
                3.984,
                3.528,
                3.017,
                2.478,
                2.237,
                2.012,
                1.737,
                1.487,
                1.259,
                1.004,
                0.865,
                0.736,
                0.615,
                0.503,
                0.366,
                0.271,
                0.183,
                0.101,
                0.000]

    # exclude first file since it's buggy (I think there are no objects detected)
    files = list(np.array(all_files)[1:])
    # redshift as xlabel
    labels = np.array(redshift[1:]).reshape((-1,1))
    print "x-axis labels:", labels
else:
    labels = np.array([snap_id]).reshape((-1,1))
    files = np.array([snap_id])


snap_dir = "/net/astrogate/export/astrodata/jgacon/filex/snapshots/"
data_cubes = [snap_dir + "snapshot_%s/cubex/temperature_restricted_density.fits" % f for f in files]
catalogue_cubes = [snap_dir + "snapshot_%s/cubex/%s/ids_%s.fits" % (f, opt, f) for f in files]
object_catalogues = [snap_dir + "snapshot_%s/cubex/%s/snap.cat" % (f, opt) for f in files]
output_avgs = "export/%s_avgs.csv" % opt
output_mask_objs= "export/%s_objs_%s.csv"
output_clus = "export/%s_clus.csv" % opt

nparams = 13
data_avgs = np.empty((len(catalogue_cubes), nparams))

bins = np.array([1e10, 1e5, 1e4, 1e3, 200, 100, 30]).astype(int)
data_clus = np.empty((len(object_catalogues), bins.size))

# ENUMERATE AND ZIP
for i, (cat, data_cube, obj_cat) in enumerate(zip(catalogue_cubes, data_cubes, object_catalogues)):
  eng = CubeProps(cat, data_cube)
  den = eng.avg()

  #ids, counts = eng.size_hist()
  #plt.figure()
  #plt.hist(counts, bins=30)
  #plt.ylabel("#Objects of size $N$")
  #plt.xlabel("Object size $N$")
  #plt.savefig("hist_%s.png" % all_files[i])

  # averages & volumes of each labelled object
  density_avgs = den["Avgs"].reshape((-1,1))
  vols = den["Vols"].reshape((-1,1))
  data_objs = np.hstack((density_avgs, vols))
  np.savetxt(output_mask_objs % (opt, "snap_%d" % (i+1)), data_objs)

  # volumes are in pixels
  data_avgs[i,:] = np.array([den["OverallAvg"], den["OverallVol"],
                             den["VoidAvg"],    den["VoidVol"],
                             den["FilAvg"],     den["FilVol"],
                             den["HaloAvg"],    den["HaloVol"],
                             den["VoidVar"],
                             den["FilVar"],
                             den["HaloVar"],
                             den["FilVolAvg"], den["FilVolVar"]])

  # fraction of occupied cube: npixels(cluster)/npixels(all filaments)
  #data_clus[i,:] = eng.cluster_sizes(bins, obj_cat)/den["FilVol"]

# attach labels (usually redshift)
print data_avgs.shape
try:
    data_avgs = np.hstack((labels, data_avgs))
else:
    pass
#data_clus = np.hstack((labels, data_clus))

# store
np.savetxt(output_avgs, data_avgs)
#np.savetxt(output_clus, data_clus)
