import sys
import numpy as np
from src.cubeprops import CubeProps
from src.snap_info import SnapInfo

mode = ["production"]

# Snapshot location
fid = "012_z003p017"
f_thres = 3
h_thres = 20
if len(sys.argv) < 5:
    print "Usage: python main.py <snap_id> <mode> <filament_threshold> <halo_threshold>"
    print "To use certain snapshot give a string of form:", fid
    print "Or to use %s just type 0" % fid
    print "Default mode is", mode
    print "Default filament/halo thresholds are %d / %d" % (f_thres, h_thres)
else:
    if not sys.argv[1] == "0":
      fid = sys.argv[1]

if len(sys.argv) >= 3:
  mode = [sys.argv[2]]

if len(sys.argv) >= 5:
  f_thres = int(sys.argv[3])
  h_thres = int(sys.argv[4])

thres = "f%d_h%d_v100" % (f_thres, h_thres)

print "Options"
print "-- snapshot: snap_%s" % fid
print "-- mode:", mode[0]
print "-- f_thres:", f_thres
print "-- h_thres:", h_thres
print ""

folder = "/net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0752/snapshot_" + fid + "/"
fname = "snap_" + fid + ".0.hdf5"
file = folder + fname
#baryon_density = SnapInfo(file).Omega_b # in units of critical density
baryon_density = 1

if "production" in mode or "load" in mode:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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

    # Cubex threshold options:
    # f: filament detection threshold
    # h: halo detection threshold

    files = list(np.array(all_files)[1:])
    #data_cube = "./data/snap.fits" # old cube

    # snapshot number as xlabel
    #xlabels = [f[1:3] for f in files]
    # redshift as xlabel
    xlabels = np.array(redshift[1:])
    print xlabels

    if "production" in mode:
      ref_folder = "/net/astrogate/export/astrodata/jgacon/EAGLE/RefL0025N0752/"
      data_cubes = [ref_folder + "snapshot_%s/cubex/%s/snap.fits" % (f, thres) for f in files]
      catalogue_cubes = [ref_folder + "snapshot_%s/cubex/%s/ids_%s.fits" % (f, thres, f) for f in files]
      object_catalogues = [ref_folder + "snapshot_%s/cubex/%s/snap.cat" % (f, thres) for f in files]
      output_lss = "export_lss_%s.csv" % thres
      output_cs = "export_cs_%s.csv" % thres

      nparams = 13
      data = np.empty((len(catalogue_cubes), nparams))

      bins = np.array([1e10, 1e5, 1e4, 1e3, 200, 100, 30]).astype(int)
      csizes = np.empty((len(object_catalogues), bins.size))

      # ENUMERATE AND ZIP
      for i, (cat, data_cube, obj_cat) in enumerate(zip(catalogue_cubes, data_cubes, object_catalogues)):
          eng = CubeProps(cat, data_cube, mult_factor=baryon_density)
          den = eng.avg()
          ids, counts = eng.size_hist()

          plt.figure()
          plt.hist(counts, bins=30)
          plt.ylabel("#Objects of size $N$")
          plt.xlabel("Object size $N$")
          plt.savefig("hist_%s.png" % all_files[i])

          # volumes are in pixels
          data[i,:] = np.array([den["OverallAvg"], den["OverallVol"],
                                den["VoidAvg"],    den["VoidVol"],
                                den["FilAvg"],     den["FilVol"],
                                den["HaloAvg"],    den["HaloVol"],
                                den["VoidVar"],
                                den["FilVar"],
                                den["HaloVar"],
                                den["FilVolAvg"], den["FilVolVar"]])

          # fraction of occupied cube: npixels(cluster)/npixels(all filaments)
          cfrac[i,:] = eng.cluster_sizes(bins, obj_cat)/den["FilVol"]

      np.savetxt(output_lss, data)
      np.savetxt(output_cs, cfrac)

    elif "load" in mode:
      data = np.genfromtxt("export_lss_%s.csv" % thres)
      cfrac = np.genfromtxt("export_cs_%s.csv" % thres)

    labels = ["Overall", "Void", "Filament", "Halo"]
    style = ["bo", "r^", "gs", "c*"]
    avg_ind = [0, 2, 4, 6]
    var_ind = [8, 9, 10]
    vol_ind = [1, 3, 5, 7]
    filvor_var = [11, 12]
    x = np.arange(data.shape[0])

    plt.figure()
    plt.title("Volume fractions of different cluster sizes")
    for j, b in enumerate(bins):
      plt.plot(x, cfrac[:,j], "o", label=("< %d" % b))

    idx = np.linspace(0, x.size - 1, x.size/4.).astype("int")
    plt.xticks(x[idx], xlabels[idx])
    plt.xlabel("Redshift")
    plt.ylabel("Volume fraction")
    plt.savefig("volume_fractions_%s.png" % thres)

    plt.figure(figsize=(12,8))
    plt.title("Density averages [in $\Omega_b$] (%s)" % thres)
    for ind, stl, label in zip(avg_ind, style, labels):
        plt.semilogy(x, data[:,ind], stl, label=label)
    # avg density
    plt.axhline(y=1, color="k", linestyle="--")
    # filament threshold
    plt.axhline(y=f_thres, color="g", linestyle="--")
    # halo threshold
    plt.axhline(y=h_thres, color="c", linestyle="--")
    plt.legend(loc="best")
    plt.xticks(x[::3], xlabels[::3])
    plt.xlabel("Redshift")
    plt.ylabel("Density constrast")
    plt.savefig("densities_%s.png" % thres)

    plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    ax.set_yscale("log", nonposy="clip")
    #plt.title("Renormalised Density averages [in $\Omega_b$] (%s)" % thres)
    plt.title("Renormalised Density Averages [in $\Omega_b$]")
    for ind, var, stl, label in zip(avg_ind, var_ind, style, labels):
        #ax.semilogy(x, data[:,ind]/data[:,0], stl, label=label)
        #plt.loglog(x, data[:,ind]/data[:,0], stl, label=label)

        # renormalising the data means we need to renormalise the stddev
        # (its actually linear, calculate it!)
        density_normed = data[:,ind]/data[:,0]
        std_normed = np.sqrt(data[:,var])/data[:,0]
        plt.errorbar(x, density_normed, yerr=std_normed)

    # avg density
    #plt.axhline(y=1, color="k", linestyle="--")
    # filament threshold
    #plt.axhline(y=f_thres, color="g", linestyle="--")
    # halo threshold
    #plt.axhline(y=h_thres, color="c", linestyle="--")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2*box.height, box.width, box.height*0.8])
    #plt.legend(loc="center left", bbox_to_anchor=(1,0.5))
    plt.legend(loc="lower center", bbox_to_anchor=(0, -0.3, 1, 0.2) , mode="expand", ncol=4)
    #plt.legend(loc="best")
    #plt.xticks(x, xlabels)
    idx = np.linspace(0, x.size - 1, x.size/4.).astype("int")
    plt.xticks(x[idx], xlabels[idx])
    plt.xlabel("Redshift")
    plt.ylabel("Density constrast/Total density")
    plt.savefig("densities_tight.png")

    plt.figure(figsize=(12,8))
    plt.subplot(211)
    plt.title("Vol-avg filament overdensity (wrt baryons, thres: %d)" % f_thres)
    plt.semilogy(x, data[:,4], "gs")
    #plt.axhline(y=f_thres, color="g", linestyle="--", label="Threshold")
    plt.ylabel("Overdensity")
    plt.xlabel("Redshift")
    plt.xticks(x, xlabels)
    y = np.linspace(data[:,4].min(), data[:,4].max(), 5)
    plt.yticks(y, [str(yl) for yl in y])

    plt.subplot(212)
    plt.title("Vol-avg halo overdensity (wrt baryons, thres: %d)" % h_thres)
    plt.semilogy(x, data[:,6], "c*")
    #plt.axhline(y=f_thres, color="c", linestyle="--", label="Threshold")
    plt.ylabel("Overdensity")
    plt.xlabel("Redshift")
    #plt.xticks(x, xlabels)
    plt.xticks(x[::3], xlabels[::3])
    y = np.linspace(data[:,6].min(), data[:,6].max(), 5)
    plt.yticks(y, [str(yl) for yl in y])

    plt.tight_layout()
    plt.savefig("closeup_%s.png" % thres)

    plt.figure(figsize=(12,8))
    plt.title("Volumes [in #pixels] (%s)" % thres)
    for ind, stl, label in zip(vol_ind, style, labels):
        plt.semilogy(x, data[:,ind], stl, label=label)
    plt.legend(loc="best")
    #plt.xticks(x, xlabels)
    plt.xticks(x[::3], xlabels[::3])
    plt.xlabel("Redshift")
    plt.ylabel("Volumes [in #pixels]")
    plt.savefig("volumes_%s.png" % thres)

#    plt.show()

if "export" in mode:
    import matplotlib.pyplot as plt

    #files = ["v01m30", "v03m30", "v03m100", "v05m30", "v05m100", "v10m30", "v20m30"]
    files = ["test"]
    data_cube = "./data/snap.fits"
    catalogue_cubes = ["./data/objects_id/snap_" + ext + ".fits" for ext in files]

    output = "export.csv"

    nparams = 8
    data = np.empty((len(catalogue_cubes), nparams))

    for i, cat in enumerate(catalogue_cubes):
        eng = CubeProps(cat, data_cube, mult_factor=baryon_density)
        den = eng.avg()

        data[i,:] = np.array([den["OverallAvg"], den["OverallVol"], den["VoidAvg"],
                              den["VoidVol"], den["FilAvg"], den["FilVol"],
                              den["HaloAvg"], den["HaloVol"]])

    np.savetxt(output, data)

    labels = ["Overall", "Void", "Filament", "Halo"]
    style = ["bo", "r^", "gs", "c*"]
    avg_ind = [0, 2, 4, 6]
    vol_ind = [1, 3, 5, 7]
    x = np.arange(data.shape[0])

    plt.figure()
    plt.title("Density averages [in $\\rho_{crit}$]")
    for ind, stl, label in zip(avg_ind, style, labels):
        plt.plot(x, data[:,ind], stl, label=label)
    plt.legend(loc="best")
    plt.xticks(x, files)
    plt.savefig("densities.png")

    plt.figure()
    plt.title("Volumes [in #pixels]")
    for ind, stl, label in zip(vol_ind, style, labels):
        plt.plot(x, data[:,ind], stl, label=label)
    plt.legend(loc="best")
    plt.xticks(x, files)
    plt.savefig("volumes.png")

    plt.show()


if "simple" in mode:
    # fits file with filament complexes
    ref_folder = "/net/astrogate/export/astrodata/jgacon/EAGLE/RefL0025N0752/"
    data_cube = ref_folder + "snapshot_%s/cubex/%s/snap.fits" % (fid, thres)
    catalogue_cube = ref_folder + "snapshot_%s/disperse/%s/id5.fits" % (fid, thres)

    eng = CubeProps(catalogue_cube, data_cube, mult_factor=baryon_density)
    eng.sanity_check(density_thres=f_thres)
    den = eng.avg()
    print "Sum(Density)/Vol of snapshot:", den["OverallAvg"], "Vol:", den["OverallVol"]
    print "Sum(Density)/Vol of void:", den["VoidAvg"], "Vol:", den["VoidVol"]
    print "Sum(Density)/Vol of filaments:", den["FilAvg"], "Vol:", den["FilVol"]
    print "Sum(Density)/Vol of haloes:", den["HaloAvg"], "Vol:", den["HaloVol"]

if "full" in mode:
    # fits file with filament complexes
    catalogue_cube = "./data/objects_id/snap_v3m30.fits"
    data_cube = "./data/snap.fits"

    eng = CubeProps(catalogue_cube, data_cube, mult_factor=baryon_density)
    eng.sanity_check(density_thres=3)
    den = eng.avg()

    print "Sum(Density)/Vol of all complexes:", den["Avgs"]
    print "Vol of all complexes:", den["Vols"]
    print "Sum(Density)/Vol of snapshot:", den["OverallAvg"], "Vol:", den["OverallVol"]
    print "Sum(Density)/Vol of void:", den["VoidAvg"], "Vol:", den["VoidVol"]
    print "Sum(Density)/Vol of filaments:", den["FilAvg"], "Vol:", den["FilVol"]
