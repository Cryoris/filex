import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import scipy.stats as st
import seaborn as sns
from snapid2z import snapid2z

def sizes(options):
    """
        @brief Visualise the sizes of the filaments throughout the redshifts
        @param options: dict, of form "mode": "saveas".
            mode: "err", "dist", "hist"
            saveas: plt.savefig(saveas)
            Example: {"err": "errplot.png", "dist": "distplot.png"}
    """
    # import data
    pixels = dict() # volumes are given in #pixels
    snap_mask = "/net/astrogate/export/astrodata/jgacon/filex/processing/" \
                "export/f8_h50_v100_objs_snap_%d.csv"
    snap_ids = np.arange(2,28+1)
    z = snapid2z(snap_ids)
    print z

    for id in snap_ids:
        snap = snap_mask % (id - 1) # fix: snap number one too low in filename
        pixels[id] = np.genfromtxt(snap)[1:-1,1] # row 2 contains volumes
                                                 # rm void & halo volumes

    # visualise
    if "err" in options.keys():
        nums = np.array([pixels[id].size for id in snap_ids])
        avgs = np.array([np.mean(pixels[id]) for id in snap_ids])
        mods = np.array([st.mode(pixels[id])[0][0] for id in snap_ids])
        meds = np.array([np.median(pixels[id]) for id in snap_ids])
        stds = np.array([np.std(pixels[id]) for id in snap_ids])

        print mods
        print mods.shape

        plt.figure()
        plt.title("Sizes of filaments as function of redshift")
        plt.xlabel("Redshift $z$")
        plt.xticks(snap_ids[::3], z[::3])

        plt.ylabel("Size in #pixels")

        plt.errorbar(snap_ids, avgs, yerr=stds, label="Mean")
        plt.plot(snap_ids, mods, "g", label="Mode")
        plt.plot(snap_ids, meds, "c", label="Median")
        plt.legend(loc="best")

        plt.twinx()
        plt.ylabel("#Filaments", color="r")
        plt.tick_params("y", colors="r")

        plt.plot(snap_ids, nums, "r--")

        plt.savefig(options["err"])

    if "dist" in options.keys():
        targets = np.array([5,10,15,20,25])
        plt.figure()
        plt.title("Volume distribution of filaments")
        plt.xlabel("Volume $V$ in #pixels")
        plt.ylabel("#Element with $V$ / Total #Elements")
        plt.xlim([0,1000])
        for target in targets:
            sns.kdeplot(pixels[target], label="$z$ = %f" % snapid2z(target))
        plt.legend(loc="best")
        plt.savefig(options["dist"])

    if "dist_inter" in options.keys():
        default = snap_ids[-1]
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        sns.kdeplot(pixels[int(default - 2)], ax=ax)
        plt.xlim([0, 1000])
        plt.ylim([0, 0.01])

        axcolor = 'lightgoldenrodyellow'
        axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
        sid = Slider(axfreq, "ID", 2, 28, valinit=default, valstep=1)

        def update(val):
            id = sid.val

            print id
            ax.clear()
            ax.set_xlim([0,1000])
            ax.set_ylim([0, 0.01])
            sns.kdeplot(pixels[int(id)], ax=ax)
            fig.canvas.draw_idle()
        sid.on_changed(update)

        plt.show()


    if "hist" in options.keys():
        conc = None
        for id, vols in pixels.iteritems():
            data = np.empty((vols.size, 2))
            data[:,0] = id
            data[:,1] = vols

            if conc is None:
                conc = data
            else:
                conc = np.vstack((conc, data))

        plt.figure()
        plt.hist2d(conc[:,0], conc[:,1], bins=(snap_ids.size, 1000))
        plt.ylim([100,400])
        plt.savefig(options["hist"])
