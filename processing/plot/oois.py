import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from ooi_data import * # define ooi datasets
filaments = ["f%d" % id for id in np.arange(1,4+1)] # f5 is bad data


ooi_set = devo
#ooi_set = evo
def oois(options, verbose=False):

    # for every snapshot extract the data for the oois
    snap_mask = "/net/astrogate/export/astrodata/jgacon/filex/processing/" \
                "export/f8_h50_v100_objs_snap_%d.csv"

    filament_data = dict()

    for i, num in enumerate(ooi_set.snap_nums):
        snap = snap_mask % (num - 1) # fix: snap numer one too low in filename
        snap_data = np.genfromtxt(snap)

        if verbose:
            print "Snapnumber", num

        for fil in filaments:
            labels_list = ooi_set.get(fil)

            if verbose:
                print "Handling", fil
                print "Labels list:", labels_list

            if len(labels_list) - 1 >= i:    # check if for the filament evolutions
                labels = labels_list[i] # there exist labels for this snapnum
                if verbose:
                    print "Found labels:", labels
            else:
                if verbose:
                    print "No labels found."
                continue
            data = []
            for label in labels:
                if verbose:
                    print "Found local data for label", label, ":", snap_data[label,:]
                data.append(snap_data[label,:])

            if fil in filament_data.keys():
                filament_data[fil].append(data)
            else:
                filament_data[fil] = [data]
            if verbose:
                print "Now we have", filament_data[fil]


    # visualise
    # legend entries
    den_full = Line2D([], [], color="red", linestyle="", marker="o",
                      label="Avg. density")
    den_hollow = Line2D([], [], color="red", linestyle="", marker="o",
                        markerfacecolor="None", label="Density")
    vol_full = Line2D([], [], color="blue", linestyle="", marker="o",
                      label="Avg. Volume")
    vol_hollow = Line2D([], [], color="blue", linestyle="", marker="o",
                        markerfacecolor="None", label="Volume")

    for mode, saveas in options.iteritems():
        if "evo" in mode or "fit" in mode:
            for fil, fil_data in filament_data.iteritems():
                x = ooi_set.snap_nums
                vol, den = [], []
                # more than 1 filament (splitted)? --> also plot the average
                vol_mult, den_mult = [], [] # this is stored here

                den_handles = [den_full]
                vol_handles = [vol_full]

                if verbose:
                    print fil, "\nData:", data

                has_splits = False
                for i, data in enumerate(fil_data):
                    # fil_data = [data, data, ...]
                    # data = [array(den,vol), array(den,vol)]
                    den_sum, vol_sum = 0, 0
                    # if filaments have split up we have multiple labels for the
                    # same snap number. If that's the case the marker is hollow
                    # (markerfacecolor="None"), otherwise colored.
                    den_mfc = "None" if len(data) > 1 else "red"
                    vol_mfc = "None" if len(data) > 1 else "blue"
                    for d in data:
                        den.append([x[i], d[0], den_mfc])
                        vol.append([x[i], d[1], vol_mfc])

                        den_sum += d[0]
                        vol_sum += d[1]

                    if len(data) > 1:
                        has_splits = True
                        den.append([x[i], float(den_sum)/len(data), "red"])
                        vol.append([x[i], float(vol_sum)/len(data), "blue"])

                plt.figure()
                ax = plt.subplot(2,1,1)
                ax2 = plt.subplot(2,1,2)
                for rho, v in zip(den, vol):
                    ax.plot(rho[0], rho[1], "ro", markerfacecolor=rho[2])
                    ax2.plot(v[0], v[1], "bo", markerfacecolor=v[2])

                # compute fit to averages
                if "fit" in mode:
                    # get averages (they have markerfacecolor=red/blue)
                    den_avg = np.array([d[:2] for d in den if d[2] == "red"])
                    vol_avg = np.array([d[:2] for d in vol if d[2] == "blue"])


                    den_fit = np.polyfit(den_avg[:,0], den_avg[:,1], 1)
                    vol_fit = np.polyfit(vol_avg[:,0], vol_avg[:,1], 1)

                    ax.plot(den_avg[:,0], np.polyval(den_fit, den_avg[:,0]),
                            "r--")
                    ax2.plot(vol_avg[:,0], np.polyval(vol_fit, vol_avg[:,0]),
                             "b--")
                    den_handles.append(Line2D([], [], color="red",
                                              linestyle="--",
                                              label="Fit w/ %f" % den_fit[0]))
                    vol_handles.append(Line2D([], [], color="blue",
                                              linestyle="--",
                                              label="Fit w/ %f" % vol_fit[0]))

                #ax.set_xlabel("Snapnumber")
                ax.set_ylabel("Density")
                if has_splits:
                    den_handles.append(den_hollow)
                    vol_handles.append(vol_hollow)

                ax.legend(handles=den_handles)
                ax2.legend(handles=vol_handles)
                ax2.set_xlabel("Snapnumber")
                ax2.set_ylabel("Volume")

                plt.tight_layout()
                plt.savefig(fil + "_" + saveas)
