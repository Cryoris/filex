import numpy as np
import matplotlib.pyplot as plt

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
    if "evo" in options.keys():
        for fil, fil_data in filament_data.iteritems():
            x = ooi_set.snap_nums
            vol, den = [], []

            if verbose:
                print fil, "\nData:", data

            for i, data in enumerate(fil_data):
                # fil_data = [data, data, ...]
                # data = [array(den,vol), array(den,vol)]
                for d in data:
                    # if filaments have split up we have multiple labels for the
                    # same snap number
                    den.append((x[i], d[0]))
                    vol.append((x[i], d[1]))

            plt.figure()
            ax = plt.subplot(1,1,1)
            ax2 = ax.twinx()
            for rho, v in zip(den, vol):
                ax.plot(rho[0], rho[1], "ro")
                ax2.plot(v[0], v[1], "bo")

            ax.set_xlabel("Snapnumber")
            ax.set_ylabel("Density")
            ax2.set_ylabel("Volume")

            plt.savefig("b" + fil + ".png")
