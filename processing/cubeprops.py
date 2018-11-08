import numpy as np
from astropy.io import fits

def safe_divide(a, b):
    """
        Element-wise division of two arrays a/b, where 0 is returned for
        those elements where b is 0
        @param a: numerator
        @param b: denominator
        @return: a/b if b != 0, otherwise 0
    """
    return np.divide(a, b, out=np.zeros_like(a), where=(b!=0))

class CubeProps:
    def __init__(self, catalogue_cube, data_cube, mult_factor=1):
        """
            @param catalogue_cube: FITS datacube with cluster ids
            @param data_cube: FITS catalogue (of same size as catalogue_cube)
                   containing the densities
            @param mult_factor: factor with which data_cube will be multiplied.
                   of use, if data_cube is density contrast and we want absolute density
        """

        # Load data
        self.data_file = data_cube
        self.catalogue_file = catalogue_cube
        self.dcube = mult_factor * fits.open(data_cube)[0].data.flatten()
        self.cat = fits.open(catalogue_cube)[0].data.flatten().astype("int")

        # Store shape of input data
        self.npixels = np.array([self.dcube.shape])[-1] # number of pixels in fits

        # Object ID data
        self.object_ids = np.unique(self.cat)
        self.nobjects = self.object_ids.max() + 1  # this is not just .size since
                                                   # some ids may be left out, e.g.
                                                   # ids = [0, 1, 2, 4]

        # Print some stats
        self.print_stats()

    def print_stats(self):
        """
            Print all kinds of values contained in class
        """
        print ""
        print "FITS Files"
        print "-- Catalogue file:", self.catalogue_file
        print "-- Datacube file:", self.data_file
        print "-- Number of pixels:", self.npixels
        print ""
        print "Catalogue info"
        print "-- Number of objects:", self.object_ids.size
        print "-- Maximal index:", self.nobjects
        print ""
        print "Data info"
        print "-- Minimal value:", self.dcube.min()
        print "-- Maximal value:", self.dcube.max()
        print "-- Total number of values:", self.dcube.size
        print ""

    def sanity_check(self, density_thres):
        """
            Check the relation
              \sum_{pixel p \in filaments} density_p >= density_thres * #{pixels in filaments}
            And print (Fulfilled?, expression left hand side, expression right hand side)
            @param density_thres: density threshold used to compute pixels that belong in filaments
            @return: None
        """
        npixels_in_filament = 0
        sum_density_p = 0
        for id, val in zip(self.cat, self.dcube):
            if id > 0:
                npixels_in_filament += 1
                sum_density_p += val

        lhs = sum_density_p
        rhs = density_thres*npixels_in_filament

        if lhs >= rhs:
            print "Passed sanity check! :)"
            print "Sum of densities in filament =", lhs
            print "Density threshold * Number of pixels in a filament =", rhs

    def size_hist(self):
        """
          Return count of how many objects with given objects id exist
          @return: object_ids, counts
        """
        # skip id 0 (= void)
        counts = [np.sum(self.cat == id) for id in self.object_ids[1:]]
        return self.object_ids, counts
          
    def avg(self):
        """
            Compute density average for voids and filaments
            @return: dictionary with density averages for:
                    * each complex ["Avgs"]
                    * voids ["VoidAvg"]
                    * filaments ["FilAvg"]
                    * haloes ["HaloAvg"]
                    * overall ["OverallAvg"]
                    as well as volumes:
                    * replace "Avg" <- "Vol"
        """

        data = {}

        sums = np.zeros(self.nobjects) # sum of values in one object
        pixels = np.zeros(self.nobjects) # number of pixels and object contains
        for id, val in zip(self.cat, self.dcube):
            """
            # For the case that sth is really wrong, check this:
            if id == 0:
              print "v,", val
            elif id == self.cat.max():
              print "h,", val
            else:
              print "f,", val
            """
            sums[id] += val
            pixels[id] += 1


        data["Avgs"] = safe_divide(sums, pixels) # avgs for all pixels
        data["VoidAvg"] = data["Avgs"][0] # avg for void
        data["FilAvg"] = safe_divide(np.sum(sums[1:-1]), np.sum(pixels[1:-1])) # avg for filaments
        data["HaloAvg"] = data["Avgs"][-1] # avg for haloes
        data["OverallAvg"] = safe_divide(np.sum(sums), np.sum(pixels)) # total

        # second run for the variance
        varsums = np.zeros(self.nobjects) # sum of values in one object
        means = np.zeros_like(varsums) 
        means[0] = data["VoidAvg"]
        means[-1] = data["HaloAvg"]
        means[1:-1] = data["FilAvg"]
        for id, val in zip(self.cat, self.dcube):
            varsums[id] += (val - means[id])**2

        data["Vars"] = safe_divide(varsums, pixels) # vars for all objects
        data["VoidVar"] = data["Vars"][0] # avg for void
        data["FilVar"] = safe_divide(np.sum(varsums[1:-1]), np.sum(pixels[1:-1])) # vars for filaments
        data["HaloVar"] = data["Vars"][-1] # avg for haloes

        data["Vols"] = pixels # number of pixels for every object
        data["VoidVol"] = pixels[0] # for void
        data["FilVol"] = np.sum(pixels[1:-1]) # for filaments
        data["FilVolAvg"] = np.mean(pixels[1:-1]) # for filaments
        data["FilVolVar"] = np.var(pixels[1:-1]) # for filaments
        data["HaloVol"] = data["Avgs"][-1] # avg for haloes
        data["OverallVol"] = np.sum(pixels) # total

        return data

    def cluster_sizes(self, snap_cat, bins):
        """
          Compute sizes of different clusters in pixels
        """
        pass

