from numpy import inf, array

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

def snapid2z(id):
    """
        @brief Convert EAGLE snapshot ID to redshift
        @param id: EAGLE snapshot id (int from 1 to 28)
        @return: numpy.array(redshift), if id is valid. -numpy.inf otherwise
    """
    print id
    try:
        res = array(redshift)[id - 1]
    except:
        res = -inf

    return res
