def get_zero_cross(absc, ords):
    """

    :param absc: Abscissa (x)
    :param ords: Ordinates (y)

    :return: Return x of 1st zero-crossing. First we find the first place x-values
    between which ordinates change sign and then we interpolate between them.
    """
    parities = np.array([ords[ii] * ords[ii - 1] for ii in range(1, len(ords))])
    crossings = np.where(parities<0)[0]

    if len(crossings)==0:
        print("No zero corssings found")
        return
    else:
        i = crossings[0]
        return np.interp(0, [ords[i-1], ords[i]], [absc[i-1], absc[i]])