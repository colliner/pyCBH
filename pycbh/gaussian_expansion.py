import numpy as np
import math
'''
distance r valued on Gaussian basis
  exp(−(r − r0)**2 /σ**2)
  where r0 takes values at 20 locations linearly
  placed between 0 and 4, and the
  width σ = 0.5
'''

def gaussian_expansion(d,
                       n_centers=20,
                       sigma=0.5,
                       min_d=0.0,
                       max_d=4.0,
                       centers=None):
    """
    Computes the Gaussian expansion values for a given distance.

    Parameters:
    - d (float): The distance for which the Gaussian expansion is computed.
    - n_centers (int): The number of Gaussian centers.
    - sigma (float): The width of the Gaussian function.
    - min_d (float): The minimum distance for the Gaussian centers.
    - max_d (float): The maximum distance for the Gaussian centers.
    - centers (list or array-like): The pre-defined Gaussian centers. If None, centers are linearly spaced between min_d and max_d.

    Returns:
    - tuple: A tuple containing two elements:
        - list of float: The Gaussian expansion values for the given distance.
        - list: The Gaussian centers used for the expansion.
    """

    if centers is None or None in centers:
        centers = np.linspace(min_d, max_d, n_centers)
    return [math.exp(-(d - x)**2 / sigma**2) for x in centers], list(centers)


if __name__ == '__main__':
    centers = [None]

    distances = [0.5, 0.7, 1.4, 4.0]

    for dist in distances:
        exp_dist, centers = gaussian_expansion(dist, centers=centers)
        print('\nGauExp dist for {}A :\n  {}\n'.format(
            dist, [round(x, 5) for x in exp_dist]))
