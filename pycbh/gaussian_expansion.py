import numpy as np
import math

'''
distance r valued on Gaussian basis
  exp(−(r − r0)**2 /σ**2)
  where r0 takes values at 20 locations linearly
  placed between 0 and 4, and the
  width σ = 0.5
'''

def gaussian_expansion(d, n_centers=20, sigma=0.5, min_d=0.0, max_d=4.0, centers=[None]):
  if None in centers:
    centers = np.linspace(min_d, max_d, n_centers)
    #print('{} centers between {} and {} :\n  {}'.format(n_centers,min_d,max_d,centers))
  return [ math.exp(-(d - x)**2 / sigma**2) for x in centers ], centers

if __name__=='__main__':
  centers=[None]

  distances = [0.5, 0.7, 1.4, 4.0]

  for dist in distances:
    exp_dist, centers = gaussian_expansion(dist, centers=centers)
    print('\nGauExp dist for {}A :\n  {}\n'.format(dist,[round(x,5) for x in exp_dist]))

