import numpy as np

def spin_orbit_grad(ndims, spin, elstates, sostates, coeff, derivmatrix):
  """ Calculate gradients of spin-orbit states """

  proj = int(2*spin + 1)
  grad = {}

  for state in sostates:
     grad[str(state)] = np.zeros(ndims)
     coeff_name = 'c' + str(state)
     c = coeff[coeff_name].reshape(len(elstates),proj)
     for i in range(1,len(elstates)+1):
        for j in range(i+1, len(elstates)+1):
            grad[str(state)] += np.sum(2.0*np.real(c[i-1].conjugate()*c[j-1]))*derivmatrix[i,j]
        grad[str(state)] += np.sum(np.real(c[i-1].conjugate()*c[i-1]))*derivmatrix[i,i]

  return grad


