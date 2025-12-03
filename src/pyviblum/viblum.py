import numpy as np
import random
import math
import itertools
import scipy.special
from scipy import constants
import pyviblum.dictionaries as dicts
import pyviblum.version as version
import pyviblum.read as read
import pyviblum.grad as grad


def print_welcome():

  print (" ")
  print ('Welcome to')
  print("  ___________                      ___________")
  print(" /__________/                      \__________\\")
  print("  ___________   p y  V i b  L u m  ___________ ")
  print(" /__________/                      \__________\\")
  print (" ")
  print ('        Code for simulation of vibronic emission')
  print ('       in noncentrosymmetric lanthanide complexes')
  print(" ")
  print ('                   Authors:')
  print ('               Vsevolod Dergachev')
  print ('               Liviu Chibotaru')
  print ('               Sergey Varganov')
  print (" ")
  print ('*'*55)
  print ('*'*7 + ' '*14 + 'Version: ' + version.get_version() + ' '*15 + '*'*7)
  print ('*'*55)
  print (" ")


def print_input(**myprint):

  print (" ")
  print ('-'*55)
  print ('-'*7 + ' '*17 + 'Input: ' + ' '*17 + '-'*7)
  print ('-'*55)
  mypaths = myprint.get("paths")
  print ('Path:') 
  print ('   to coordinates: '     + str(mypaths["coordinates"])) 
  print ('   to hessian: '         + str(mypaths["hessian"])) 
  print ('   to electronic data: ' + str(mypaths["electronic"])) 
  print ('   to gradients: '       + str(mypaths["gradients"])) 
  print ('   to couplings: '       + str(mypaths["couplings"])) 
  print (" ")
  myparams = myprint.get("params")
  print ("Simulation parameters:")
  print ("   linewidth: {:.2f}".format(myparams["line_width"]))
  print ("   min energy: {:.2f}".format(myparams["min_energy"]))
  print ("   max energy: {:.2f}".format(myparams["max_energy"]))
  print ("   step: {:.5f}".format(myparams["step"]))
  print ('   max number of excitations per mode: ' + str(myparams["num_exct"]))
  print ('   number of modes: '       + str(myparams["num_modes"]))
  print ('-'*55)
  print (" ")
  print (" ")
  print ('-'*55)
  print ("                     Simulation")
  print ('-'*55)


def simulate_spectrum(**mysim):
  """ Calculates vibronic intensities """

  # Get constants dictionary
  myconsts = dicts.get_constants()


  # Stack spin-orbit gradients to create a matrix
  grads = mysim.get("sogradients")
  grads = np.stack([grads[k] for k in grads])       # shape (nstates, 3*natoms)
  
  
  # Unpack vibdata
  myvib = mysim.get("vibdata")
  mass =  myvib["masses"]                           # shape (3*natoms)
  modes = myvib["eigvecs"]                          # shape (3*natoms, 3*natoms-6)
  evals = myvib["eigvals"]                          # shape (3*natoms-6)
  sqrtm = np.sqrt(mass)

 
  # Mass-weight the gradients
  mwgrads = np.divide(grads,sqrtm)                  # shape (nstates, 3*natoms)


  # Project gradients on normal modes  
  qgrads = np.dot(mwgrads,modes)                    # shape (nstates, 3*natoms-6)


  # Get equilibrium positions
  qeq = -1.0*np.divide(qgrads,evals)                # shape (nstates, 3*natoms-6)


  # Calculate adiabatic energies
  energies = mysim.get("energies")                  # shape (nstates)
  shift = -0.5*np.sum(np.divide(np.square(qgrads),evals),axis=1)
  energies += shift
  energies -= energies[-1] 


  # Calculate Huang-Rhys factors
  S = 0.5*np.sqrt(evals)*np.square(qeq[-1] - qeq)   # shape (nstates, 3*natoms-6)


  # Prepare to calculate vibronic intensities
  strengths = mysim.get("strengths")                # shape (nstates)
  myparams = mysim.get("simparams")
  sigma = myparams["line_width"]
  minenergy = myparams["min_energy"]
  maxenergy = myparams["max_energy"]
  step = myparams["step"]
  nex = myparams["num_exct"]
  nmodes = myparams["num_modes"]

  sigma *= myconsts["cm-1_to_ev"]
  sigmasquare = np.square(sigma)


  # Calculate vibronic spectrum
  npoints = int(math.ceil((maxenergy - minenergy)/step))
  energy = np.array([minenergy + float(i)*step for i in range(npoints)])
  intensity = np.zeros(npoints)
  vibronic = np.zeros(npoints)
  zpl = np.zeros_like(intensity)
  tmp = np.zeros_like(intensity)
  nstates = len(qeq[:,0])


  # Loop over spin-orbit states
  for state in range(nstates-1):

     # sort Huang-Rhys (HR) factors in descending order and keep only specified number of modes
     idxS = S[state].argsort()[::-1]
     Ssort = S[state][idxS][0:nmodes]
     exponent = np.exp(-1.0*np.sum(Ssort))
     evalssort = evals[idxS][0:nmodes]
     freqsort = np.sqrt(np.abs(evalssort))*myconsts["to_cm-1"]
     omegasort = 2.0*freqsort*myconsts["pi"]*myconsts["c"]
    
     # energy gap between two states and corresponding osc. strength
     deltaenergy = abs(energies[state])*myconsts["hartree_to_ev"]
     osc = strengths[state]

     print (" ")
     print ("Transition to state", state+1, "at {:.3f}".format(deltaenergy), "eV",\
             "with oscillator strength of {:.3e}".format(osc))
     print (" ")
     print ('-'*38)
     print ("|  Mode  |  Frequency  |  HR factor  |")
     print ('-'*38)
     
     for mode in range(nmodes):
       print ("|"," ", str(mode+1), " "*2, "|", " "*3, "{:.1f}".format(freqsort[mode]), " "*2, "|",\
              "{:.3e}".format(Ssort[mode]), " ", "|")
       print ('-'*38)


     # Loop over excitations    
     k = 0
     for quanta in itertools.product(range(0,nex+1),repeat=nmodes):
        k += 1
        prod_quanta_fact = np.prod(scipy.special.factorial(quanta))
        prod_S_power_quanta = np.prod(np.power(Ssort,quanta))
        total_product = exponent*prod_S_power_quanta/prod_quanta_fact
        vibsum = myconsts["hbar_ev"]*np.dot(quanta,omegasort)
        rho = np.divide(1.0,np.sqrt(2.0*myconsts["pi"]*sigmasquare))*\
              np.exp(-1.0*np.divide(np.square(deltaenergy - (vibsum + energy)),2.0*sigmasquare))
        tmp += (3.0*myconsts["pi"]*osc/deltaenergy)*total_product*rho
        if k == 1: zpl += tmp
    
     intensity += tmp
     idxS = []
     tmp = np.zeros(npoints)

  vibronic = intensity - zpl
  data = np.array(list(zip(energy, vibronic, intensity)))
  np.savetxt('Spectrum.txt', data,
              fmt='%.7e',
              delimiter='                  ',
              header='Scanning frequency (eV)      Vibronic intensity (a.u.)      Total intensity (a.u.)',
              comments='')


  return None


def collect_data_for_simulation(**simdict):
  """ Prepare all data """

  # Print welcome
  print_welcome()

  # Print input data
  print_input(**simdict)

  # Unpack paths to coordinates and Hessian
  coordpath = simdict.get("paths")["coordinates"]
  hessianpath = simdict.get("paths")["hessian"]


  # Read, mass-weight, and diagonalize molecular Hessian
  myhessian = simdict.get("format")
  vibdata = read.normal_modes(coordpath, hessianpath, myhessian)


  # Unpack path to electronic data
  elecpath = simdict.get("paths")["electronic"]


  # Read electronic data
  states_electronic = simdict.get("elstates")
  states_spin_orbit = simdict.get("sostates")
  mult = simdict.get("spin")
  elecdata = read.electronic(elecpath, states_spin_orbit, states_electronic, mult)


  # Unpack paths to electronic gradients and nacs
  gradspath = simdict.get("paths")["gradients"]
  nacspath = simdict.get("paths")["couplings"]


  # Read electronic gradients and nacs
  numatoms = vibdata["natoms"]
  grads_and_nacs = read.gradients_and_nacs(states_electronic, numatoms, gradspath, nacspath)


  # Calculate gradients of spin-orbit states
  numdims = vibdata["ndims"]
  coeff = elecdata["coefficients"]
  spin_orbit_grads = grad.spin_orbit_grad(numdims, mult, states_electronic, states_spin_orbit, coeff, grads_and_nacs)


  # Wrap data and call simulation
  to_sim = {
    "vibdata": vibdata,
    "sogradients": spin_orbit_grads,
    "energies": elecdata["energies"],
    "strengths": elecdata["strengths"],
    "simparams": simdict.get("params")
  }


  # Run simulation
  simulate_spectrum(**to_sim)

  print (" ")
  print ("Vibronic spectrum has been saved to Spectrum.txt")
  print (" ")
  print ('-'*55)
  print ('-'*7 + ' '*9 + 'pyVibLum ends normally ' + ' '*9 + '-'*7)
  print ('-'*55)

  return None

