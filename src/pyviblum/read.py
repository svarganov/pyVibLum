import os
from pathlib import Path
import numpy as np
from pyviblum import dictionaries as dicts
import math
import itertools
import io
import re
from scipy import constants


def normal_modes(geomfile, hessfile, hessian):

  # read atomic coordinates
  f = open(geomfile, 'r')
  lines = f.readlines()
  f.close()

  atoms = []
  positions = []
  n = 0

  natoms = int(lines[n])
  n+=1

  for i in range(natoms):
    n+=1
    line = lines[n].split()
    atoms.append(line[0])


  # read, mass-weight, and diagonalize Hessian
  const_dict = dicts.get_constants()

  mass_dict = dicts.get_atomic_masses()
  mass = np.asarray([mass_dict[atom]*const_dict["amu_to_au"] for atom in atoms for i in range(3)])
    
  f = open(hessfile)
  lines = f.readlines()
  f.close()
  n = 0

  # initialize Hessian
  ndims = 3*natoms
  h = np.zeros(shape=(ndims,ndims))

  if hessian == "gamess":
# Read GAMESS-format molecular Hessian
      while n < len(lines):
         line = lines[n]
         if "$HESS" in line:
            n+=1
            for i in range(ndims):   
               for j in range(1,int(math.ceil(float(ndims)/5))+1):
                  n+=1
                  line = lines[n]
                  if j != int(math.ceil(float(ndims)/5)) :
                    last_number = 5
                  else:
                    last_number = 5 - 5*int(math.ceil(float(ndims)/5)) + ndims
                  for l in range(1,last_number+1):
                      k = j*5 + l - 5
                      h[i][k-1] = line[(15*(l-1)+5):(15*l+5)]
         n+=1


  elif hessian == "orca":
# Read ORCA-format molecular Hessian
      while n < len(lines):
        line = lines[n]
        if "$hessian" in line:
          break
        n+=1
      n+=2

      number_of_five_column_blocks = int(math.ceil(float(ndims)/5))

      for i in range(number_of_five_column_blocks):
         if (i*5 <= ndims):
            column = 5
         else:
            column = 5 - 5*int(math.ceil(float(ndims)/5)) + ndims
     
         if (column == 5):
           for j in range(ndims):
               n+=1
               line = lines[n]
               array = line.split()
               h[j][i*5:i*5 + 5] = array[1:6]
         elif (column == 1):
            for j in range(ndims):
               n+=1
               array = line.split()
               h[j][i*5] = array[1]
         elif (column == 2):
            for j in range(ndims):
               n+=1
               array = line.split()
               h[j][i*5:i*5 + 2] = array[1:3]
         elif (column == 3):
            for j in range(ndims):
               n+=1
               array = line.split()
               h[j][i*5:i*5 + 3] = array[1:4]
         elif (column == 4):
            for j in range(ndims):
               n+=1
               array = line.split()
               h[j][i*5:i*5 + 4] = array[1:5]
         n+=1 

  elif hessian == "qchem":
# Read QChem-format molecular Hessian
      lower_triang = ndims * (ndims + 1) // 2
      
      values = []
      found_header = False 
     
      with open(hessfile, "r") as f: 
          for line in f:
             if "Cartesian Force Constants" in line:
                 found_header = True
                 break

          if not found_header:
              raise ValueError("Header 'Cartesian Force Constants' not found in file.")

          # Read numeric tokens until we get lower_triang numbers
          for line in f:
             tokens = line.strip().split()
             for token in tokens:
                values.append(float(token))
                if len(values) == lower_triang:
                  break

      if len(values) != lower_triang:
         raise ValueError(f"Only read {len(values)} of {N} required Hessian constants.")

      lower = np.array(values, dtype=float)
      k = 0
      for i in range(ndims):
         for j in range(i+1):
            h[i, j] = lower[k]
            h[j, i] = h[i, j]
            k += 1

  else:
# Not supported format or a typo 
     raise KeyError(f"The Hessian format '{hessian}' is not supported.")

  sqrtm = np.sqrt(mass)

  # build mass weighted hessian
  h_mw = np.zeros_like(h)

  for idim in range(ndims):
     h_mw[idim, :] = h[idim, :] / sqrtm

  for idim in range(ndims):
     h_mw[:, idim] = h_mw[:, idim] / sqrtm

  # symmetrize mass weighted hessian
  h_mw = 0.5 * (h_mw + h_mw.T)

  # diagonalize mass weighted hessian
  evals, modes = np.linalg.eig(h_mw)


  # sort eigenvectors
  idx = evals.argsort()[::1]
  evals = evals[idx]
  modes = modes[:, idx]

  # delete rotational and translational modes

  evals = evals[6:ndims]
  modes = modes[:,6:ndims]

  if any(val < 0 for val in evals) == True:
    print ("! negative frequency(ies) found. Are you sure you want to continue?  Exiting")
    quit()

  vib_dict = {
     "natoms"   : natoms,
     "ndims"    : ndims,
     "eigvals"  : evals,
     "eigvecs"  : modes,
     "masses"   : mass

  }

  return vib_dict


def electronic(path, sostates, elstates, spin):
    """ Extract SO energies, coefficients, and oscillator strengths """

    with open(path, 'r', errors='ignore') as f:
        lines = f.read().splitlines()

    # -------------------- (1) Energies --------------------
    energies_dict = {}
    try:
        start_idx = next(i for i, ln in enumerate(lines) if "Total energies including SO-coupling" in ln)
    except StopIteration:
        start_idx = 0

    # Matches lines like:
    #   SO-RASSI State 1 Total energy: -XXXXX...
    # Groups:
    #   (1) = state number (integer)
    #   (2) = energy value (float)
    energy_pattern = re.compile(r"SO-RASSI\s+State\s+(\d+)\s+Total energy:\s*([\-0-9\.Ee\+]+)")

    for line in lines[start_idx:]:
        match = energy_pattern.search(line)
        if match:
            state = int(match.group(1))
            if state in sostates:
                energies_dict[state] = float(match.group(2))
        # Stop when we reach the coefficient section
        if "Complex eigenvectors in basis of non-so eigenstates" in line:
            break

    # Fill missing energies with None (will convert to float array; absent entries become np.nan)
    for state in sostates:
        energies_dict.setdefault(state, None)

    # Convert energies to NumPy array in the same order as sostates
    # Note: None will become np.nan when cast to float
    energies = np.array([energies_dict[state] for state in sostates], dtype=float)


    # -------------------- (2) Coefficients --------------------
    # Generate allowed MS projections from spin multiplicity S
    S = spin
    ms_values = []
    cur = -S
    while cur <= S + 1e-8:
        ms_values.append(round(cur, 1))
        cur += 1.0

    # Initialize coefficient containers as nested dicts filled with 0+0j
    coeffs_tmp = {
        f"c{state}": {sfs: {ms: 0.0+0.0j for ms in ms_values} for sfs in elstates}
        for state in sostates
    }

    # --- Regex patterns used for parsing coefficient tables ---

    # Detects header lines like:
    #   "SFS  S  Ms      1       2       3       4"
    # Captures everything after "Ms" so we can extract the SO-state numbers (1,2,3,4)
    state_header_pattern = re.compile(r"SFS\s+S\s+Ms(.*)")

    # Matches data rows like:
    #   "  1   1.50  -1.50   ( 0.010, 0.000) ( -0.002, 0.001) ..."
    # Groups:
    #   (1) = SFS index
    #   (2) = S value (float)
    #   (3) = Ms projection (float)
    #   (4) = rest of the line containing the coefficient pairs
    row_pattern = re.compile(r"^\s*(\d+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+(.*)$")

    # Matches complex-number pairs like "( 0.123, -0.456 )"
    # Groups:
    #   (1) = real part
    #   (2) = imaginary part
    complex_pattern = re.compile(r"\(\s*([\-0-9.]+)\s*,\s*([\-0-9.]+)\s*\)")

    # Parse the coefficient blocks
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Energy (au)"):
            if i + 1 < len(lines):
                header_line = lines[i + 1]
                match = state_header_pattern.search(header_line)
                if match:
                    # Extract state numbers from header (e.g., 1 2 3 4)
                    state_numbers = list(map(int, re.findall(r"\d+", match.group(1))))

                    j = i + 2
                    while j < len(lines):
                        line = lines[j]
                        if not line.strip():
                            j += 1
                            continue

                        # End of this block when encountering new header or section
                        if (
                            line.strip().startswith("Energy (au)")
                            or ("SFS  S     Ms" in line and j != i + 1)
                            or line.strip().startswith("++ ")
                        ):
                            break

                        match_row = row_pattern.match(line)
                        if match_row:
                            sfs = int(match_row.group(1))
                            ms = float(match_row.group(3))
                            comps = complex_pattern.findall(match_row.group(4))  # list of (real, imag)
                            for col_idx, st_no in enumerate(state_numbers):
                                if st_no in sostates and col_idx < len(comps):
                                    re_part, im_part = comps[col_idx]
                                    coeffs_tmp[f"c{st_no}"][sfs][ms] = complex(float(re_part), float(im_part))
                        j += 1
                    i = j
                    continue
        i += 1

    # Convert nested structure into NumPy 1D arrays
    coeffs = {}
    for key, sfs_dict in coeffs_tmp.items():
        arr = []
        for sfs in elstates:
            for ms in ms_values:
                arr.append(sfs_dict[sfs][ms])
        coeffs[key] = np.array(arr, dtype=np.complex128)


    # -------------------- (3) Oscillator strengths  -------------------------
    strengths_map = {}
    in_dipole_section = False   # True only while we're within the Dipole strengths section
    in_table = False            # True only after the "From   To   Osc. strength" header within Dipole
    # Matches transition lines like:
    #   "  1    53    3.25E-04"
    # Groups:
    #   (1) = From-state number
    #   (2) = To-state number
    #   (3) = oscillator strength
    row2_pattern = re.compile(r"^\s*(\d+)\s+(\d+)\s+([0-9.Ee\+\-]+)")

    for line in lines:
        # Enter Dipole section
        if "++ Dipole transition strengths (SO states):" in line:
            in_dipole_section = True
            in_table = False
            continue

        # If we're inside the Dipole section and we hit the start of another "++ ..." section,
        # it means Dipole section ended; stop parsing strengths entirely to avoid velocity overwrite.
        if in_dipole_section and line.startswith("++ "):
            break

        # Detect the Dipole table header "From   To   Osc. strength"
        if in_dipole_section and ("From" in line and "Osc. strength" in line):
            in_table = True
            continue

        # Collect Dipole rows while inside the Dipole table only
        if in_dipole_section and in_table:
            match = row2_pattern.match(line)
            if match:
                from_st = int(match.group(1))
                to_st = int(match.group(2))
                osc = float(match.group(3))
                strengths_map[(from_st, to_st)] = osc
            # allow blank or separator lines to pass; the loop will break when a new "++ " header appears

    # Create a strengths list aligned with provided SO states
    target_to = sostates[-1]
    strengths = [strengths_map.get((state, target_to), 0.0) for state in sostates[:-1]] + [0.0]

    # Return unified structure
    return {"energies": energies, "strengths": strengths, "coefficients": coeffs}


 
def gradients_and_nacs(elstates, natoms, gradspath, nacspath):
  """ Read gradients of and NACs between the CASSCF states """

  grads_and_nacs = np.zeros(shape=(len(elstates)+1, len(elstates)+1, 3*natoms))
  tmp = []


######## Gradients

  if os.path.exists(gradspath) and os.path.isdir(gradspath): 
       for filename in os.listdir(gradspath):
          indx = int(Path(filename).stem[len("gradient"):])
          filepath = os.path.join(gradspath, filename)
          f  = open(filepath)
          lines = f.readlines()
          f.close()
          n = 0
          while n < len(lines):
            line = lines[n]
            if "Molecular gradients" in line:
              break
            n += 1
          n += 7
          for j in range(natoms):
            n += 1
            line = lines[n].split()
            data = [float(num) for num in line[1:4]]
            for x in data: tmp.append(x)
          grads_and_nacs[indx,indx] = np.array(tmp)
          tmp = []
  else:
       raise NotADirectoryError(f"Directory containing gradients not found")


######## NACs

  if os.path.exists(nacspath) and os.path.isdir(nacspath):
    
    for filename in os.listdir(nacspath):
       filepath = os.path.join(nacspath, filename)
       f  = open(filepath)
       lines = f.readlines()
       f.close()
       n = 0

       while n < len(lines):
           line = lines[n]
           if "Lagrangian multipliers" in line:
             line_split = line.split()
             state_indx = [int(t.replace("/", "")) for t in line_split[-2:]]
           if "Energy difference" in line:
              DeltaE = float(line.split()[2])
           if "Total derivative coupling" in line:
             line_split = line.split()
             if "(divided" in line_split:
               mult_factor = float(line_split[-2].replace(")", ""))
             else:
               mult_factor = 1.0
             break
           n += 1

       n += 7
       for k in range(natoms):
          n += 1
          line = lines[n].split()
          data = [float(num) for num in line[1:4]]
          for x in data: tmp.append(x)
       grads_and_nacs[state_indx[0], state_indx[1]] = np.array(tmp)*DeltaE*mult_factor
       tmp = []

  else:
    raise NotADirectoryError(f"Directory containing NACs not found")

  return grads_and_nacs
