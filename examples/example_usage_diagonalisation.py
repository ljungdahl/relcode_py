import numpy as np
from matplotlib import pyplot as plt
import sys
# We need to import modules from the relcode_py repository.
# We can add the relcode_py repo to our python path, or we can add it manually in any particular script.
# In this example we do the latter. The path to where the relcode_py repo is located will likely be
# different from what is used in this example!
# If you are using some integrated development environment there is probably a better workflow
# for you for addding the relcode_py repo to your environment.

relcode_py_repo_path = "/home/anton/lx04/relcode_py/"

sys.path.append(relcode_py_repo_path)

from fortran_output_analysis.diagonalisation import Diagonalisation, compute_diag_cross_section
from fortran_output_analysis.constants_and_parameters import g_eV_per_Hartree

# NOTE(anton): The Diagonalisation-class is not so fleshed out as OnePhoton or TwoPhotons, ie it basically only loads
# in the raw data and parses the channel indices file.
# Initialise the container for the diagCIS-data from Fortran:
neon_diag = Diagonalisation("neon")

data_dir = "/home/anton/lx04/relcode_examples/neon_example/output/"

# Read Fortran output data from file.
# NOTE(anton): For plotting resonances in absorption cross-section the channel data is not used.
# But it's good to load this anyway since we can reason about what the Fortran program is outputting using the
# channel data. Right now the script fails if we haven't loaded channel data!
# Loading channel data entails some parsing of the "diag_channel_indices.dat"-file, and there might be errors here...
neon_diag.load_channel_data(data_dir+"diag_channel_indices.dat")

# For computing the absorption cross section we do need the eigenvalues and the matrix elements.
neon_diag.load_eigenvalues(data_dir+"diag_eigenvalues.dat")
neon_diag.load_matrix_elements(data_dir+"diag_matrix_elements.dat")

# NOTE(anton): Define some energy range and a zero-initialised cross_section array.
num_omegas = 3000
omega_range_ev = np.linspace(30.0, 70.0, num_omegas)
cross_section = np.zeros(len(omega_range_ev))

# Compute cross section for each omega using the output from diagonalisation.
for i in range(len(cross_section)):
    omega = omega_range_ev[i]
    # Now we can just pass the container for all the data, in this case neon.
    cross_section[i] = compute_diag_cross_section(neon_diag, omega) # TODO(anton): add this as a method to the Diag class.

# Example plot of absorption cross section with a line marking 2s1/2 threshold
plt.figure(1)
x = omega_range_ev
y = cross_section
plt.plot(x, y, linewidth=2)
threshold_2s = 1.93583100454*g_eV_per_Hartree
plt.axvline(threshold_2s, linewidth=2, color='black', linestyle='dashed')
plt.title("Absorption cross section for %s near $s_{1/2}$ threshold" % (neon_diag.name))
plt.xlabel("Photon energy [eV]")
plt.ylabel("Cross section [Mbarn]")
plt.grid(linestyle="dotted")
plt.show()