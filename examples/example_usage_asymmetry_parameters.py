import numpy as np
from matplotlib import pyplot as plt
import sys

# We need to import modules from the relcode_py repository.
# We can add the relcode_py repo to our python path, or we can add it manually in any particular script.
# In this example we do the latter. The path to where the relcode_py repo is located will likely be
# different from what is used in this example!
# If you are using some integrated development environment there is probably a better workflow
# for you for addding the relcode_py repo to your environment.

relcode_py_repo_path = "/home/jsorngard/Mirrors/atomlx04/Repo/relcode_py"

sys.path.append(relcode_py_repo_path)

from fortran_output_analysis.twophotons import TwoPhotons
from fortran_output_analysis.common_utility import kappa_from_l_and_j

#
# Defining paths to the Fortran output data:
#

data_dir = "/home/jsorngard/Mirrors/atomlx04/Repo/relcode_examples/argon_example/"
twophoton_data_dir = data_dir+"second_photon/"
path_abs = twophoton_data_dir + "m_elements_abs_-2_2.dat"
path_emi = twophoton_data_dir + "m_elements_emi_-2_2.dat"

# Create instance of the TwoPhotons class.
# The Fortran data can then be added (i.e loaded from file) to this instance.
# The argument to the constructor is just a string that will be stored in two_photons.name.
# It can be used as a plot title or such.
two_photons = TwoPhotons("Ar sidebands 16/18")

#We read in the XUV photon energies here by specifying the path.
two_photons.add_omega(data_dir + "pert_-2_2/omega.dat")
omega = two_photons.omega_eV  # We can also access it like two_photons.omega_Hartree if we want it in a.u.

# Then we can add matrix elements to the TwoPhotons instance.
# We need to give it the path to the output file, and if it's absorption or emission as "abs" or "emi".
# Then we also give it the hole kappa (ie -2 for p3/2 here).
# Last argument is the principal quantum number n, for the hole, mainly used for labels in plots.
two_photons.add_matrix_elements(path_abs, "abs", -2, 3)  # absorption
two_photons.add_matrix_elements(path_emi, "emi", -2, 3)  # emission

# To add the output from a different hole, in this case 3p1/2 we change the paths and kappa accordingly.
# Note that we overwrite the path variables used previously, this is just a matter of convenience.
path_abs = twophoton_data_dir + "m_elements_abs_1_2.dat"
path_emi = twophoton_data_dir + "m_elements_emi_1_2.dat"

# We can use the helper functions if we don't know the maps (l,j) <-> kappa by heart
kappa_3p1half = kappa_from_l_and_j(l=1, j=1/2)
print("l=1, j=1/2 -> kappa: ", kappa_3p1half)  # Should be 1 for l=1, j=1/2
n = 3
two_photons.add_matrix_elements(path_abs, "abs", kappa_3p1half, n)  # absorption
two_photons.add_matrix_elements(path_emi, "emi", kappa_3p1half, n)  # emission

# Then we can get the asymmetry parameters like so:
b2 = two_photons.get_asymmetry_parameter(2, kappa_3p1half)
b4 = two_photons.get_asymmetry_parameter(4, kappa_3p1half)

# And plot them
plt.plot(omega, b2, label="$\\beta_2$")
plt.plot(omega, b4, label="$\\beta_4$")
plt.legend()
plt.xlabel("XUV photon energy [eV]")
plt.ylabel("Asymmetry parameter value")
plt.show()