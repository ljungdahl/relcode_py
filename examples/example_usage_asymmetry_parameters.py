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
from fortran_output_analysis.common_utility import kappa_from_l_and_j, coulomb_phase

#
# Defining paths to the Fortran output data:
#

data_dir = "/home/jsorngard/Mirrors/atomlx04/Repo/relcode_examples/argon_example/full/"
twophoton_data_dir = data_dir+"second_photon/"
path_abs = twophoton_data_dir + "m_elements_abs_-2_2.dat"
path_emi = twophoton_data_dir + "m_elements_emi_-2_2.dat"

# Some physical parameters of the system
Z = 18
binding_energy_p1half = 15.9 #eV
binding_energy_p3half = 15.7 #eV

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
hole_kappa = kappa_from_l_and_j(1,3/2)
two_photons.add_matrix_elements(path_abs, "abs", hole_kappa, 3)  # absorption
two_photons.add_matrix_elements(path_emi, "emi", hole_kappa, 3)  # emission

# We can calculate the possible final kappas with this function
final_kappas = two_photons.final_kappas(hole_kappa, only_reachable = False)

# Then we can retrieve the value of the coupled matrix elements like so
M_abs = [two_photons.get_coupled_matrix_element(hole_kappa, "abs", kf) for kf in final_kappas]
M_emi = [two_photons.get_coupled_matrix_element(hole_kappa, "emi", kf) for kf in final_kappas]


# The energy of the final electron is
# omega - binding_energy + omega_IR for absorption
# omega - binding_energy - omega_IR for emission
# which means that we must match up the indices of the matrix elements from absorption and emission
# so that the same index corresponds to the same final energy
steps_per_IR_photon = 64
electron_kinetic_energy = (omega - binding_energy_p3half)[steps_per_IR_photon:-steps_per_IR_photon]
M_abs = M_abs[:-2*steps_per_IR_photon]
M_emi = M_emi[2*steps_per_IR_photon:]

#We multiply in the coulomb phase of the final photoelectron as well
coul_phase = coulomb_phase(1, electron_kinetic_energy, Z)
M_abs *= coul_phase
M_emi *= coul_phase

b2_p3half = two_photons.get_asymmetry_parameter(2, hole_kappa, M_abs, M_abs)
b4_p3half = two_photons.get_asymmetry_parameter(4, hole_kappa, M_abs, M_abs)