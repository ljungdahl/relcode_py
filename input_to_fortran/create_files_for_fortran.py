import os
import numpy as np
from input_to_fortran.list_of_user_input_variables import g_user_input_params_list, \
    g_bool_parameters, g_string_parameters
from input_to_fortran.create_knotpoint_sequence import get_knotpoint_sequence_from_params, write_box_parameters_to_file

g_eV_per_Hartree = 27.211396641308

##################################################################################################
#
##################################################################################################
def write_string_to_file(file, string):
    file.write(string + "\n")
    return


def write_integer_var_comment_and_value(file, var_str, value):
    comment_str = "# " + var_str
    val_str = "%i" % value
    write_string_to_file(file, comment_str)
    write_string_to_file(file, val_str)
    return


def write_double_var_comment_and_value(file, var_str, value):
    comment_str = "# " + var_str
    val_str = "%.5fd0" % value
    write_string_to_file(file, comment_str)
    write_string_to_file(file, val_str)
    return

def create_folder_if_it_doesnt_exist(path):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
            print("Created directory ", path)
        except OSError as error:
            print(error)

def is_folder_empty(path):
    empty = False
    if os.path.isdir(path):
        if not os.listdir(path):
            empty = True

    return empty

##################################################################################################
#
##################################################################################################
def create_atom_parameters_file(parsed_vars_dict, generated_input_path):

    charge_Z = parsed_vars_dict["nuclear_charge_Z"]
    list_of_orbital_tuples = get_atom_orbitals_list(charge_Z)

    var_list = g_user_input_params_list

    filename = generated_input_path + "/" + "atom_parameters.input"
    file = open(filename, "w")

    orbital_counter = 0

    var = "nuclear_charge_Z"
    write_double_var_comment_and_value(file, var, parsed_vars_dict[var])


    num_orbitals_comment = "# num orbitals"
    write_string_to_file(file, num_orbitals_comment)
    num_orbitals = "%i" % len(list_of_orbital_tuples)
    write_string_to_file(file, num_orbitals)

    var ="number_of_holes"
    write_integer_var_comment_and_value(file, var, parsed_vars_dict[var])

    var = "last_kappa"
    write_integer_var_comment_and_value(file, var, parsed_vars_dict[var])

    # Write all the included orbitals in a sequence here
    write_string_to_file(file, "# Orbitals for atom (n, l, 2j, occupation number)")
    for orbital_tuple in list_of_orbital_tuples:
        write_str = "%i %i %i %i" % (orbital_tuple[0], orbital_tuple[1], orbital_tuple[2], orbital_tuple[3])
        write_string_to_file(file, write_str)

    file.close()

    #print("Wrote to %s" % filename)
    return


##################################################################################################
#
##################################################################################################
def create_run_parameters_file(parsed_vars_dict, generated_input_path):
    filename = generated_input_path + "/" + "run_parameters.input"
    file = open(filename, "w")

    is_any_param_true = False
    for param in g_bool_parameters:
        bool_val = parsed_vars_dict[param]
        value = 0
        if bool_val:
            value = 1
            is_any_param_true = True

        comment_str = param
        write_integer_var_comment_and_value(file, comment_str, value)

    if not is_any_param_true:
        print("\nWARNING! Only running ground state calculation!\n")

    if parsed_vars_dict["run_two_photons"] and not parsed_vars_dict["run_one_photon"]:
        print("\nERROR! run_two_photons set to true while run_one_photon set to false.")
        print("two photons require one photon data.")
        raise Exception("Can't run two photons without running first photon calculation.")

    if parsed_vars_dict["run_diagonalise_CIS"] and \
            (parsed_vars_dict["run_one_photon"] or parsed_vars_dict["run_two_photons"]):
        print("\nWARNING! run_diagonalise_CIS set to true will skip one- and two-photon calculations.\n")

    file.close()
    #print("Wrote to %s" % filename)
    return


##################################################################################################
#
##################################################################################################
def create_file_io_parameters_file(parsed_vars_dict, current_workdir, generated_input_path):
    filename = generated_input_path + "/" + "file_io_parameters.input"

    default_dir_path = current_workdir + "/output"
    # Create folders if necessary
    write_path = parsed_vars_dict["path_to_output_folder"]
    if write_path != "default":
        create_folder_if_it_doesnt_exist(write_path)
    else:
        create_folder_if_it_doesnt_exist(default_dir_path)

    read_path = parsed_vars_dict["path_to_previous_output"]

    if len(read_path) > 1:
        if read_path == "default":
            read_path = default_dir_path

        create_folder_if_it_doesnt_exist(read_path)
        # If the folder is empty and we still try to read from it, things will break.
        # So if the folder is given by user but nothing is in it, we will force this option to zero
        # and give a warning.
        if is_folder_empty(read_path):
            print("WARNING! Empty read folder supplied. Ignoring reading from previous calculation this run.")
            parsed_vars_dict["path_to_previous_output"] = "0"

    # Check that exp en file actually exists if it's supposed to be used.
    exp_en_file = parsed_vars_dict["path_to_experimental_energies"]
    if exp_en_file != "0":
        file_exists = os.path.exists(exp_en_file)
        if not file_exists:
            print("FATAL ERROR: Experimental energies file %s doesn't exist." % exp_en_file)
            raise Exception("Invalid experimental energies file.")



    file = open(filename, "w")
    for param in g_string_parameters:
        val_str = parsed_vars_dict[param]
        if val_str == "default":
            val_str = default_dir_path
            create_folder_if_it_doesnt_exist(default_dir_path)

        comment_str = "# %s" % param
        write_string_to_file(file, comment_str)
        write_string_to_file(file, val_str)

    file.close()



    #print("Wrote to %s" % filename)
    return


##################################################################################################
#
##################################################################################################
def create_knotpoint_sequence_and_box_parameters_file(parsed_vars_dict, generated_input_path):
    print("\nGenerating knotpoint sequence:")
    knotsequence, start_imag_coord = get_knotpoint_sequence_from_params(parsed_vars_dict)
    print("\n")
    knotsequence_filename = generated_input_path + "/" + "knotpoint_sequence.dat"
    np.savetxt(knotsequence_filename, knotsequence, delimiter="    ", fmt='%1.13e')
    #print("Wrote to %s" % knotsequence_filename)

    file_params_filename = generated_input_path + "/" + "box_parameters.input"
    file_params = open(file_params_filename, "w")
    write_box_parameters_to_file(file_params, parsed_vars_dict, start_imag_coord)
    file_params.close()
    #print("Wrote to %s" % file_params_filename)

    return


##################################################################################################
#
##################################################################################################
def create_photon_sequence_and_parameters_file(parsed_vars_dict, generated_input_path):
    # We compute the step size as delta_omega = omega_IR/fraction.
    # Note that we translate all input values to atomic units since this is what is appropriate
    # for the Fortran computations.
    fraction = parsed_vars_dict["first_photon_step_fraction"]
    omega_IR_eV = parsed_vars_dict["second_photon_energy"]
    omega_IR_au = omega_IR_eV/g_eV_per_Hartree

    delta_omega = omega_IR_au/fraction

    start_omega_eV = parsed_vars_dict["first_photon_energy_start"]
    end_omega_eV = parsed_vars_dict["first_photon_energy_end"]
    in_start_omega_au = start_omega_eV/g_eV_per_Hartree
    in_end_omega_au = end_omega_eV/g_eV_per_Hartree

    # Compute start and end points adhering to step size delta_omega
    tmp_start = 0.0
    count = 0
    while tmp_start < in_start_omega_au:
        tmp_start += delta_omega
        count += 1

    # Go back one step for starting energy to assure we start slightly before selected energy.
    start_energy_au = tmp_start - delta_omega
    num_start_steps = count - 1  # Withdraw one count since we're taking one step back

    tmp_end = 0.0
    count = 0
    while tmp_end < in_end_omega_au:
        tmp_end += delta_omega
        count += 1

    # We will end up in one step above the chosen energy (or exactly hitting it)
    num_end_points = count

    total_photon_points = num_end_points-num_start_steps

    photon_range_au = np.zeros(total_photon_points)
    for i in range(total_photon_points):
        photon_range_au[i] = start_energy_au + i*delta_omega
        #print("%1.13e" % photon_range_au[i])

    photon_range_filename = generated_input_path + "/" + "photon_range.dat"
    np.savetxt(photon_range_filename, photon_range_au, fmt='%1.13e')
    #print("Wrote to %s" % photon_range_filename)


    photon_params_filename = generated_input_path + "/" + "photon_parameters.input"
    file = open(photon_params_filename, "w")

    num_photons = len(photon_range_au)
    write_integer_var_comment_and_value(file, "number of photons", num_photons)
    write_integer_var_comment_and_value(file, "fraction of omega_IR for step size", fraction)
    write_double_var_comment_and_value(file, "second_photon_energy (eV)", omega_IR_eV)

    file.close()
    #print("Wrote to %s" % photon_params_filename)

    print("\n")
    print("For first photon interval between %.3f eV and %.3f eV, and a step size of omega_IR/%i,"
          " we have %i points." % (start_omega_eV, end_omega_eV, fraction, num_photons))
    print("omega_IR = %.3f" % omega_IR_eV)
    print("\n")

def create_generation_complete_file_for_fortran_validation(generated_input_path):
    filename = generated_input_path + "/" + "generation_complete.input"

    file = open(filename, "w")
    write_string_to_file(file, "# generation complete - fortran checks if this file exists")

    return


def remove_previous_generation_complete_file(generated_input_path):
    filename = generated_input_path + "/" + "generation_complete.input"
    file_exists = os.path.exists(filename)
    if(file_exists):
        os.remove(filename)

    return


def get_atom_orbitals_list(charge_Z):
    # The Fortran program uses a list of orbitals in a (n, l, 2j, occupation_number) format.
    # This function provides that kind of list to the generator script.
    # You can extend this whatever way necessary, for example when calculation some ion.
    # The Fortran program calculates the long range charge as
    # lr_charge = charge_Z-sum(occupation_numbers)+1.0
    # So specifying Z and orbitals with proper occupation numbers should suffice for running ions.

    # He = 1s^2
    helium_orbitals = [
        (1, 0, 1, 2),
    ]

    # Ne = [He] 2s^2 2p^6
    neon_orbitals = helium_orbitals + [
        (2, 0, 1, 2),
        (2, 1, 1, 2),
        (2, 1, 3, 4),
    ]

    # Ar = [Ne] 3s^2 3p^6
    argon_orbitals = neon_orbitals + [
        (3, 0, 1, 2),
        (3, 1, 1, 2),
        (3, 1, 3, 4)
    ]

    # Kr = [Ar] 3d^10 4s^2 4p^6
    krypton_orbitals = argon_orbitals + [
        (3, 2, 3, 4),
        (3, 2, 5, 6),
        (4, 0, 1, 2),
        (4, 1, 1, 2),
        (4, 1, 3, 4)
    ]

    # Xe = [Kr] 4d^10 5s^2 5p^6
    xenon_orbitals = krypton_orbitals + [
        (4, 2, 3, 4),
        (4, 2, 5, 6),
        (5, 0, 1, 2),
        (5, 1, 1, 2),
        (5, 1, 3, 4)
    ]

    # Rn = [Xe] 4f^14 5d^10 6s^2 6p^6
    # Since we now ad 4f before n=5 orbitals we
    # just add first two rows from xenon orbitals first.
    radon_orbitals = krypton_orbitals + [
        (4, 2, 3, 4),
        (4, 2, 5, 6),
        (4, 3, 5, 6),
        (4, 3, 7, 8),
        (5, 0, 1, 2),
        (5, 1, 1, 2),
        (5, 1, 3, 4),
        (5, 2, 3, 4),
        (5, 2, 5, 6),
        (6, 0, 1, 2),
        (6, 1, 1, 2),
        (6, 1, 3, 4)
    ]

    if charge_Z == 2.0:
        return helium_orbitals
    elif charge_Z == 10.0:
        return neon_orbitals
    elif charge_Z == 18.0:
        return argon_orbitals
    elif charge_Z == 36.0:
        return krypton_orbitals
    elif charge_Z == 54.0:
        return xenon_orbitals
    elif charge_Z == 86.0:
        return radon_orbitals
    else:
        print("\nFATAL ERROR: Nuclear charge Z = %.1f not supported!" % charge_Z)
        print("You can add your own set of orbitals to get_atom_orbitals_list() in file %s \n"
              % __file__)
        raise Exception("FATAL ERROR: Nuclear charge Z = %.1f not supported!" % charge_Z,
        "You can add your own set of orbitals to get_atom_orbitals_list() in file %s"
        % __file__)


