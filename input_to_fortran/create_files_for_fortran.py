import numpy as np
from input_to_fortran.list_of_user_input_variables import g_user_input_params_list, \
    g_bool_parameters, g_string_parameters
from input_to_fortran.create_knotpoint_sequence import get_knotpoint_sequence_from_params, write_box_parameters_to_file

g_eV_per_Hartree = 27.211396641308

# This is a list of possible orbitals for an atomic system
# The tuples represent (n, l, 2j, occupation_number)
g_list_of_orbital_tuples = [
    (1, 0, 1, 2),
    (2, 0, 1, 2),
    (2, 1, 1, 2),
    (2, 1, 3, 4),
    (3, 0, 1, 2),
    (3, 1, 1, 2),
    (3, 1, 3, 4),
    (3, 2, 3, 4),
    (3, 2, 5, 6),
    (4, 0, 1, 2),
    (4, 1, 1, 2),
    (4, 1, 3, 4),
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


##################################################################################################
#
##################################################################################################
def create_atom_parameters_file(parsed_vars_dict, generated_input_path):
    global g_list_of_orbital_tuples

    var_list = g_user_input_params_list

    filename = generated_input_path + "/" + "atom_parameters.input"
    file = open(filename, "w")

    orbital_counter = 0
    # We loop through the list since it should be in the same order as the user input .txt-file.
    for var in var_list:
        if var == "nuclear_charge_Z":
            write_double_var_comment_and_value(file, var, parsed_vars_dict[var])

        elif var == "highest_occupied_orbital":
            value = parsed_vars_dict[var]

            for orbital_tuple in g_list_of_orbital_tuples:
                orbital_counter += 1
                if orbital_tuple == value:
                    break
            num_orbitals_comment = "# num orbitals"
            write_string_to_file(file, num_orbitals_comment)
            num_orbitals = "%i" % orbital_counter
            write_string_to_file(file, num_orbitals)

        elif var == "number_of_holes":
            write_integer_var_comment_and_value(file, var, parsed_vars_dict[var])

        elif var == "last_kappa":
            write_integer_var_comment_and_value(file, var, parsed_vars_dict[var])

    # Write all the included orbitals in a sequence here
    write_string_to_file(file, "# Orbitals for atom (n, l, 2j, occupation number)")
    for i in range(orbital_counter):
        orbital_tuple = g_list_of_orbital_tuples[i]
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

    file = open(filename, "w")
    for param in g_string_parameters:
        val_str = parsed_vars_dict[param]
        if val_str == "default":
            val_str = current_workdir + "/output/"
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





