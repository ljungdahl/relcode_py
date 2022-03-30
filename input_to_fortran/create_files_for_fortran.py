from input_to_fortran.list_of_user_input_variables import get_list_of_user_input_vars

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

def write_string_to_file(file, string):
    file.write(string+"\n")
    return

def write_integer_var_comment_and_value(file, var_str, value):
    comment_str = "# " + var_str
    val_str = "%i" % value
    write_string_to_file(file, comment_str)
    write_string_to_file(file, val_str)
    return

def write_double_var_comment_and_value(file, var_str, value):
    comment_str = "# " + var_str
    val_str = "%.1fd0" % value
    write_string_to_file(file, comment_str)
    write_string_to_file(file, val_str)
    return

def generate_atom_parameters_file(parsed_vars_dict, current_working_dir, generated_input_path):
    global g_list_of_orbital_tuples

    var_list = get_list_of_user_input_vars()

    filename = generated_input_path+"/"+"atom_parameters.input"
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
        tuple = g_list_of_orbital_tuples[i]
        write_str = "%i %i %i %i" % (tuple[0], tuple[1], tuple[2], tuple[3])
        write_string_to_file(file, write_str)


    file.close()
    return


