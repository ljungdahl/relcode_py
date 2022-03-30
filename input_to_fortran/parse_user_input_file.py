import numpy as np
from itertools import islice  # Slicing when reading lines from Fortran files.

from input_to_fortran.list_of_user_input_variables import get_list_of_user_input_vars

g_parsed_variables_dict = {}

g_input_file_path = ""
g_line_counter = 0

def raise_error_if_parse_fail(line):
    global g_input_file_path
    global g_line_counter
    print("\n")
    print("ERROR:")
    print("Error parsing ", line, " to key value pair.")
    print("From line %i in file %s" % (g_line_counter, g_input_file_path))
    print("Make sure that user input parameters are properly set and that the variable name exists"
          " in the relcode_py repository.")
    print("\n")
    raise Exception("Couldn't parse key value pair from input file.")

def debug_line_parse(line, key, value):
    print("from line:")
    print(line)
    print("parsed key value pair:", key, value)
    print("type(key):", type(key), "type(value):", type(value))
    print("\n")


def parse_string_to_key_value_pair(line):
    # This is the function that decides how
    key = ""
    value = ""
    key_val_split = line.split("=")
    key_str = key_val_split[0]
    val_str = key_val_split[-1]

    if key_str == "" or val_str == "" or key_str not in get_list_of_user_input_vars():
        raise_error_if_parse_fail(line)

    # Key is just the string of the variable name.
    key = key_str

    # Value can be some different datatypes which we handle here.
    # Most will be integers.
    if key_str == "nuclear_charge_Z":
        value = float(val_str)
    elif key_str == "highest_occupied_orbital":
        val_str = val_str.strip(")")
        val_str = val_str.strip("(")
        val_str = val_str.split(",")
        value = tuple(map(int, val_str))
    else:
        value = int(val_str)

    return key, value


def parse_user_input_file(file_path):
    global g_input_file_path
    global g_line_counter
    global g_parsed_variables_dict
    g_input_file_path = file_path

    print("Parsing user input file %s " % (file_path))

    file = open(file_path, "r")
    for line in islice(file, 1, None):
        g_line_counter += 1
        # Remove trailing newline chars
        line[-1].replace("\n","")
        # We skip the line if it starts with a hash (a comment),
        # or if it's empty
        if line[0] == "#" or len(line) == 0 or not line.strip():
            continue

        # remove whitespace
        line = line.replace(" ", "")

        # Split by comments by just taking the first element in the split string
        line = line.split("#")[0]

        # If we now have a zero length string this line was a comment only line, and we skip.
        if(len(line) == 0):
            continue

        # Else we proceed by saving it as a key value pair.
        key, value = parse_string_to_key_value_pair(line)

        g_parsed_variables_dict[key] = value
        #debug_line_parse(line, key, value)

    file.close()
    for key, value in g_parsed_variables_dict.items():
        print(key,": ", value)

    return g_parsed_variables_dict
