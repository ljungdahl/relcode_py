########################################################################################################################
#
# This code parses the user input file and generates appropriate input files to be read by the Fortran code.
# Typically the Fortran code expects to read the input files from a subfolder "generated_input/" in the same
# folder as the Fortran executable is called from.
# The file with user input is typically in the same folder as the executable is called from,
# and is named "relcode_input_parameters.txt" - The variables in the .txt-files must correspond to those
#
########################################################################################################################
print("\n")
print("===== START Python log output =====")
print("\n")
# First we handle paths to make sure we can import the modules
import sys
import os

# This is the path where this file is located, it should be relative
# to the input_to_fortran/-subfolder so we can import modules.
relcode_py_path = os.path.dirname(os.path.abspath(__file__))

# This is the folder from which this script is called.
# It is here that the user input file is located.
current_workdir_path = os.path.abspath(os.getcwd())

# Update sys.path with relcode_py_path so we can import local modules
sys.path.append(relcode_py_path)

# Import "local" modules
from input_to_fortran.parse_user_input_file import parse_user_input_file
from input_to_fortran.create_files_for_fortran import generate_atom_parameters_file

input_file_name = current_workdir_path+"/"+"relcode_input_parameters.txt"

# Create directories for Fortran file if they don't exist.
generated_input_path = current_workdir_path+"/generated_input"
if not os.path.exists(generated_input_path):
    try:
        os.mkdir(generated_input_path)
        print("Created directory ", generated_input_path)
    except OSError as error:
        print(error)

#
#
#

parsed_vars_dict = parse_user_input_file(input_file_name)
generate_atom_parameters_file(parsed_vars_dict, current_workdir_path, generated_input_path)


print("\n")
print("===== END Python log output =====")
print("\n")