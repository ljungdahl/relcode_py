# This is a hardcoded list of strings corresponding to the parameters in the input file.
# If this is not a 1-to-1 correspondence things will probably not work.
# If additional input parameters are added to the input file they should also be added here.
# Then they also need to be handled in the proper way after parsing, of course.
def get_list_of_user_input_vars():
    user_input_params_list = [
        "nuclear_charge_Z",
        "highest_occupied_orbital",
        "number_of_holes",
        "last_kappa"
    ]
    return user_input_params_list

