import numpy as np
import re  # Regular expressions
from itertools import islice  # Slicing when reading lines from Fortran files.
from fortran_output_analysis.constants_and_parameters import g_eV_per_Hartree, g_omega_IR
from fortran_output_analysis.common_utility import l_from_kappa, l_to_str, \
     wigner_eckart_phase, wigner3j_numerical2, j_from_kappa, \
     j_from_kappa_int, IonHole


# ==================================================================================================
#
# ==================================================================================================
class IonisationPath:
    def __init__(self, intermediate_kappa, final_kappa, file_col_idx, col_idx):
        self.kappa_intermediate = intermediate_kappa
        self.l_intermediate = l_from_kappa(intermediate_kappa)
        self.j_intermediate = j_from_kappa(intermediate_kappa)
        self.kappa_final = final_kappa
        self.l_final = l_from_kappa(final_kappa)
        self.j_final = j_from_kappa(final_kappa)
        self.name_intermediate = "" + l_to_str(self.l_intermediate) + ("_{%i/2}" % (j_from_kappa_int(intermediate_kappa)))
        self.name_final = "" + l_to_str(self.l_final) + ("_{%i/2}" % (j_from_kappa_int(final_kappa)))

        # We don't store zero columns from the input file,
        # so the index in the file and the index used to retrieve data will differ.
        self.file_column_index = file_col_idx
        self.column_index = col_idx

    def full_name(self):
        return self.name_intermediate + " to " + self.name_final


# ==================================================================================================
#
# ==================================================================================================
class MatrixElements:
    # This class handles the fortran output data
    # for each channel described by a two-photon (XUV+IR) interaction:
    # hole - intermediate - final
    # It also handles the summation over intermediate states required
    # for a "measurable" description of the channel.
    def __init__(self, path, hole_kappa, hole_n, abs_or_emi):
        self.is_initialised = False
        self.path = path
        self.hole = IonHole(hole_kappa, hole_n)
        self.ionisation_paths = {}
        file = open(path, "r")

        # Sets up what the available ionisation paths are, and what columns in the file it is corresponding to.
        parse_first_line_from_fortran_matrix_element_output_file(file, self.hole, self.ionisation_paths)
        self.number_of_ionisation_paths = len(self.ionisation_paths)

        # We use this list to only pick the non-zero data points after parsing raw data below.
        self.path_col_indices = []
        for ionisation_path in self.ionisation_paths.values():
            self.path_col_indices.append(ionisation_path.file_column_index)

        # These will be filled with parsing the rest of the file.
        self.raw_data_real = np.zeros(1)
        self.raw_data_imag = np.zeros(1)
        # For any energy the matrix element is printed at each breakpoint from the Fortran program.
        # So here we can control what breakpoints we choose.
        self.breakpoint_step = 5

        self.raw_data_real, self.raw_data_imag = \
            parse_matrix_element_raw_data_from_fortran_output_file(file,
                                                                   self.number_of_ionisation_paths,
                                                                   self.path_col_indices,
                                                                   self.breakpoint_step)
        file.close()

        # Errors for this are catched outside this class.
        abs_or_emi_string = "emission"
        if(abs_or_emi == "abs"):
            abs_or_emi_string = "absorption"

        self.name = self.hole.name + " " + abs_or_emi_string + " matrix elements."

        print("Added:")
        print(self.name)
        self.is_initialised = True

    def get_ionisation_path(self, intermediate_kappa, final_kappa):
        if self.is_initialised:
            key_tuple = (self.hole.kappa, intermediate_kappa, final_kappa)
            ionisation_path = self.ionisation_paths[key_tuple]
            column_index = ionisation_path.column_index
            z = self.raw_data_real[:, column_index] + 1j * self.raw_data_imag[:, column_index]
            name = self.hole.name + "$ to $" + ionisation_path.full_name()
            return z, name
        else:
            raise ValueError("MatrixElements not initialised!")

    def get_ionisation_path_summed_over_intermediate(self, final_kappa, mj):

        paths_summed_over_intermediate_states = sum_over_intermediate_states_including_3j_symbols(self.hole.kappa,
                                                                                                  final_kappa,
                                                                                                  self.raw_data_real,
                                                                                                  self.raw_data_imag,
                                                                                                  self.ionisation_paths,
                                                                                                  mj)

        for state_tuple in paths_summed_over_intermediate_states.keys():

            if final_kappa == state_tuple[1]:
                final_label = l_to_str(l_from_kappa(final_kappa)) + ("_{%i/2}" % (j_from_kappa_int(final_kappa)))
                name = self.hole.name + "$  to  $" + final_label
                return paths_summed_over_intermediate_states[state_tuple], name

        # If channel doesn't exist we return None and handle error outside this function.
        return None

    def get_summed_channels(self, mj):
        # This function just gets the hole- and final-kappas present, ie implies the sum over all intermediate.
        # It returns a list of the unique tuples (hole_kappa, final_kappa).
        # It can be used as input to getting all the data for these channels summed over intermediate states.
        hole_kappa = self.hole.kappa
        channels = []
        for path_tuple in self.ionisation_paths.keys():
            final_kappa = path_tuple[2]

            j_hole = j_from_kappa_int(hole_kappa)
            j_final = j_from_kappa_int(final_kappa)
            mjj = int(2 * mj)
            if (np.abs(mjj) > j_hole or np.abs(mjj) > j_final):
                continue

            if (hole_kappa, final_kappa) not in channels:
                channels.append((hole_kappa, final_kappa))

        return channels


# ==================================================================================================
#
# ==================================================================================================
class TwoPhotons:
    def __init__(self, atom_name):
        self.atom_name = atom_name
        self.matrix_elements_abs = {}
        self.matrix_elements_emi = {}
        self.eV_per_Hartree = g_eV_per_Hartree
        self.omega_path = ""
        self.omega_Hartree = 0
        self.omega_eV = 0

    def add_omega(self, path):
        self.omega_path = path
        self.omega_Hartree = np.loadtxt(path)
        self.omega_eV = self.omega_Hartree * self.eV_per_Hartree
        #print(self.omega_eV.shape)

    # Add matrix elements after absorption or emission of an IR photon
    # for a particular hole (from XUV absorption),
    # using the output data from the Fortran program.
    def add_matrix_elements(self, path, abs_or_emi, hole_kappa, hole_n):
        if (abs_or_emi == "abs"):
            if (path.find("abs") == -1):
                raise ValueError("Specified path not compliant with abs/emi convention! Got: path.find('abs') == -1")
            else:
                self.matrix_elements_abs[hole_kappa] = MatrixElements(path, hole_kappa, hole_n, abs_or_emi)

        elif (abs_or_emi == "emi"):
            if (path.find("emi") == -1):
                raise ValueError("Specified path not compliant with abs/emi convention! Got: path.find('emi') == -1")
            else:
                self.matrix_elements_emi[hole_kappa] = MatrixElements(path, hole_kappa, hole_n, abs_or_emi)

        else:
            raise ValueError("Need to specify emission ('emi') or absorption ('abs') when adding matrix elements data!")

        return

    # This is a method to get the data that is summed over intermediate states and over m values.
    def get_matrix_elements_for_hole_to_final(self, hole_kappa, final_kappa, abs_or_emi, mj):

        retval = np.zeros(1)
        if (abs_or_emi == "abs"):
            retval, name = \
                self.matrix_elements_abs[hole_kappa].get_ionisation_path_summed_over_intermediate(final_kappa, mj)

        elif (abs_or_emi == "emi"):
            retval, name = \
                self.matrix_elements_emi[hole_kappa].get_ionisation_path_summed_over_intermediate(final_kappa, mj)
        else:
            raise ValueError("Need to specify emission ('emi') or absorption ('abs') when getting matrix element data!")

        if len(retval) < 2:
            raise ValueError("Couldn't get matrix elements for hole_kappa %i and final_kappa %i, \n"
                             "get_ionisation_path_summed_over_intermediate_and_mj returned None." % (hole_kappa, final_kappa))
        else:
            return retval, name

    def get_hole_to_final_channels(self, hole_kappa, abs_or_emi, mj):
        channels = []
        if (abs_or_emi == "abs"):
            channels = self.matrix_elements_abs[hole_kappa].get_summed_channels(mj)

        elif (abs_or_emi == "emi"):
            channels = self.matrix_elements_emi[hole_kappa].get_summed_channels(mj)
        else:
            raise ValueError("Need to specify emission ('emi') or absorption ('abs') when getting matrix element data!")

        return channels

    def get_matrix_element_for_intermediate_resolved_channel(self, hole_kappa, intermediate_kappa, final_kappa, abs_or_emi):
        retval = np.zeros(1)
        if (abs_or_emi == "abs"):
            retval, name = \
                self.matrix_elements_abs[hole_kappa].get_ionisation_path(intermediate_kappa, final_kappa)

        elif (abs_or_emi == "emi"):
            retval, name = \
                self.matrix_elements_emi[hole_kappa].get_ionisation_path(intermediate_kappa, final_kappa)
        else:
            raise ValueError("Need to specify emission ('emi') or absorption ('abs') when getting matrix element data!")

        if len(retval) < 2:
            raise ValueError("Couldn't get matrix elements for hole_kappa %i and final_kappa %i, \n"
                             "get_ionisation_path_summed_over_intermediate_and_mj returned None." % (
                             hole_kappa, final_kappa))
        else:
            return retval, name








def parse_first_line_from_fortran_matrix_element_output_file(file, in_hole, ionisation_paths):
    # The first line contains information about what is in each column of the output file.
    # They are the kappas for the hole - intermediate - final channels, according to:
    # <hole> <offset_don't_care> <intermediate1> <intermediate2> <intermediate3> <final1> <final2> ... <final9>
    first_line = file.readline().rstrip()
    # split using regex - The kappas are separated by spaces.
    split_first_line = re.split("\s+", first_line)
    # For some reason we get an extra first element that is an empty string.
    # We discard this
    split_first_line = split_first_line[1:]

    # Catch input error here:
    hole_kappa = int(split_first_line[0])
    assert hole_kappa == in_hole.kappa, \
        "Mismatch between hole kappa read from file, and the input to Channels constructor"

    #print(split_first_line)

    intermediate_kappas_str = split_first_line[2:5]
    #print(intermediate_kappas_str)

    final_kappas_str = split_first_line[5:]
    #print(final_kappas_str)

    final_col_index = 0
    raw_data_col_index = 0
    intermediate_stride = 3  # There are three possible final states (including zero for non-channel) for each intermediate.
    intermediate_index = 0
    for intermediate_kappa_str in intermediate_kappas_str:

        intermediate_kappa = int(intermediate_kappa_str)

        i = intermediate_index
        for j in range(3):  # 3 finals per kappa.
            final_index = j + i * intermediate_stride
            final_kappa = int(final_kappas_str[final_index])
            if (intermediate_kappa != 0 and final_kappa != 0):
                ionisation_paths[(hole_kappa, intermediate_kappa, final_kappa)] = \
                    IonisationPath(intermediate_kappa, final_kappa, final_col_index, raw_data_col_index)
                raw_data_col_index += 1

            final_col_index += 1 # note how this is inceremented even though we have a zero column.

        intermediate_index += 1


def parse_matrix_element_raw_data_from_fortran_output_file(file,
                                                           number_of_ionisation_paths,
                                                           non_zero_col_indices,
                                                           breakpoint_step):
    N = 9  # we know we have 9 columns with real/imag pairs
    real_line_np = np.zeros(N, dtype=np.double)
    imag_line_np = np.zeros(N, dtype=np.double)

    real_dynamic = []
    imag_dynamic = []

    # Read rest of lines with actual matrix elements
    for line in islice(file, 1, None, breakpoint_step):  # Skip first line and get every breakpoint_step:th

        line = line.replace(" ", "")  # remove whitespace
        line = line.split(")(")  # split by parentheses
        line[0] = line[0].replace("(", "")  # remove stray parenthesis from first element.
        line[-1] = line[-1].replace(")\n", "")  # remove crap from last element
        # print(line)

        np_index = 0
        for real_imag_pair in line:
            real, imag = real_imag_pair.split(",")
            real_line_np[np_index] = np.double(real)
            imag_line_np[np_index] = np.double(imag)
            np_index += 1

        # print(real_line_np)

        real_dynamic.append(np.copy(real_line_np))
        imag_dynamic.append(np.copy(imag_line_np))

    raw_data_real = np.array(real_dynamic)
    raw_data_imag = np.array(imag_dynamic)

    return raw_data_real[:, non_zero_col_indices], raw_data_imag[:, non_zero_col_indices]


def sum_over_intermediate_states_including_3j_symbols(hole_kappa, final_kappa,
                                                      raw_real, raw_imag,
                                                      ionisation_paths, mj):
    # This function gives back a dictionary containing the two-photon matrix elements (as a function of energy).
    # They are from an initial state (hole) to the final continuum state, summed over all possible intermediate
    # continuum states. We have to provide the mj value to correctly describe the states. If the mj value is not
    # valid for any of the channels we throw an error.
    # Note that this assumes linearly polarised light, so mj will be the same for initial-intermediate-final states.

    paths_summed = {}
    # First we check what has to be summed over.
    # For each unique final state we check what the intermediates are.
    unique_states = []
    all_j_numbers = []

    # NOTE(anton): This is kind of weird code right now, I did it more generally before but at least this
    # checks if the hole and final kappas are states that actually exist.
    for path_tuple in ionisation_paths.keys():
        if (hole_kappa, final_kappa) not in unique_states:
            unique_states.append((hole_kappa, final_kappa))

        for kappa in path_tuple:
            j = j_from_kappa(kappa)
            all_j_numbers.append(j)

    max_mj = np.max(all_j_numbers)
    #print("max mj: ", max_mj)

    N = raw_real.shape[0]
    summed_real = np.zeros((N, len(unique_states)), dtype=np.double)
    summed_imag = np.zeros((N, len(unique_states)), dtype=np.double)
    # For each unique hole and final state we sum over the intermediate j.

    if np.abs(mj) <= max_mj:
        for path_tuple in ionisation_paths.keys():
            ionisation_path = ionisation_paths[path_tuple]
            col_idx = ionisation_path.column_index

            #hole_kappa = path_tuple[0]
            intermediate_kappa = path_tuple[1]
            #final_kappa = path_tuple[2]

            j_hole = j_from_kappa_int(hole_kappa)
            j_intermediate = j_from_kappa_int(intermediate_kappa)
            j_final = j_from_kappa_int(final_kappa)
            mjj = int(2 * mj)
            if (np.abs(mjj) > j_hole or np.abs(mjj) > j_final or np.abs(mjj) > j_intermediate):
                # This will result in a zero 3j symbol for sure, so skip this.
                continue

            #
            unique_state_index = -1
            for unique_state in unique_states:
                unique_state_index += 1
                if (hole_kappa, final_kappa) == unique_state:
                    break

            # (-1)^(j-m) * 3jsymbol
            wigner3j_hole_intermediate = wigner_eckart_phase(intermediate_kappa, mj)*wigner3j_numerical2(j_hole, j_intermediate, mjj)
            wigner3j_intermediate_final = wigner_eckart_phase(final_kappa, mj)*wigner3j_numerical2(j_intermediate, j_final, mjj)

            w3j_total = np.double(wigner3j_hole_intermediate)*np.double(wigner3j_intermediate_final)

            # If the 3j symbol is zero we skip the term, if we somehow didn't catch this earlier.
            if(w3j_total != None and np.abs(w3j_total) > 1e-12):
                summed_real[:, unique_state_index] += w3j_total*raw_real[:, col_idx]
                summed_imag[:, unique_state_index] += w3j_total*raw_imag[:, col_idx]
                #print(path_tuple, ("%i/2" % mjj), w3j_total)

            #print(unique_state_index, path_tuple)

    else:
        # This is an invalid mj provided as input and it won't compute any non-zero 3j-symbol.
        # So we throw an error on this input.
        raise ValueError("Invalid mj-value on input to function. No non-zero 3j-symbols will be computed.")

    unique_state_index = 0
    for unique_state in unique_states:
        z = summed_real[:, unique_state_index] + 1j*summed_imag[:, unique_state_index]
        paths_summed[unique_state] = z
        unique_state_index += 1

    return paths_summed


# ==================================================================================================
#
# ==================================================================================================
def XUV_to_abs_or_emi(x, abs_or_emi):
    xlabel_shift = ""
    if(abs_or_emi == "abs"):
        x = x + g_omega_IR
        xlim = xlim + g_omega_IR
        xlabel_shift = " + IR"
    elif(abs_or_emi == "emi"):
        x = x - g_omega_IR
        xlim = xlim - g_omega_IR
        xlabel_shift = " - IR"
    else:
        raise ValueError("abs_or_emi needs to be 'abs' or 'emi' on entry")

    return x, xlabel_shift


# ==================================================================================================
# ## OLD UNUSED
# ==================================================================================================
# class ParticleState:
#     def __init__(self, kappa):
#         self.kappa = kappa
#         self.l = l_from_kappa(kappa)
#         self.j = j_from_kappa(kappa)
#         self.name = "\epsilon " + l_to_str(self.l) + ("_{%i/2}" % (j_from_kappa_int(kappa)))