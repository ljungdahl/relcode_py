import numpy as np
from fortran_output_analysis.constants_and_parameters import g_eV_per_Hartree
from fortran_output_analysis.common_utility import l_from_kappa, l_to_str, \
    wigner_eckart_phase, wigner3j_numerical2, j_from_kappa, \
    j_from_kappa_int, IonHole, load_raw_data, convert_rate_to_cross_section


# ==================================================================================================
#
# ==================================================================================================
class FinalState:
    def __init__(self, kappa, pcur_col_idx):
        self.kappa = kappa
        self.l = l_from_kappa(kappa)
        self.j = j_from_kappa(kappa)
        self.name = l_to_str(self.l) + ("_{%i/2}" % (j_from_kappa_int(kappa)))
        self.pcur_column_index = pcur_col_idx


# ==================================================================================================
#
# ==================================================================================================
class OnePhoton:
    def __init__(self, atom_name):
        self.name = atom_name
        self.channels = {}
        self.num_channels = 0

    def add_channels_for_hole(self, path_to_pcur_all, hole_kappa, n_qn):
        self.channels[hole_kappa] = Channels(path_to_pcur_all, hole_kappa, n_qn)
        self.num_channels += 1

    def get_channel_labels_for_hole(self, hole_kappa):
        channel_labels = []
        channel = self.channels[hole_kappa]
        hole_name = channel.hole.name
        for final_state_key in channel.final_states.keys():
            final_state = channel.final_states[final_state_key]
            channel_labels.append(hole_name + " to " + final_state.name)

        return channel_labels

    def get_omega_eV(self):
        assert (len(self.channels) > 0)
        first_kappa = list(self.channels.keys())[0]
        first_channel = self.channels[first_kappa]
        omega_eV = first_channel.raw_data[:, 0] * g_eV_per_Hartree  # omega energies in atomic units in file.
        return omega_eV

    def get_omega_Hartree(self):
        assert (len(self.channels) > 0)
        first_kappa = list(self.channels.keys())[0]
        first_channel = self.channels[first_kappa]
        omega_eV = first_channel.raw_data[:, 0] # omega energies in atomic units in file.
        return omega_eV

    def get_partial_cross_section(self, hole_kappa, final_kappa, divide_omega=True):
        # Depending on conventions when creating the dipole elements in the Fortran program we might
        # have to divide or multiply by the photon energy (omega) when calculating cross sections.
        # Usually it is correct to divide by omega, and that is default behaviour of this function.
        hole = self.channels[hole_kappa]
        rate = hole.get_rate_for_channel(final_kappa)
        omega = self.get_omega_Hartree()
        cross_section = convert_rate_to_cross_section(rate, omega, divide_omega)
        return cross_section

    def get_total_cross_section(self, divide_omega=True):
        # Use omega to initialise total cs array:
        omega = self.get_omega_Hartree()
        N = len(omega)
        total_cs = np.zeros(N)
        for channel_key in self.channels.keys():
            channel = self.channels[channel_key]
            hole_kappa = channel.hole.kappa
            for final_key in channel.final_states.keys():
                final_state = channel.final_states[final_key]
                final_kappa = final_state.kappa
                cs = self.get_partial_cross_section(hole_kappa, final_kappa, divide_omega)
                total_cs += cs

        return total_cs




# ==================================================================================================
#
# ==================================================================================================
class Channels:
    def __init__(self, path_to_pcur, hole_kappa, n_qn):
        self.path_to_pcur = path_to_pcur
        self.hole = IonHole(hole_kappa, n_qn)
        self.final_states = {}
        self.raw_data = load_raw_data(path_to_pcur)
        self.add_final_states()

    def add_final_states(self):
        kappa_hole = self.hole.kappa
        # One can convince oneself that the following is true for a given hole_kappa.
        possible_final_kappas = np.array([-kappa_hole, kappa_hole+1, -(-kappa_hole+1)])
        # It is possible that one of the final kappas are zero, so we need to handle this.
        # NOTE(anton): The pcur-files have three columns, one for each possible final kappa.
        # If there is no possibility for one of them the column is zero, and I
        # think the convention is that the zero column is the left-most (lowest index) then.
        # So if the kappas are sorted by ascending absolute value we should get this, since
        # if kappa = 0 the channel is closed.
        sort_kappas_idx = np.argsort(np.abs(possible_final_kappas))
        possible_final_kappas = possible_final_kappas[sort_kappas_idx]
        # This is for getting the data from the pcur files. The first column is the photon energy.
        pcur_column_index = 1
        for kappa in possible_final_kappas:
            if (kappa != 0):
                self.final_states[kappa] = FinalState(kappa, pcur_column_index)

            pcur_column_index += 1

    def get_rate_for_channel(self, final_kappa):
        state = self.final_states[final_kappa]
        column_index = state.pcur_column_index
        rate = self.raw_data[:, column_index]
        return rate


# ==================================================================================================
#
# ==================================================================================================
class ChannelsOld:
    def __init__(self, atom_name):
        self.name = atom_name
        self.ion_holes = {}
        self.eV_per_Hartree = 27.211396641307999
        self.num_holes = 0

    def add_ion_hole(self, path, kappa, n_qn):
        self.ion_holes[kappa] = IonHole(path, kappa, n_qn)
        self.num_holes += 1

    def get_omega_eV(self):
        assert (len(self.ion_holes) > 0)
        first_kappa = list(self.ion_holes.keys())[0]
        first_hole = self.ion_holes[first_kappa]
        omega_eV = first_hole.raw_data[:, 0] * self.eV_per_Hartree  # omega energies in atomic units in file.
        return omega_eV

    def print_channels(self):
        for hole_kappa in self.ion_holes:
            hole = self.ion_holes[hole_kappa]
            for final_state_kappa in hole.final_states:
                final_state = hole.final_states[final_state_kappa]
                channel_name = hole.name + " -> " + final_state.name
                print(channel_name)

    def get_channel_label_lj(self, l_hole, j_hole, l_final, j_final):
        hole_kappa = kappa_from_l_and_j(l_hole, j_hole)
        final_kappa = kappa_from_l_and_j(l_final, j_final)
        return get_channel_label(hole_kappa, final_kappa)

    def get_channel_label_kappa(self, hole_kappa, final_kappa):
        hole = self.ion_holes[hole_kappa]
        final_state = hole.final_states[final_kappa]
        return hole.name + " \\rightarrow " + final_state.name

    def get_channel_label_for_filename(self, hole_kappa, final_kappa):
        hole = self.ion_holes[hole_kappa]
        final_state = hole.final_states[final_kappa]
        l_hole = hole.l
        n_hole = hole.n
        l_final = final_state.l
        hole_str = str(n_hole) + l_to_str(l_hole) + ("%i" % (j_from_kappa_int(hole_kappa)))
        final_str = l_to_str(l_final) + ("%i" % (j_from_kappa_int(final_kappa)))
        return hole_str + "_to_" + final_str

    def get_cross_section_for_channel(self, hole_kappa, final_kappa):
        hole = self.ion_holes[hole_kappa]
        final_state = hole.final_states[final_kappa]
        omega_au = hole.raw_data[:, 0]
        rate = hole.raw_data[:, final_state.pcur_column_index]
        cs = convert_rate_to_cross_section(rate, omega_au, divide=True)
        return cs

    def compute_wigner3j_for_channel(self, hole_kappa, final_kappa, mj):
        mjj = int(2 * mj)
        hole = self.ion_holes[hole_kappa]
        final = hole.final_states[final_kappa]
        j_hole = j_from_kappa_int(hole_kappa)
        j_final = j_from_kappa_int(final_kappa)
        if (mjj > j_hole or mjj > j_final):
            return
        # print("j_hole, j_final, mj: %i/2 %i/2 %i/2" % (j_hole, j_final, mjj))
        K = 1
        q = 0
        w3j = wigner_3j(j_final / 2, K, j_hole / 2, -mjj / 2, q, mjj / 2)
        # print(w3j, sympy_to_num(w3j))
        return sympy_to_num(w3j)

    def get_wigner_eckart_phase(self, final_kappa, mj):
        return np.power(-1.0, (j_from_kappa(final_kappa) - mj))

    def get_total_cross_section(self):
        assert (self.num_holes > 0)
        omega = self.get_omega_eV()
        tot_cs = np.zeros(len(omega))
        for hole_kappa in self.ion_holes.keys():
            hole = self.ion_holes[hole_kappa]
            for final_kappa in hole.final_states.keys():
                part_cs = self.get_cross_section_for_channel(hole_kappa, final_kappa)
                tot_cs += part_cs

        return tot_cs
