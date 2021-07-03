import math
from Bio.PDB import PPBuilder


def secondary_structure(conformation):
    # Load the structure

    # Setting the different regions of the Ramachandran plot
    rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                      (-180, 50, 80, 130, 'E', 'blue'),
                      (-100, -180, 100, 60, 'P', 'green'),
                      (-100, 50, 100, 130, 'P', 'green'),
                      (-180, -120, 180, 170, 'H', 'red'),
                      (0, -180, 180, 360, 'L', 'yellow')]

    # Calculate PSI and PHI
    ppb = PPBuilder()  # PolyPeptideBuilder
    rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }
    for chain in conformation:
        for pp in ppb.build_peptides(chain):

            phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
            for i, residue in enumerate(pp):

                # Convert radians to degrees and remove first and last value that are None
                if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                    rama.setdefault(chain.id, [[], [], []])
                    rama[chain.id][0].append(residue)
                    rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
                    rama[chain.id][2].append(math.degrees(phi_psi[i][1]))

    secondary_structure_index = []

    for chain_id in rama:
        for residue, phi, psi in zip(*rama[chain_id]):
            ss_class = None
            for x, y, width, height, ss_c, color in rama_ss_ranges:
                if x <= phi < x + width and y <= psi < y + height:
                    ss_class = ss_c
                    break
            secondary_structure_index.append(ss_class)

    return secondary_structure_index
