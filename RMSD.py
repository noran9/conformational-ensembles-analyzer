import numpy as np
from Bio.PDB import Superimposer

def get_average_rmsd(structure):
    # Superimpose all models to the first model, fragment-by-fragment (sliding window)
    super_imposer = Superimposer()
    structure_rmsd_fragments = []  # RMSD, no_models X no_fragments X fragment_size
    window_size = 5
    ref_model = [atom for atom in structure[0].get_atoms() if atom.get_name() == "CA"]  # CA of the first model

    # Iterate all models
    for i, model in enumerate(structure):
        if i > 0:
            model_rmsd = []  # RMSD, no_fragment X fragment_size
            alt_model = [atom for atom in model.get_atoms() if atom.get_name() == "CA"]  # coords of the model

            # Iterate fragments
            for start in range(len(ref_model) - window_size):
                end = start + window_size
                ref_fragment = ref_model[start:end]
                alt_fragment = alt_model[start:end]

                # Calculate rotation/translation matrices
                super_imposer.set_atoms(ref_fragment, alt_fragment)
                # print(super_imposer.rms, super_imposer.rotran)

                # Rotate-translate coordinates
                alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
                alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
                alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]

                # Calculate RMSD
                # https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
                ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
                dist = ref_fragment_coord - alt_fragment_coord
                # rmsd_fragment = np.sqrt(np.sum(dist * dist) / window_size)  # Total RMSD of the fragment. Identical to super_imposer.rms
                rmsd_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment

                model_rmsd.append(rmsd_res)

            structure_rmsd_fragments.append(model_rmsd)

    # Calculate average RMSD per position
    structure_rmsd_fragments = np.array(structure_rmsd_fragments)  # no_models X no_fragments X fragment_size
    # Calculate the RMSD average for each fragments along all models
    structure_rmsd_fragments = np.average(structure_rmsd_fragments, axis=0)  # no_fragments X fragment_size
    # Pad with right zeros to reach the sequence length (no_fragments + fragment_size)
    structure_rmsd_fragments = np.pad(structure_rmsd_fragments, ((0, 0), (0, structure_rmsd_fragments.shape[0])))

    # Roll the fragments one by one (add heading zeros)
    for i, row in enumerate(structure_rmsd_fragments):
        structure_rmsd_fragments[i] = np.roll(row, i)

    # Calculate average along columns of overlapping fragments (average RMSD per residue)
    structure_rmsd_average = np.average(structure_rmsd_fragments, axis=0)
    return structure_rmsd_average
