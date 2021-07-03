from Bio.PDB import DSSP, PDBParser
from relative_asa import relative_asa
import pandas as pd
import numpy as np


def get_median_solvent_accessibility(conformation, structure_name):
    # Create a dataframe which contains all values of relative ASA for each conformations
    for model_id in range(len(conformation)):
        if model_id == 0:
            asa = relative_asa(conformation[model_id], structure_name)
            df = pd.DataFrame(asa, columns = ['model_{}'.format(model_id)])
        else:
            asa = relative_asa(conformation[model_id], structure_name)
            df['model_{}'.format(model_id)] = asa

    # Convert dataframe into numpy array
    to_np_array = df.to_numpy()

    # Calculate median solvent accessibility
    median_solvent_accessibility = []
    for row in to_np_array:
        median_solvent_accessibility.append(np.median(row))
    return(median_solvent_accessibility)