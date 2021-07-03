from Bio.PDB import DSSP


# Relative accessible surface area (ASA)
def relative_asa(conformation, structure_name):
    # Load the structure

    dssp = DSSP(conformation, "data/{}.pdb".format(structure_name), dssp="binx/dssp/mkdssp")

    asa = []
    for ss in dssp:
        # The dssp data returned for a single residue is a tuple with the index 3 indicates relative ASA
        asa.append(ss[3])

    return asa
