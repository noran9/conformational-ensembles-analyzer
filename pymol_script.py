from secondary_structure_project import secondary_structure
from Bio import PDB
from Bio.PDB import Selection
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import sys
import os
from matplotlib import cm, colors


def translate_coordinates(num_model, name_PED, translation, structure_i):
    io = PDB.PDBIO()

    for model in structure_i:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(0, atom.get_coord() + translation)

    io.set_structure(structure_i)
    io.save('./translation-files/{}_{}_translated.ent'.format(num_model, name_PED))


if not os.path.exists('./cluster-centers'):
    print('Please compute the graph for a chosen structure first. ')
    sys.exit()

if not os.path.exists('./translation-files'):
    os.makedirs('./translation-files')

with open('./cluster-centers', 'r') as f:

    structure_name_f = f.readline().strip().split(' ')

    k = f.readline().split()

    closest = f.readline().split()

    f.close()

structure_name = structure_name_f[0]
optimal_k = int(k[0])
conformation_id = 0

# Splitting the structure into conformation files
with open('./data/' + structure_name + '.pdb', 'r') as f:
    if not os.path.exists('./data/split-structure-' + structure_name):
        os.makedirs('./data/split-structure-' + structure_name)

    conformation_file = open('./data/split-structure-' + structure_name + '/conformation-' + str(conformation_id), 'w')

    while True:
        line = f.readline()
        if not line:
            break

        conformation_file.write(line)

        if 'ENDMDL' in line:
            conformation_id += 1
            conformation_file = open('./data/split-structure-' + structure_name + '/conformation-' + str(conformation_id), 'w')

ss_matrix = []
center_residue = []

for i in range(len(closest)):
    structure = PDBParser(QUIET=True).get_structure(structure_name,
                                                    "./data/split-structure-{}/conformation-{}".format(
                                                        structure_name, closest[i]))

    ss_matrix.append(secondary_structure(structure))

con = np.ones(len(ss_matrix[0]))

for i in range(len(closest)):
    for j in range(len(ss_matrix[0])):
        if con[j] == 1:
            if ss_matrix[0][j] != ss_matrix[i][j]:
                con[j] = 0

for i in range(len(ss_matrix[0])):
    if con[i] == 1:
        center_residue.append(i)

center_residue = np.array(center_residue)

longest_seq = max(np.split(center_residue, np.where(np.diff(center_residue) != 1)[0] + 1), key=len).tolist()

structure_0 = PDBParser(QUIET=True).get_structure(structure_name,
                                                  "./data/split-structure-{}/conformation-{}".format(structure_name,
                                                                                                     closest[0]))

start_index = [residue.id[1] for residue in structure_0[0]['A'] if residue.id[0] == ' '][0]

for i in range(len(closest)):
    structure_i = PDBParser(QUIET=True).get_structure(structure_name,
                                                      "./data/split-structure-{}/conformation-{}".format(
                                                          structure_name, closest[i]))

    first = structure_0[0]['A'][longest_seq[0] + start_index]["CA"].get_vector()
    second = structure_i[0]['A'][longest_seq[0] + start_index]["CA"].get_vector()
    translation = first - second

    array = translation[:]
    print('The distance between the residue from conformation', closest[0], 'and the residue from conformation', closest[i], 'is: ', array)

    translate_coordinates(closest[i], structure_name, array, structure_i)

distances = []

for i in range(len(closest)):
    row = []
    if i > 0:
        structure_i = PDBParser(QUIET=True).get_structure(structure_name,
                                                          './translation-files/{}_{}_translated.ent'.format(closest[i], structure_name))

        # coords of the model
        res1 = [residue for residue in structure_0[0]['A']]
        res2 = [residue for residue in structure_i[0]['A']]

        for j in range(len(res1)):
            if res1[j].id[0] == " " and res2[j].id[0] == " " and res1[j]["CA"].id[1] == res2[j]["CA"].id[1]:  # Exclude hetero/water residues

                row.append(np.linalg.norm(res1[j]["CA"].get_vector() - res2[j]["CA"].get_vector()))

        distances.append(row)

print("Generating PyMOL image..")
pymol.finish_launching()  # Open Pymol
for i in range(optimal_k):
    cmd.load("./translation-files/{}_{}_translated.ent".format(closest[i], structure_name),
             '{}_{}'.format(closest[i], structure_name))  # Load from file
    cmd.show("cartoon", '{}_{}'.format(closest[i], structure_name))  # Show cartoon

cmd.center()
cmd.zoom(complete=1.2)
norm = colors.Normalize(vmin=np.amin(distances), vmax=np.amax(distances))
for i in range(optimal_k - 1):
    print(i)
    for j, residue in enumerate(Selection.unfold_entities(structure_0[0], "R")):
        rgb = cm.bwr(norm(distances[i][j]))
        print(rgb)
        print(j, residue.id[1], rgb)
        cmd.set_color("col_{}".format(j), list(rgb)[:3])
        cmd.color("col_{}".format(j), "resi {}".format(residue.id[1]))

if not os.path.exists('./images'):
  os.makedirs('./images')

cmd.png("images/pymol_{}.png".format(structure_name), width=2000, height=2000, ray=1)

print("Image is generated in images/pymol_{}.png".format(structure_name))