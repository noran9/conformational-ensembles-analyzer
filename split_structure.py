import os


def split_sructure(structure_name):
    conformation_id = 0
    with open('./data/' + structure_name + '.pdb', 'r') as f:
        if not os.path.exists('./data/split-structure-' + structure_name):
            os.makedirs('./data/split-structure-' + structure_name)

        conformation_file = open('./data/split-structure-' + structure_name +
                                 '/conformation-' + str(conformation_id) + '.ent', 'w')

        while True:
            line = f.readline()
            if not line:
                break

            conformation_file.write(line)

            if 'ENDMDL' in line:
                conformation_id += 1
                conformation_file = open('./data/split-structure-' + structure_name
                                         + '/conformation-' + str(conformation_id) + '.ent', 'w')
