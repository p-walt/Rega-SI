from rdkit import Chem
import csv
# TODO RA: Make a function? 
# TODO RA: Maybe some comments/docs to explain. 

with open('Whole_molecule_corrected_reagents.csv', 'r') as csv_file:
    compounds_dict = {}
    ts_converged = {}
    molecule_sites = {}
    compound_smiles_dict = {}
    prev_calc_dict = {}
    radicals_dict = {}
    csv_reader = list(csv.reader(csv_file, delimiter=','))[1:]
    temp_dict = {}
    new_mopac_2016 = True
    # print(f'CSV READER IS {csv_reader}')
    none_list = []
    for i, smile_line in enumerate(csv_reader):
        reaction_smiles = smile_line[5]
        print(reaction_smiles)
        smiles_split = reaction_smiles.split('>>')
        reactant_smi = smiles_split[0]
        product_smi = smiles_split[1]
        reactant = Chem.MolFromSmiles(reactant_smi)
        if reactant is None:
            none_list.append((i+2, reaction_smiles))
        print(reactant)
        product = Chem.MolFromSmiles(product_smi)
        if product is None:
            none_list.append((i+2, reaction_smiles))
        print(product)
    print(none_list)
    with open('wrong_smiles5.txt', 'w+') as file:
        for wrong in none_list:
            file.write(f'{wrong}\n')
