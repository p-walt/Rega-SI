import shutil
import store_data as sd
from pathlib import Path
import re
import os
import json
from subprocess import DEVNULL, STDOUT, check_call
import subprocess
import csv
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from itertools import chain


def read_mopac_output(in_file: Path, out_file: Path, directory: Path) -> int:
    """
    Reads a MOPAC output file, processes specified sections, and writes processed coordinates
    to an output XYZ file with additional side processing such as reversing file lines.

    This function extracts atomic coordinates from a MOPAC output file, removes unnecessary
    data, and rearranges it for creating a formatted XYZ file. It also performs some basic
    error detection in the MOPAC output file. A temporary file is used during this process,
    which is deleted after usage.

    :param in_file: The input file path containing MOPAC output.
    :type in_file: Path
    :param out_file: The output file path where the processed XYZ file will be written.
    :type out_file: Path
    :param directory: The directory path where temporary files will be created and processed.
    :type directory: Path
    :return: Status code indicating the success (1) or failure (2) of the operation.
    :rtype: int
    """
    # Read am1 output line by line but reverse order and put into new temporary file out_rev.txt.
    with open(in_file) as rf, open(Path(f'{directory}/out_rev.txt'), 'w') as wf:
        for line in reversed(rf.readlines()):
            wf.write(line)

    # Create empty list to put cartesian coordinates from output file into
    a = []
    # Scan the reversed file for the keyword EIGENVALUES and collect the following lines until CARTESIAN COORDINATES
    # is reached and add them to the list a.
    with open(Path(f'{directory}/out_rev.txt')) as file:
        for line in file:
            lines = line.strip('\n')
            if 'CALCULATION IS TERMINATED' in lines:
                return 2
            elif 'Too many variables. By definition, at least one force constant is exactly zero' in lines:
                return 2
            elif 'EXCESS NUMBER OF OPTIMIZATION CYCLES' in lines:
                return 2
            elif 'TS FAILED TO LOCATE TRANSITION STATE' in lines:
                return 2
            elif 'UNABLE TO ACHIEVE SELF-CONSISTENCE' in lines:
                return 2
            elif 'Error and normal termination messages reported in this calculation' in lines:
                return 2
            elif 'Empirical Formula' in lines:
                # collect block-related lines
                while True:
                    try:
                        lines = next(file)
                    except StopIteration:
                        # there is no lines left
                        break
                    if 'CARTESIAN COORDINATES' in lines:
                        # we've reached the end of block
                        break
                    a.append(lines)
                # stop iterating over file
                break
    # Reverse the list back to the right way round.
    a.reverse()
    # Remove the lines CARTESIAN COORDINATES and EIGENVALUES and the blank lines before and after the coordinates of
    # interest.
    trimmed = a[1:-2]
    coordinates = []
    # Remove the first column from the coordinates output of MOPAC which gives the atom number which is not needed.
    # This leaves ATOM LABEL, X, Y Z
    for atom in trimmed:
        strip = atom.strip('\n')
        line_elements = re.split(' +', str(strip))
        coordinate = line_elements[1:]
        if not coordinate:
            break
        else:
            del coordinate[0]

        coordinate1 = '\t'.join(map(str, coordinate))
        coordinates.append(coordinate1)
    system_size = len(trimmed)
    # Write AM1 optimised coordinates to an xyz file
    with open(out_file, 'w+') as xyz:
        xyz.write(str(system_size) + '\n\n')
        for i in coordinates:
            xyz.write(str(i) + '\n')
    os.remove(Path(f'{directory}/out_rev.txt'))
    return 1


directory = Path(f'{str(Path("/")).join(str(Path.cwd()).split(str(Path("/")))[:-1])}/Calculation_data')
print(directory)
with open(str(Path(f'{directory}/all_calcs.txt')), 'r') as all_calcs:
    alines1 = all_calcs.readlines()
    aline2 = []
    dir_list = []
    dir_dict = {}
    results_dict = {}

    for fi in os.listdir(directory):
        d = os.path.join(directory, fi)
        if os.path.isdir(d):
            dir_list.append(d)
    for aline1 in alines1:
        st = aline1.strip('\n')
        aline2.append(st)
    with open(str(Path(f'{directory}/sites_complete.txt')), 'r') as comp_calcs:
        clines1 = comp_calcs.readlines()
        cline2 = []
        cline3 = []
        for cline1 in clines1:
            ct = cline1.strip('\n')
            ct_split = ct.split(str(Path('/')))
            cline3.append('/'.join(ct_split[1:]))
            cline2.append(ct)

    compounds_dict = {}
    compounds_di = {}
    # print(cline3)
    for line in aline2:
        stri = line.strip('\n')
        splitted = stri.split('/')
        for di in dir_list:
            # print(str(Path(f'{di}/{splitted[0]}/{splitted[1]}/mm.xyz')))
            if os.path.isfile(str(Path(f'{di}/{splitted[0]}/{splitted[1]}/mm.xyz'))):
                # print('here')
                molecule_dir = str(Path(f'{di}/{splitted[0]}'))
                if molecule_dir not in compounds_di:
                    compounds_dict[splitted[0]] = di.split(str(Path('/')))[-1]
                    compounds_di[molecule_dir] = {}
    # print(compounds_dict)
    # for compound in compounds_list:
    #     compounds_dict[compound] = {}
    for lin in aline2:
        stripp = lin.strip('\n')
        splitt = stripp.split('/')
        # print(str(Path('/')).join(splitt))
        # print(str('/'.join(splitt)))
        # print(type(splitt[-2]))
        try:
            site = int(splitt[-2])
            # print('here')
            if str('/'.join(splitt)) in cline3:
                # print('True')
                compounds_di[str(Path(f'{directory}/{compounds_dict[splitt[0]]}/{splitt[0]}'))][splitt[-2]] = True
            else:
                # print('False')
                compounds_di[str(Path(f'{directory}/{compounds_dict[splitt[0]]}/{splitt[0]}'))][splitt[-2]] = False
        except ValueError:
            continue
    print(compounds_di)
    for molecule in compounds_di.keys():
        filelist = os.listdir(molecule)
        previous_dict_path = Path(
            f'{str(Path("/")).join(molecule.split(str(Path("/")))[:-1])}/{molecule.split(str(Path("/")))[-1]}_results.json')
        with open(previous_dict_path, 'r') as prev:
            previous_dict = json.load(prev)
        for filename in filelist:
            if filename.endswith('.txt'):
                mol_name = filename.split('.')[0]
            elif filename.endswith('.smi'):
                try:
                    with open(str(Path(f'{molecule}/{filename}')), 'r') as smile:
                        mol_smi = smile.readlines()[0].strip('\n')
                except FileNotFoundError:
                    print(filename)
        mol_sites = list(compounds_di[molecule].keys())
        print(f'MOL NAME {mol_name}')
        print(f'MOL SMILES {mol_smi}')
        print(f'MOL SITES {mol_sites}')
        mole_prod_dict = {}
        for sit in mol_sites:
            if not os.path.isfile(str(Path(f'{molecule}/{sit}/prod.out'))):
                shutil.copy(Path(f'{molecule}/{sit}/am1.dat'), Path(f'{molecule}/{sit}/prod.dat'))
                print(str(Path(f'{molecule}/{sit}/am1.dat')))
                with open(str(Path(f'{molecule}/{sit}/prod.dat')), 'r') as prodin:
                    prod_lines = prodin.readlines()
                    prod_lines[0] = 'AM1 LET MMOK GEO-OK DISP CYCLES=10000 UHF\n'
                    carbon = prod_lines[-4].split()
                    carbon[1] = '1.600000\t'
                    carbon.append('\n')
                    t = " "
                    t = t.join(carbon)
                    prod_lines[-4] = t
                os.remove(Path(f'{molecule}/{sit}/prod.dat'))
                with open(str(Path(f'{molecule}/{sit}/prod.dat')), 'w+') as new:
                    for line in prod_lines:
                        new.write(str(line))

                check_call(['C:\\Program Files\\MOPAC\\MOPAC2016.exe', Path(f'{molecule}/{sit}/prod')], stdout=DEVNULL,
                           stderr=STDOUT)
                if read_mopac_output(Path(f'{molecule}/{sit}/prod.out'), Path(f'{molecule}/{sit}/prod.xyz'),
                                     Path(f'{molecule}/{sit}')) == 1:
                    check_call(['obabel', '-ixyz', Path(f'{molecule}/{sit}/prod.xyz'), '-osmi', '-O',
                                Path(f'{molecule}/{sit}/prod.smi')], stdout=DEVNULL, stderr=STDOUT)
                    with open(str(Path(f'{molecule}/{sit}/prod.smi')), 'r') as smi:
                        psmiles = smi.readlines()[0].split()[0]
                        mole_prod_dict[sit] = psmiles
            else:
                try:
                    with open(str(Path(f'{molecule}/{sit}/prod.smi')), 'r') as smi:
                        psmiles = smi.readlines()[0].split()[0]
                        mole_prod_dict[sit] = psmiles
                except FileNotFoundError:
                    mole_prod_dict[sit] = 'NO PRODUCT SMILES FOUND'

        # # print(f'MOL PROD DICT {mole_prod_dict}')
        # print('')
        # print('')
        # print(str(Path(f'{str(Path("/")).join(molecule.split(str(Path("/")))[:-1])}/{compounds_dict[molecule.split(str(Path("/")))[-1]]}')))
        # print(f'PREV MOL DICT {previous_dict}')
        molecule_dict = sd.store_data(compound_name=mol_name, smiles=mol_smi, compound_sites=mol_sites, radical='cf3',
                                      reactant_optimised=True, ts_converged=compounds_di[molecule],
                                      directory=molecule, hpc_calcs=True, prev_molecule_dict=previous_dict,
                                      product_smiles_dict=mole_prod_dict)
        results_dict[str(mol_name)] = molecule_dict
        print(results_dict)
    with open(str(Path(f'{directory}/test_results3.json')), 'w') as results:
        json.dump(results_dict, results)

with open(str(Path(f'{directory}/test_results3.json')), 'r') as results:
    results_dict = json.load(results)
    with open(str(Path(f'{directory}/test_results3.csv')), 'w') as csvfile:
        header = ['compound_name', 'smiles', 'radical', 'Site', 'Product_Smiles', 'Reaction_Smiles',
                  'hf_reactant_energy', 'hf_ts_energy', 'hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-',
                  'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding',
                  'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for com in results_dict.keys():
            for radical in results_dict[com]['radicals'].values():
                for key, value in radical.items():
                    if str.isdigit(key):
                        reaction_smiles = f'{results_dict[com]["smiles"]}>>{radical[key]["product_smiles"]}'
                        row = [results_dict[com]['compound_name'], results_dict[com]['smiles'],
                               radical['radical'], key, radical[key]['product_smiles'], reaction_smiles,
                               radical['hf_reactant_energy'], radical[key]['hf_ts_energy'],
                               radical[key]['hf_act_energy'], radical[key]['hf_ratio'], radical[key]['fa+'],
                               radical[key]['fa-'], radical[key]['fa0'], radical[key]['electron_density'],
                               radical[key]['electrostatic_potential'], radical[key]['diamagnetic_shielding'],
                               radical[key]['mulliken_charge'], radical[key]['bond_index'],
                               radical[key]['gross_population'], radical[key]['HF_C-H_Bond_Length'],
                               radical[key]['HF_TS_C-C_Bond_Length']]
                        writer.writerow(row)

# canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(mol_smi))
# mol = Chem.MolFromSmiles(canonical_smiles)
# # d = rdMolDraw2D.MolDraw2DSVG(250, 250)
# # print(mol_sites)
# # for atom in mol_sites:
# #     mol.GetAtomWithIdx(atom - 1).SetProp("_displayLabel", str(atom))
# # d.DrawMolecule(mol)
# # d.FinishDrawing()
# # molecule = d.GetDrawingText()
# rxn_smarts = '[cH1:1]>>[c:1][C:2]([F])([F])[F]'
# rxn = AllChem.ReactionFromSmarts(rxn_smarts)
# # product = rxn.RunReactants((com1, ))[0][0]
# glycerol_products = list(chain.from_iterable(rxn.RunReactants((mol, ))))
# # Chem.SanitizeMol(product)
# for product in glycerol_products:
#     print(Chem.MolToSmiles(product))
# # Draw.MolsToGridImage(product)
