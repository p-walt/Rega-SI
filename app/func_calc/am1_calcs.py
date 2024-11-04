import os
import csv
import hashlib
from rdkit.Chem import AllChem
from rdkit import Chem
import chemcoord as cc
import func_calc.add_radical as ar
import func_calc.store_data as sd
import func_calc.hf_calcs as hf
import json
from subprocess import DEVNULL, STDOUT, check_call
import shutil
import string, random
from app import pycharm
from pathlib import Path, PurePosixPath, PureWindowsPath

if pycharm == 1:
    start = os.getcwd()
    start_path = start + f'/ML-for-CH/app/Calculation_data/'
    mopac_path = "C:\\Program Files\\MOPAC\\MOPAC2016.exe"
    start_path = "ML-for-CH\\app\\Calculation_data"
else:
    mopac_path = '/app/mopac/MOPAC2016.exe'
    start_path = "/app/Calculation_data"


def calculated_previously(starting_dir, canonical_smiles, single=True, multiple=False, calculation_directory=None,
                          compound_name=None, molecule_sites_dict=None, compounds_dict=None, row=None, radical=None):
    # Previously Calculated key: 1 = Compound has been previously calculated in both AM1 and HF methods
    #                            2 = Compound has not been previosuly calculated in either AM1 or HF methods.
    #                            3 = Compound has been calculated in AM1 but not in HF method yet.
    previously_calculated = 2
    compound_yes_rad_no = False
    prev_calculated_dir = None
    if multiple and calculation_directory and molecule_sites_dict is not None and row is not None:
        for direc in os.listdir(starting_dir):
            if os.path.isdir(Path(f'{starting_dir}/{direc}')):
                if os.path.isfile(Path(f'{starting_dir}/{direc}/smiles.smi')):
                    with open(Path(f'{starting_dir}/{direc}/smiles.smi')) as sf:

                        check_smile = sf.readlines()[0]
                    if check_smile == canonical_smiles:
                        print('smiles match')
                        if os.path.isfile(Path(f'{starting_dir}/{direc}_results.json')):
                            with open(Path(f'{starting_dir}/{direc}_results.json')) as single_compound_json:
                                temp_molecule_dict = json.load(single_compound_json)
                                print(type(temp_molecule_dict))
                                if type(temp_molecule_dict) is not dict:
                                    temp_molecule_dict = list(json.load(single_compound_json).values())[0]
                            # print(temp_molecule_dict)
                            compounds_dict[str(row[0])] = temp_molecule_dict
                            mol = Chem.MolFromSmiles(canonical_smiles)
                            mol = Chem.AddHs(mol)
                            AllChem.EmbedMolecule(mol)
                            sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))

                            listsor = list(sor)
                            atom_no = []
                            for atom in listsor:
                                atom_no.append(list(atom))
                            flat_list = [item for sublist in atom_no for item in sublist]
                            # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                            aromatic_atom = []
                            for u in flat_list:
                                aromatic_atom.append(u + 1)
                            molecule_sites_dict[canonical_smiles] = aromatic_atom
                            # print(temp_molecule_dict)
                            try:
                                rads = list(temp_molecule_dict['radicals'].keys())
                            except KeyError:
                                # print(list(temp_molecule_dict.values())[0])
                                compounds_dict[str(row[0])] = list(temp_molecule_dict.values())[0]
                                rads = list(list(temp_molecule_dict.values())[0]['radicals'].keys())
                            for radicals in rads:
                                if radicals == row[2]:
                                    print('radical match')
                                    prev_calculated_dir = str(Path(f'{starting_dir}/{direc}'))
                                    for site in aromatic_atom:
                                        if os.path.isfile(Path(f'{starting_dir}/{direc}/{site}/hfoutp.out')):
                                            previously_calculated = 1
                                            return molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                    # Compound has been calculated in AM1 but not yet put through more expensive cluster calculation.
                                    previously_calculated = 3
                                    return molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                else:
                                    continue

                            compound_yes_rad_no = True
                            previously_calculated = 2
                            return molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                        else:
                            continue
                    else:
                        continue
                else:
                    for direc2 in os.listdir(Path(f'{starting_dir}/{direc}')):
                        if os.path.isdir(Path(f'{starting_dir}/{direc}/{direc2}')):
                            if os.path.isfile(Path(f'{starting_dir}/{direc}/{direc2}/smiles.smi')):
                                with open(Path(f'{starting_dir}/{direc}/{direc2}/smiles.smi')) as sf:
                                    check_smile = sf.readlines()[0]
                                if check_smile == canonical_smiles:
                                    print('smiles match')
                                    if os.path.isfile(Path(f'{starting_dir}/{direc}/{direc2}_results.json')):
                                        with open(
                                                Path(
                                                    f'{starting_dir}/{direc}/{direc2}_results.json')) as single_compound_json:
                                            temp_molecule_dict = json.load(single_compound_json)
                                            # print(temp_molecule_dict)
                                            # print(type(temp_molecule_dict))
                                            if type(temp_molecule_dict) is not dict:
                                                temp_molecule_dict = list(json.load(single_compound_json).values())[0]
                                        # print(temp_molecule_dict)
                                        # print(type(temp_molecule_dict))
                                        compounds_dict[str(row[0])] = temp_molecule_dict
                                        mol = Chem.MolFromSmiles(canonical_smiles)
                                        mol = Chem.AddHs(mol)
                                        AllChem.EmbedMolecule(mol)
                                        sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
                                        listsor = list(sor)
                                        atom_no = []
                                        for atom in listsor:
                                            atom_no.append(list(atom))
                                        flat_list = [item for sublist in atom_no for item in sublist]
                                        # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                                        aromatic_atom = []
                                        for u in flat_list:
                                            aromatic_atom.append(u + 1)
                                        molecule_sites_dict[canonical_smiles] = aromatic_atom
                                        try:
                                            rads = list(temp_molecule_dict['radicals'].keys())
                                        except KeyError:
                                            # print(list(temp_molecule_dict.values())[0])
                                            compounds_dict[str(row[0])] = list(temp_molecule_dict.values())[0]
                                            rads = list(list(temp_molecule_dict.values())[0]['radicals'].keys())
                                        for radicals in rads:
                                            if radicals == row[2]:
                                                print('radical match')
                                                prev_calculated_dir = str(Path(f'{starting_dir}/{direc}/{direc2}'))
                                                for site in aromatic_atom:
                                                    if os.path.isfile(
                                                            Path(f'{starting_dir}/{direc}/{direc2}/{site}/hfoutp.out')):
                                                        # Compound has been calculated in both AM1 and HF methods.
                                                        previously_calculated = 1
                                                        return molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                                # Compound has been calculated in AM1 but not yet put through more
                                                # expensive cluster calculation.
                                                previously_calculated = 3
                                                return molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                            else:
                                                continue

                                        compound_yes_rad_no = True
                                        previously_calculated = 2
                                        return molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                    else:
                                        continue
                                else:
                                    # print('moving on')
                                    continue
                            else:
                                continue
                        else:
                            continue

            else:
                continue
        if previously_calculated == 2:
            radical_calculated = False
            return molecule_sites_dict, compounds_dict, previously_calculated, radical_calculated, prev_calculated_dir
    elif single is True and compound_name is not None and molecule_sites_dict is None:
        aromatic_atom = []
        compounds_dict = {}
        for direc in os.listdir(starting_dir):

            if os.path.isdir(Path(f'{starting_dir}/{direc}')):
                if os.path.isfile(Path(f'{starting_dir}/{direc}/smiles.smi')):
                    with open(Path(f'{starting_dir}/{direc}/smiles.smi')) as sf:
                        check_smile = sf.readlines()[0]
                    if check_smile == canonical_smiles:
                        print('smiles match')

                        if os.path.isfile(Path(f'{starting_dir}/{direc}_results.json')):
                            with open(Path(f'{starting_dir}/{direc}_results.json')) as single_compound_json:
                                temp_molecule_dict = json.load(single_compound_json)
                                if type(temp_molecule_dict) is not dict:
                                    temp_molecule_dict = list(json.load(single_compound_json).values())[0]
                            # print(temp_molecule_dict)
                            compounds_dict[str(compound_name)] = temp_molecule_dict
                            mol = Chem.MolFromSmiles(canonical_smiles)
                            mol = Chem.AddHs(mol)
                            AllChem.EmbedMolecule(mol)
                            sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
                            listsor = list(sor)
                            atom_no = []
                            for atom in listsor:
                                atom_no.append(list(atom))
                            flat_list = [item for sublist in atom_no for item in sublist]
                            # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                            aromatic_atom = []
                            for u in flat_list:
                                aromatic_atom.append(u + 1)
                            # print(temp_molecule_dict.values())
                            try:
                                rads = list(temp_molecule_dict['radicals'].keys())
                            except KeyError:
                                # print(list(temp_molecule_dict.values())[0])
                                compounds_dict[str(compound_name)] = list(temp_molecule_dict.values())[0]
                                rads = list(list(temp_molecule_dict.values())[0]['radicals'].keys())
                            for radicals in rads:
                                if radicals == radical:
                                    print('radical match')
                                    prev_calculated_dir = str(Path(f'{starting_dir}/{direc}'))
                                    for site in aromatic_atom:
                                        if os.path.isfile(Path(f'{starting_dir}/{direc}/{site}/hfoutp.out')):
                                            previously_calculated = 1
                                            return aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                        # Compound has been calculated in AM1 but not yet put through more expensive cluster calculation.
                                    previously_calculated = 3
                                    return aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir
                                else:
                                    continue

                            compound_yes_rad_no = True
                            previously_calculated = 2
                            return aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir

                        else:
                            continue
                    else:
                        continue
                else:
                    for direc2 in os.listdir(Path(f'{starting_dir}/{direc}')):
                        if os.path.isdir(Path(f'{starting_dir}/{direc}/{direc2}')):
                            if os.path.isfile(Path(f'{starting_dir}/{direc}/{direc2}/smiles.smi')):
                                with open(Path(f'{starting_dir}/{direc}/{direc2}/smiles.smi')) as sf:
                                    check_smile = sf.readlines()[0]
                                if check_smile == canonical_smiles:
                                    print('smiles match')
                                    if os.path.isfile(Path(f'{starting_dir}/{direc}/{direc2}_results.json')):
                                        with open(
                                                Path(
                                                    f'{starting_dir}/{direc}/{direc2}_results.json')) as single_compound_json:
                                            temp_molecule_dict = json.load(single_compound_json)
                                            if type(temp_molecule_dict) is not dict:
                                                temp_molecule_dict = list(json.load(single_compound_json).values())[0]
                                        # print(temp_molecule_dict)
                                        compounds_dict[str(compound_name)] = temp_molecule_dict
                                        mol = Chem.MolFromSmiles(canonical_smiles)
                                        mol = Chem.AddHs(mol)
                                        AllChem.EmbedMolecule(mol)
                                        sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
                                        listsor = list(sor)
                                        atom_no = []
                                        for atom in listsor:
                                            atom_no.append(list(atom))
                                        flat_list = [item for sublist in atom_no for item in sublist]
                                        # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                                        aromatic_atom = []
                                        for u in flat_list:
                                            aromatic_atom.append(u + 1)
                                        try:
                                            rads = list(temp_molecule_dict['radicals'].keys())
                                        except KeyError:
                                            # print(list(temp_molecule_dict.values())[0])
                                            compounds_dict[str(compound_name)] = list(temp_molecule_dict.values())[0]
                                            rads = list(list(temp_molecule_dict.values())[0]['radicals'].keys())
                                        for radicals in rads:
                                            if radicals == radical:
                                                print('radical match')
                                                prev_calulated_dir = str(Path(f'{starting_dir}/{direc}/{direc2}'))
                                                for site in aromatic_atom:
                                                    if os.path.isfile(
                                                            Path(f'{starting_dir}/{direc}/{direc2}/{site}/hfoutp.out')):
                                                        previously_calculated = 1
                                                        return aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calulated_dir
                                                    # Compound has been calculated in AM1 but not yet put through more expensive cluster calculation.
                                                previously_calculated = 3
                                                return aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calulated_dir
                                            else:
                                                continue

                                        compound_yes_rad_no = True
                                        previously_calculated = 2
                                        return aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, prev_calculated_dir

                                    else:
                                        continue
                                else:
                                    # print('moving on')
                                    continue
                            else:
                                continue
                        else:
                            continue

            else:
                continue
        if previously_calculated == 2:
            radical_calculated = False
            return aromatic_atom, compounds_dict, previously_calculated, radical_calculated, prev_calculated_dir


def calculate_am1_energies_list(input_file, starting_directory, csv_dir, random_string, hpc_calcs):
    with open(f'{input_file}') as csv_file:
        compounds_dict = {}
        ts_converged = {}
        molecule_sites = {}
        compound_smiles_dict = {}
        prev_calc_dict = {}
        radicals_dict = {}
        csv_reader = list(csv.reader(csv_file, delimiter=','))
        temp_dict = {}
        new_mopac_2016 = True
        print(f'CSV READER IS {csv_reader}')
        for row in csv_reader:
            print(row)
            pos_charged = False
            reagent_optimised = True
            # Hash SMILES string for directory name
            canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(row[1]))
            compound_smiles_dict[canonical_smiles] = row[0]
            radicals_dict[canonical_smiles] = row[2]
            hash_object = hashlib.md5(canonical_smiles.encode())
            md5_hash = hash_object.hexdigest()
            molecule_sites_dict, compounds_dict, previously_calculated, compound_yes_rad_no, previously_calculated_compound_dir = calculated_previously(
                starting_dir=starting_directory, canonical_smiles=canonical_smiles, single=False, multiple=True,
                calculation_directory=csv_dir, compound_name=None, molecule_sites_dict=molecule_sites,
                compounds_dict=compounds_dict, row=row, radical=None)
            # print(f'COMPOUNDS_DICT {compounds_dict}')
            prev_calc_dict[str(row[0])] = previously_calculated
            # print(f'PREVIOUSLY CALCULATED {previously_calculated}')
            if previously_calculated == 2:

                with open(Path(f'{csv_dir}/test.dat'), 'w') as test:
                    test.write('AM1 LET MMOK GEO-OK PRECISE\n')
                    test.write('untitled.xyz\n\n')
                    test.write('C    0.000000  1    0.000000  1    0.000000  1     0   0   0\n')
                    test.write('H    1.070000  1    0.000000  1    0.000000  1     1   0   0\n')
                    test.write('H    1.070001  1  109.471386  1    0.000000  1     1   2   0\n')
                    test.write('H    1.070000  1  109.471408  1  239.999575  1     1   2   3\n')
                    test.write('H    1.070004  1  109.471340  1  119.999985  1     1   2   3\n')
                check_call([mopac_path, Path(f'{csv_dir}/test')], stdout=DEVNULL, stderr=STDOUT)
                if os.path.isfile(Path(f'{csv_dir}/test.arc')):
                    with open(Path(f'{csv_dir}/test.arc')) as tst:
                        for li in tst:
                            lin = li.strip('\n')
                            if 'TOTAL ENERGY' in lin:
                                new_mopac_2016 = False
                os.remove(Path(f'{csv_dir}/./test.arc'))
                os.remove(Path(f'{csv_dir}/./test.out'))
                os.remove(Path(f'{csv_dir}/./test.dat'))
                dispatcher_ts = {'cf3': ar.cf3, 'cf2h': ar.cf2h, 'ipr': ar.ipr}
                dispatcher_reactant = {'cf3': ar.cf3_reactant, 'cf2h': ar.cf2h_reactant, 'ipr': ar.ipr_reactant}
                if os.path.exists(Path(f'{csv_dir}/{md5_hash}')) is False:
                    print(f'MD5 HASH IS {md5_hash}')
                    os.mkdir(Path(f'{csv_dir}/{md5_hash}'))
                    with open(Path(f'{csv_dir}/{md5_hash}/smiles.smi'), 'w') as smiles_file:
                        smiles_file.write(str(canonical_smiles))
                    mol = Chem.MolFromSmiles(canonical_smiles)
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                    xyz = Chem.MolToXYZBlock(mol)
                    for item in str(canonical_smiles):
                        if item == '+':
                            pos_charged = True
                        else:
                            pass
                    try:
                        with open(Path(f'{csv_dir}/{md5_hash}/{row[0]}.txt'), 'w+') as mole:
                            mole.write(row[0])
                    except FileNotFoundError:
                        pass
                    # Site of reaction (SoR) is an aromatic carbon with a hydrogen bonded. Any already functionalised carbons in
                    # a ring are ignored.
                    sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
                    listsor = list(sor)
                    atom_no = []
                    for atom in listsor:
                        atom_no.append(list(atom))
                    flat_list = [item for sublist in atom_no for item in sublist]
                    # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                    aromatic_atom = []
                    for u in flat_list:
                        aromatic_atom.append(u + 1)
                    with open(Path(f'{csv_dir}/{md5_hash}/tmp.xyz'), 'w+') as v:
                        v.write(str(xyz))
                    if not aromatic_atom:
                        reagent_optimised = False
                        temp_dict[str(row[0])] = [row[0], canonical_smiles, aromatic_atom, row[2], reagent_optimised,
                                                  ts_converged, str(Path(f'{csv_dir}/{md5_hash}/')),
                                                  compound_yes_rad_no, str(Path(f'{csv_dir}/{md5_hash}'))]
                        continue
                    elif os.path.getsize(Path(f'{csv_dir}/{md5_hash}/tmp.xyz')) == 0:
                        reagent_optimised = False
                        temp_dict[str(row[0])] = [row[0], canonical_smiles, aromatic_atom, row[2], reagent_optimised,
                                                  ts_converged, str(Path(f'{csv_dir}/{md5_hash}/')),
                                                  compound_yes_rad_no, str(Path(f'{csv_dir}/{md5_hash}'))]
                    # compound = Species(name=row[0], charge=0, mult=1, atoms=xyz_file_to_atoms('tmp.xyz'))
                    # xyz_file_to_atoms('tmp.xyz')
                    # os.system(
                    #     'obabel tmp.xyz -O mm.xyz --minimize --sd --crit 1e-7 --steps 9999 --ff MMFF94 > /dev/null 2>&1')
                    check_call(['obabel', Path(f'{csv_dir}/{md5_hash}/tmp.xyz'), '-O',
                                Path(f'{csv_dir}/{md5_hash}/mm.xyz'), '--minimize', '--sd', '--crit',
                                '1e-7', '--steps', '9999', '--ff', 'MMFF94'], stdout=DEVNULL, stderr=STDOUT)
                    if os.path.getsize(Path(f'{csv_dir}/{md5_hash}/mm.xyz')) != 0:
                        check_call(
                            ['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/mm.xyz'), '-omopin', '-O',
                             Path(f'{csv_dir}/{md5_hash}/am1.dat')], stdout=DEVNULL,
                            stderr=STDOUT)
                    else:
                        check_call(
                            ['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/tmp.xyz'), '-omopin', '-O',
                             Path(f'{csv_dir}/{md5_hash}/am1.dat')], stdout=DEVNULL,
                            stderr=STDOUT)
                    os.mkdir(Path(f'{csv_dir}/{md5_hash}/reagent'))
                    os.mkdir(Path(f'{csv_dir}/{md5_hash}/fukui'))
                    if os.path.getsize(Path(f'{csv_dir}/{md5_hash}/mm.xyz')) != 0:
                        shutil.copy(Path(f'{csv_dir}/{md5_hash}/mm.xyz'), Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'))
                    else:
                        shutil.copy(Path(f'{csv_dir}/{md5_hash}/tmp.xyz'), Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'))
                    shutil.copy(Path(f'{csv_dir}/{md5_hash}/am1.dat'),
                                Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'))

                    dispatcher_reactant[row[2]](Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'),
                                                Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'),
                                                aromatic_atom[0], Path(f'{csv_dir}/{md5_hash}/reagent/'),
                                                pos_charged, reagent_optimised, new_mopac_2016, multiple=True)
                    if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.arc')):
                        print('Reagent optimised')
                        reagent_optimised = True
                    else:
                        if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.den')):
                            os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.den'))
                        else:
                            pass
                        if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.out')):
                            os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.out'))
                        else:
                            pass
                        if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.res')):
                            os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.res'))
                        else:
                            pass
                        if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat')):
                            os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'))

                        check_call(['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'), '-omopin',
                                    '-O', Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat')], stdout=DEVNULL, stderr=STDOUT)
                        # os.system('obabel -ixyz mm.xyz -omopin -O am1.dat > /dev/null 2>&1')
                        if len(aromatic_atom) != 1:

                            dispatcher_reactant[row[2]](Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'),
                                                        Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'),
                                                        aromatic_atom[1], Path(f'{csv_dir}/{md5_hash}/reagent/'),
                                                        pos_charged, reagent_optimised, new_mopac_2016, multiple=True)
                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.arc')):
                                print('Reagent optimised')
                                reagent_optimised = True
                            else:
                                print('Reagent not optimised')
                                reagent_optimised = False
                                if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.den')):
                                    os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.den'))
                                else:
                                    pass
                                if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.out')):
                                    os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.out'))
                                else:
                                    pass
                                if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.res')):
                                    os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.res'))
                                else:
                                    pass
                                if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat')):
                                    os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'))
                                check_call(['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'),
                                            '-omopin', '-O', Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat')],
                                           stdout=DEVNULL, stderr=STDOUT)
                                # os.system('obabel -ixyz mm.xyz -omopin -O am1.dat > /dev/null 2>&1')
                                dispatcher_reactant[row[2]](Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'),
                                                            Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'),
                                                            aromatic_atom[0], Path(f'{csv_dir}/{md5_hash}/reagent/'),
                                                            pos_charged, reagent_optimised, new_mopac_2016,
                                                            multiple=True)
                        else:
                            print('Reagent not optimised')
                            reagent_optimised = False
                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.den')):
                                os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.den'))
                            else:
                                pass
                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.out')):
                                os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.out'))
                            else:
                                pass
                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.res')):
                                os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.res'))
                            else:
                                pass
                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat')):
                                os.remove(Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'))
                            check_call(['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'),
                                        '-omopin', '-O', Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat')], stdout=DEVNULL,
                                       stderr=STDOUT)
                            # os.system('obabel -ixyz mm.xyz -omopin -O am1.dat > /dev/null 2>&1')
                            dispatcher_reactant[row[2]](Path(f'{csv_dir}/{md5_hash}/reagent/am1.dat'),
                                                        Path(f'{csv_dir}/{md5_hash}/reagent/mm.xyz'),
                                                        aromatic_atom[0], Path(f'{csv_dir}/{md5_hash}/reagent/'),
                                                        pos_charged, reagent_optimised, new_mopac_2016, multiple=True)
                    # compound = cc.Cartesian.read_xyz('mm.xyz', start_index=1)
                    # connection_table = compound.get_bonds()
                    if os.path.getsize(Path(f'{csv_dir}/{md5_hash}/mm.xyz')) != 0:
                        check_call(
                            ['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/mm.xyz'), '-omopin', '-O',
                             Path(f'{csv_dir}/{md5_hash}/am1.dat')], stdout=DEVNULL,
                            stderr=STDOUT)
                    else:
                        check_call(
                            ['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/tmp.xyz'), '-omopin', '-O',
                             Path(f'{csv_dir}/{md5_hash}/am1.dat')], stdout=DEVNULL,
                            stderr=STDOUT)
                    # os.system('obabel -ixyz mm.xyz -omopin -O am1.dat > /dev/null 2>&1')
                    for site in aromatic_atom:
                        if os.path.isdir(Path(f'{csv_dir}/{md5_hash}/{site}')):
                            pass
                        else:

                            os.mkdir(Path(f'{csv_dir}/{md5_hash}/{site}'))
                            if os.path.getsize(Path(f'{csv_dir}/{md5_hash}/mm.xyz')) != 0:
                                shutil.copy(Path(f'{csv_dir}/{md5_hash}/mm.xyz'),
                                            Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'))
                            else:
                                shutil.copy(Path(f'{csv_dir}/{md5_hash}/tmp.xyz'),
                                            Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'))
                            shutil.copy(Path(f'{csv_dir}/{md5_hash}/am1.dat'),
                                        Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat'))
                            # os.system('cp ../am1.dat ../mm.xyz .')
                            dispatcher_ts[row[2]](Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat'),
                                                  Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'), site,
                                                  Path(f'{csv_dir}/{md5_hash}/{site}/'), pos_charged,
                                                  False, True, new_mopac_2016)

                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/am1.arc')) is True:
                                ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/am1.out'),
                                                     Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                     Path(f'{csv_dir}/{md5_hash}/{site}'))
                                if ar.mopac_distance_checks(Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'), f'{row[2]}',
                                                            site) is False:

                                    os.remove(Path(f'{csv_dir}/{md5_hash}/{site}/am1.out'))
                                    os.remove(Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat'))
                                    check_call(['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'),
                                                '-omopin', '-O',
                                                Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat')],
                                               stdout=DEVNULL, stderr=STDOUT)
                                    dispatcher_ts[row[2]](Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat'),
                                                          Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'),
                                                          site, Path(f'{csv_dir}/{md5_hash}/{site}/'),
                                                          pos_charged, True, False, new_mopac_2016)
                                    if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/am1.arc')) is \
                                            True:

                                        ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/am1.out'),
                                                             Path(f'{csv_dir}/{md5_hash}/{site}/constrained.xyz'),
                                                             Path(f'{csv_dir}/{md5_hash}/{site}'))
                                        if ar.mopac_distance_checks(
                                                Path(f'{csv_dir}/{md5_hash}/{site}/constrained.xyz'), f'{row[2]}',
                                                site) is True:
                                            check_call(['obabel', '-ixyz',
                                                        Path(f'{csv_dir}/{md5_hash}/{site}/constrained.xyz'),
                                                        '-omopin', '-O', Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat')],
                                                       stdout=DEVNULL, stderr=STDOUT)
                                            with open(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'), 'r+') as relaxed:
                                                lines = relaxed.readlines()
                                                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF CYCLES=10000\n'

                                            os.remove(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'))
                                            with open(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'), 'w+') as file:
                                                for line in lines:
                                                    file.write(str(line))
                                            check_call([mopac_path, Path(f'{csv_dir}/{md5_hash}/{site}/ts2')],
                                                       stdout=DEVNULL, stderr=STDOUT)
                                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.arc')) is True:

                                                ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.out'),
                                                                     Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                                     Path(f'{csv_dir}/{md5_hash}/{site}'))
                                                ts_converged[site] = ar.mopac_freq_check(
                                                    Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'), pos_charged, row[2],
                                                    site, Path(f'{csv_dir}/{md5_hash}/{site}'), new_mopac_2016,
                                                    multiple=True)

                                            else:
                                                ar.tweak_distance(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'),
                                                                  f'{row[2]}', Path(f'{csv_dir}/{md5_hash}/{site}/'))
                                                if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.arc')):
                                                    ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.out'),
                                                                         Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                                         Path(f'{csv_dir}/{md5_hash}/{site}'))
                                                    ts_converged[site] = ar.mopac_freq_check(
                                                        Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                        pos_charged, row[2], site, Path(f'{csv_dir}/{md5_hash}/{site}'),
                                                        new_mopac_2016, multiple=True)
                                                else:
                                                    ts_converged[site] = False
                                                    print('Did not converge')

                                        else:
                                            ts_converged[site] = False
                                            print('Did not converge')
                                    else:
                                        ts_converged[site] = False
                                        print('Did not converge')
                                else:
                                    ts_converged[site] = ar.mopac_freq_check(
                                        Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'), pos_charged, row[2], site,
                                        Path(f'{csv_dir}/{md5_hash}/{site}'), new_mopac_2016, multiple=True)

                            else:
                                os.remove(Path(f'{csv_dir}/{md5_hash}/{site}/am1.out'))
                                os.remove(Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat'))
                                check_call(['obabel', '-ixyz', Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'),
                                            '-omopin', '-O', Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat')],
                                           stdout=DEVNULL, stderr=STDOUT)
                                dispatcher_ts[row[2]](Path(f'{csv_dir}/{md5_hash}/{site}/am1.dat'),
                                                      Path(f'{csv_dir}/{md5_hash}/{site}/mm.xyz'), site,
                                                      Path(f'{csv_dir}/{md5_hash}/{site}/'),
                                                      pos_charged, True, False, new_mopac_2016)
                                if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/am1.arc')) is True:

                                    ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/am1.out'),
                                                         Path(f'{csv_dir}/{md5_hash}/{site}/constrained.xyz'),
                                                         Path(f'{csv_dir}/{md5_hash}/{site}'))
                                    if ar.mopac_distance_checks(Path(f'{csv_dir}/{md5_hash}/{site}/constrained.xyz'),
                                                                f'{row[2]}', site) is True:
                                        check_call(['obabel', '-ixyz',
                                                    Path(f'{csv_dir}/{md5_hash}/{site}/constrained.xyz'), '-omopin',
                                                    '-O', Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat')], stdout=DEVNULL,
                                                   stderr=STDOUT)
                                        with open(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'), 'r+') as relaxed:
                                            lines = relaxed.readlines()
                                            if new_mopac_2016 is False:
                                                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF CYCLES=10000\n'
                                            else:
                                                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF CYCLES=10000\n'

                                        os.remove(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'))
                                        with open(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'), 'w+') as file:
                                            for line in lines:
                                                file.write(str(line))
                                        check_call([mopac_path, Path(f'{csv_dir}/{md5_hash}/{site}/ts2')],
                                                   stdout=DEVNULL,
                                                   stderr=STDOUT)
                                        if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.arc')) is True:

                                            ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.out'),
                                                                 Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                                 Path(f'{csv_dir}/{md5_hash}/{site}'))
                                            ts_converged[site] = ar.mopac_freq_check(
                                                Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                pos_charged, row[2], site, Path(f'{csv_dir}/{md5_hash}/{site}'),
                                                new_mopac_2016, multiple=True)

                                        else:
                                            ar.tweak_distance(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.dat'), f'{row[2]}',
                                                              Path(f'{csv_dir}/{md5_hash}/{site}/'))
                                            if os.path.isfile(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.arc')):
                                                ar.read_mopac_output(Path(f'{csv_dir}/{md5_hash}/{site}/ts2.out'),
                                                                     Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'),
                                                                     Path(f'{csv_dir}/{md5_hash}/{site}'))
                                                ts_converged[site] = ar.mopac_freq_check(
                                                    Path(f'{csv_dir}/{md5_hash}/{site}/ts.xyz'), pos_charged, row[2],
                                                    site, Path(f'{csv_dir}/{md5_hash}/{site}'),
                                                    new_mopac_2016, multiple=True)
                                            else:
                                                ts_converged[site] = False
                                                print('Did not converge')

                                    else:
                                        ts_converged[site] = False
                                        print('Did not converge')

                                else:
                                    ts_converged[site] = False
                                    print('Did not converge')
                    # data for this compound and its sites are temporarily stored in temp_dict for use later in the
                    # store data function once the HPC calculations have been finished
                    temp_dict[str(row[0])] = [row[0], canonical_smiles, aromatic_atom, row[2], reagent_optimised,
                                              ts_converged, str(Path(f'{csv_dir}/{md5_hash}/')),
                                              compound_yes_rad_no, str(Path(f'{csv_dir}/{md5_hash}'))]
                    print(temp_dict)
            elif previously_calculated == 3:
                new_molecule_directory = previously_calculated_compound_dir
                print('ADDING HPC CALCS MULTIPLE')
                for item in str(canonical_smiles):
                    if item == '+':
                        pos_charged = True
                    else:
                        pass
                mol = Chem.MolFromSmiles(canonical_smiles)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
                listsor = list(sor)
                atom_no = []
                for atom in listsor:
                    atom_no.append(list(atom))
                flat_list = [item for sublist in atom_no for item in sublist]
                # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                aromatic_atom = []
                for u in flat_list:
                    aromatic_atom.append(u + 1)
                print(f'AROMATIC ATOM IS {aromatic_atom}')
                if os.path.isfile(str(Path(f'{new_molecule_directory}/reagent/am1.arc'))):
                    reagent_optimised = True
                else:
                    reagent_optimised = False
                for site in aromatic_atom:
                    if os.path.isfile(str(Path(f'{new_molecule_directory}/{site}/am1_converged.txt'))):
                        ts_converged[site] = True
                    else:
                        ts_converged[site] = False
                if hpc_calcs:
                    hf_jobs = str(Path(f'{str(input_file).split(".")[0]}/hf_jobs.txt'))
                    print(f'HF JOBS ARE {hf_jobs}')
                    if os.path.isfile(hf_jobs) is True:
                        with open(hf_jobs, 'r') as fil:
                            lins = fil.readlines()
                        for li in lins:
                            li = str(Path(li)).split(str(Path('/')))

                            # removes the previous hf_jobs.txt file if other compounds calculated in the directory before.
                            if new_molecule_directory.split(str(Path("/")))[-1] == li[-3]:
                                print('REMOVING FILE')
                                os.remove(hf_jobs)
                                break
                        with open(hf_jobs, 'a') as fi:
                            for site in aromatic_atom:
                                print(f'NEW MOLECULE DIRECTORY {new_molecule_directory}')
                                print(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                                if os.path.isfile(Path(f'{new_molecule_directory}/{site}/hf.sdf')):
                                    fi.write(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                            if os.path.isfile(Path(f'{new_molecule_directory}/reagent/hf.sdf')):
                                fi.write(str(Path(f'{new_molecule_directory}/reagent/hf.sdf\n')))
                            fi.write(str(Path(f'{new_molecule_directory}/fukui/hf.sdf\n')))
                    else:
                        with open(hf_jobs, 'w') as fil:
                            for site in aromatic_atom:
                                print(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                                if os.path.isfile(Path(f'{new_molecule_directory}/{site}/hf.sdf')):
                                    fil.write(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                            if os.path.isfile(Path(f'{new_molecule_directory}/reagent/hf.sdf')):
                                fil.write(str(Path(f'{new_molecule_directory}/reagent/hf.sdf\n')))
                            fil.write(str(Path(f'{new_molecule_directory}/fukui/hf.sdf\n')))
                # hf.run_pyreact(f'-d a --nwchem --mult 2 --FC --NoEigen -F UHF --nProc 8 --TS --Rega --array {hf_jobs}')
                temp_dict[str(row[0])] = [row[0], canonical_smiles, aromatic_atom, row[2], reagent_optimised,
                                          ts_converged, previously_calculated_compound_dir,
                                          compound_yes_rad_no, previously_calculated_compound_dir]
                print(temp_dict)
            else:
                new_molecule_directory = previously_calculated_compound_dir
                for item in str(canonical_smiles):
                    if item == '+':
                        pos_charged = True
                    else:
                        pass
                mol = Chem.MolFromSmiles(canonical_smiles)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                sor = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
                listsor = list(sor)
                atom_no = []
                for atom in listsor:
                    atom_no.append(list(atom))
                flat_list = [item for sublist in atom_no for item in sublist]
                # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
                aromatic_atom = []
                for u in flat_list:
                    aromatic_atom.append(u + 1)
                if os.path.isfile(str(Path(f'{new_molecule_directory}/reagent/am1.arc'))):
                    reagent_optimised = True
                else:
                    reagent_optimised = False
                for site in aromatic_atom:
                    if os.path.isfile(str(Path(f'{new_molecule_directory}/{site}/am1_converged.txt'))):
                        ts_converged[site] = True
                    else:
                        ts_converged[site] = False
                molecule_sites[canonical_smiles] = aromatic_atom
                temp_dict[str(row[0])] = [row[0], canonical_smiles, aromatic_atom, row[2], reagent_optimised,
                                          ts_converged, previously_calculated_compound_dir,
                                          compound_yes_rad_no, previously_calculated_compound_dir]

        print(f'CSV DIR IS {csv_dir}')
        hf_jobs = str(Path(f'{str(input_file).split(".")[0]}/hf_jobs.txt'))
        print(hf_jobs)
        if hpc_calcs:
            hf.run_pyreact(f'-d a --nwchem --mult 2 --FC --NoEigen -F UHF --nProc 8 --TS --Rega --array {hf_jobs}')
            print('finished pyReact')
        # print(csv_reader)
        for row in csv_reader:
            # print('here')
            # print(temp_dict)
            if str(row[0]) in compounds_dict:
                prev_calcs = compounds_dict[str(row[0])]
            else:
                prev_calcs = None
            if str(row[0]) in temp_dict:
                if temp_dict[str(row[0])][7] is True:
                    molecule_dict = sd.store_data(temp_dict[str(row[0])][0], temp_dict[str(row[0])][1],
                                                  temp_dict[str(row[0])][2], temp_dict[str(row[0])][3],
                                                  temp_dict[str(row[0])][4], temp_dict[str(row[0])][5],
                                                  temp_dict[str(row[0])][6], hpc_calcs, prev_calcs)
                    # print(molecule_dict)
                else:
                    # print(temp_dict[str(row[0])][0], temp_dict[str(row[0])][1],
                    #       temp_dict[str(row[0])][2], temp_dict[str(row[0])][3],
                    #       temp_dict[str(row[0])][4], temp_dict[str(row[0])][5],
                    #       temp_dict[str(row[0])][6], hpc_calcs, prev_calcs)

                    molecule_dict = sd.store_data(temp_dict[str(row[0])][0], temp_dict[str(row[0])][1],
                                                  temp_dict[str(row[0])][2], temp_dict[str(row[0])][3],
                                                  temp_dict[str(row[0])][4], temp_dict[str(row[0])][5],
                                                  temp_dict[str(row[0])][6], hpc_calcs, prev_calcs)
                    # print(molecule_dict)
                compounds_dict[temp_dict[str(row[0])][0]] = molecule_dict
                # print(compounds_dict)
                with open(Path(f'{temp_dict[str(row[0])][8]}_results.json'), 'w') as single_compound_json:
                    json.dump(molecule_dict, single_compound_json)
                molecule_sites[temp_dict[str(row[0])][1]] = temp_dict[str(row[0])][2]
            else:
                continue
        with open(Path(f'{csv_dir}_results.json'), 'w') as json_file:
            json.dump(compounds_dict, json_file)
        with open(Path(f'{csv_dir}_results.csv'), 'w') as csvfile:
            if hpc_calcs:
                header = ['compound_name', 'smiles', 'radical', 'Site', 'hf_reactant_energy', 'hf_ts_energy',
                          'hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0']
            else:
                header = ['compound_name', 'smiles', 'radical', 'Site', 'am1_reactant_energy', 'am1_ts_energy',
                          'am1_activation_energy', 'am1_ratio']
            writer = csv.writer(csvfile)
            writer.writerow(header)
            for com in compounds_dict.keys():
                for radical in compounds_dict[com]['radicals'].values():
                    for key, value in radical.items():
                        if str.isdigit(key):
                            if hpc_calcs:
                                row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                       radical['radical'], key,
                                       radical['hf_reactant_energy'], radical[key]['hf_ts_energy'],
                                       radical[key]['hf_act_energy'], radical[key]['hf_ratio'],
                                       radical[key]['fa+'],
                                       radical[key]['fa-'], radical[key]['fa0']]
                            else:
                                row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                       radical['radical'], key,
                                       radical['am1_reactant_energy'], radical[key]['am1_ts_energy'],
                                       radical[key]['am1_act_energy'], radical[key]['am1_ratio']]
                            writer.writerow(row)
    return compounds_dict, molecule_sites, compound_smiles_dict, radicals_dict


def calculate_am1_energies_single(compound_name, smiles, radical, calc_directory, new_molecule_directory, hpc_calcs):
    compounds_dict = {}
    ts_converged = {}
    pos_charged = False
    reagent_optimised = True
    new_mopac_2016 = True
    # Hash SMILES string for directory name
    canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    aromatic_atom, compounds_dict, previously_calculated, compound_yes_rad_no, previously_calculated_compound_dir = calculated_previously(
        starting_dir=calc_directory,
        canonical_smiles=canonical_smiles,
        single=True, multiple=False,
        calculation_directory=None,
        compound_name=compound_name,
        molecule_sites_dict=None,
        compounds_dict=compounds_dict,
        row=None, radical=radical)
    if previously_calculated == 2:
        dispatcher_ts = {'cf3': ar.cf3, 'cf2h': ar.cf2h, 'ipr': ar.ipr}
        dispatcher_reactant = {'cf3': ar.cf3_reactant, 'cf2h': ar.cf2h_reactant, 'ipr': ar.ipr_reactant}
        with open(Path(f'{calc_directory}/test.dat'), 'w') as test:
            test.write('AM1 LET MMOK GEO-OK PRECISE\n')
            test.write('untitled.xyz\n\n')
            test.write('C    0.000000  1    0.000000  1    0.000000  1     0   0   0\n')
            test.write('H    1.070000  1    0.000000  1    0.000000  1     1   0   0\n')
            test.write('H    1.070001  1  109.471386  1    0.000000  1     1   2   0\n')
            test.write('H    1.070000  1  109.471408  1  239.999575  1     1   2   3\n')
            test.write('H    1.070004  1  109.471340  1  119.999985  1     1   2   3\n')

        check_call([mopac_path, Path(f'{calc_directory}/test')], stdout=DEVNULL, stderr=STDOUT)
        if os.path.isfile(Path(f'{calc_directory}/test.arc')):
            with open(Path(f'{calc_directory}/test.arc')) as tst:
                for li in tst:
                    lin = li.strip('\n')
                    if 'TOTAL ENERGY' in lin:
                        new_mopac_2016 = False
        os.remove(Path(f'{calc_directory}/test.arc'))
        os.remove(Path(f'{calc_directory}/test.out'))
        os.remove(Path(f'{calc_directory}/test.dat'))
        if os.path.exists(new_molecule_directory) is False:

            os.mkdir(str(new_molecule_directory))
            with open(Path(f'{new_molecule_directory}/smiles.smi'), 'w') as smiles_file:
                smiles_file.write(str(canonical_smiles))
            mol = Chem.MolFromSmiles(canonical_smiles)
            mol = Chem.AddHs(mol)

            AllChem.EmbedMolecule(mol)
            xyz = Chem.MolToXYZBlock(mol)
            for item in smiles:
                if item == '+':
                    pos_charged = True
                else:
                    pass
            with open(Path(f'{new_molecule_directory}/{compound_name}.txt'), 'w+') as mole:
                mole.write(compound_name)
            # Site of reaction (SoR) is an aromatic carbon with a hydrogen bonded. Any already functionalised carbons in
            # a ring are ignored.
            SoR = mol.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
            listSoR = list(SoR)
            atom_no = []
            for atom in listSoR:
                atom_no.append(list(atom))
            flat_list = [item for sublist in atom_no for item in sublist]
            # aromatic_atom is the actual site of reaction atom number rather than atom index which starts at 0
            aromatic_atom = []
            for u in flat_list:
                aromatic_atom.append(u + 1)
            with open(Path(f'{new_molecule_directory}/tmp.xyz'), 'w+') as v:
                v.write(str(xyz))
            if not aromatic_atom:
                reagent_optimised = False

                molecule_dict = sd.store_data(compound_name, canonical_smiles, aromatic_atom, radical,
                                              reagent_optimised, ts_converged, Path(f'{new_molecule_directory}/'))
                compounds_dict[str(compound_name)] = molecule_dict
                with open(f'{new_molecule_directory}_results.json', 'w') as json_file:
                    json.dump(compounds_dict, json_file)
                return compounds_dict, aromatic_atom
            print(f'AROMATIC ATOM {aromatic_atom}')
            # compound = Species(name=row[0], charge=0, mult=1, atoms=xyz_file_to_atoms('tmp.xyz'))
            # xyz_file_to_atoms('tmp.xyz')
            check_call(
                ['obabel', Path(f'{new_molecule_directory}/tmp.xyz'), '-O', Path(f'{new_molecule_directory}/mm.xyz'),
                 '--minimize', '--sd', '--crit', '1e-7', '--steps', '9999', '--ff',
                 'MMFF94'], stdout=DEVNULL, stderr=STDOUT)
            if os.path.getsize(Path(f'{new_molecule_directory}/mm.xyz')) != 0:
                check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/mm.xyz'), '-omopin', '-O',
                            Path(f'{new_molecule_directory}/am1.dat')], stdout=DEVNULL, stderr=STDOUT)
            else:
                check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/tmp.xyz'), '-omopin', '-O',
                            Path(f'{new_molecule_directory}/am1.dat')], stdout=DEVNULL, stderr=STDOUT)
            os.mkdir(Path(f'{new_molecule_directory}/reagent'))
            os.mkdir(Path(f'{new_molecule_directory}/fukui'))
            if os.path.getsize(Path(f'{new_molecule_directory}/mm.xyz')) != 0:
                shutil.copy(Path(f'{new_molecule_directory}/mm.xyz'), Path(f'{new_molecule_directory}/reagent/mm.xyz'))
            else:
                shutil.copy(Path(f'{new_molecule_directory}/tmp.xyz'), Path(f'{new_molecule_directory}/reagent/mm.xyz'))
            shutil.copy(Path(f'{new_molecule_directory}/am1.dat'), Path(f'{new_molecule_directory}/reagent/am1.dat'))
            # os.system('cp ../mm.xyz .')
            # check_call(['obabel', 'tmp.xyz', '-O', 'mm.xyz',
            #             '--minimize', '--sd', '--crit', '1e-7', '--steps', '9999', '--ff', 'MMFF94'], stdout=DEVNULL
            #            , stderr=STDOUT)
            dispatcher_reactant[radical](Path(f'{new_molecule_directory}/reagent/am1.dat'),
                                         Path(f'{new_molecule_directory}/reagent/mm.xyz'), aromatic_atom[0],
                                         Path(f'{new_molecule_directory}/reagent/'), pos_charged, reagent_optimised,
                                         new_mopac_2016, multiple=False)
            if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.arc')):
                print('Reagent optimised')
                reagent_optimised = True
            else:
                if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.den')):
                    os.remove(Path(f'{new_molecule_directory}/reagent/am1.den'))
                else:
                    pass
                if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.out')):
                    os.remove(Path(f'{new_molecule_directory}/reagent/am1.out'))
                else:
                    pass
                if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.res')):
                    os.remove(Path(f'{new_molecule_directory}/reagent/am1.res'))
                else:
                    pass
                if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.dat')):
                    os.remove(Path(f'{new_molecule_directory}/reagent/am1.dat'))
                    check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/reagent/mm.xyz'), '-omopin', '-O',
                                Path(f'{new_molecule_directory}/reagent/am1.dat')], stdout=DEVNULL, stderr=STDOUT)
                if len(aromatic_atom) != 1:

                    dispatcher_reactant[radical](Path(f'{new_molecule_directory}/reagent/am1.dat'),
                                                 Path(f'{new_molecule_directory}/reagent/mm.xyz'),
                                                 aromatic_atom[1], Path(f'{new_molecule_directory}/reagent/'),
                                                 pos_charged, reagent_optimised, new_mopac_2016, multiple=False)
                    if os.path.isfile('am1.arc'):
                        print('Reagent optimised')
                        reagent_optimised = True
                    else:
                        print('Reagent not optimised')
                        reagent_optimised = False
                        if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.den')):
                            os.remove(Path(f'{new_molecule_directory}/reagent/am1.den'))
                        else:
                            pass
                        if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.out')):
                            os.remove(Path(f'{new_molecule_directory}/reagent/am1.out'))
                        else:
                            pass
                        if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.res')):
                            os.remove(Path(f'{new_molecule_directory}/reagent/am1.res'))
                        else:
                            pass
                        if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.dat')):
                            os.remove(Path(f'{new_molecule_directory}/reagent/am1.dat'))
                        check_call(
                            ['obabel', '-ixyz', Path(f'{new_molecule_directory}/reagent/mm.xyz'), '-omopin', '-O',
                             Path(f'{new_molecule_directory}/reagent/am1.dat')], stdout=DEVNULL, stderr=STDOUT)
                        # os.system('obabel -ixyz mm.xyz -omopin -O am1.dat > /dev/null 2>&1')
                        dispatcher_reactant[radical](Path(f'{new_molecule_directory}/reagent/am1.dat'),
                                                     Path(f'{new_molecule_directory}/reagent/mm.xyz'),
                                                     aromatic_atom[0], Path(f'{new_molecule_directory}/reagent/'),
                                                     pos_charged, reagent_optimised, new_mopac_2016, multiple=False)
                else:
                    print('Reagent not optimised')
                    reagent_optimised = False
                    if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.den')):
                        os.remove(Path(f'{new_molecule_directory}/reagent/am1.den'))
                    else:
                        pass
                    if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.out')):
                        os.remove(Path(f'{new_molecule_directory}/reagent/am1.out'))
                    else:
                        pass
                    if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.res')):
                        os.remove(Path(f'{new_molecule_directory}/reagent/am1.res'))
                    else:
                        pass
                    if os.path.isfile(Path(f'{new_molecule_directory}/reagent/am1.dat')):
                        os.remove(Path(f'{new_molecule_directory}/reagent/am1.dat'))
                    check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/reagent/mm.xyz'), '-omopin', '-O',
                                Path(f'{new_molecule_directory}/reagent/am1.dat')], stdout=DEVNULL, stderr=STDOUT)
                    # os.system('obabel -ixyz mm.xyz -omopin -O am1.dat > /dev/null 2>&1')
                    dispatcher_reactant[radical]('am1.dat', 'mm.xyz', aromatic_atom[0],
                                                 Path(f'{new_molecule_directory}/reagent/'),
                                                 pos_charged, reagent_optimised, new_mopac_2016, multiple=False)
            for site in aromatic_atom:
                if os.path.isdir(site):
                    pass
                else:

                    os.mkdir(Path(f'{new_molecule_directory}/{site}'))
                    if os.path.getsize(Path(f'{new_molecule_directory}/mm.xyz')) != 0:
                        shutil.copy(Path(f'{new_molecule_directory}/mm.xyz'),
                                    Path(f'{new_molecule_directory}/{site}/mm.xyz'))
                    else:
                        shutil.copy(Path(f'{new_molecule_directory}/tmp.xyz'),
                                    Path(f'{new_molecule_directory}/{site}/mm.xyz'))
                    shutil.copy(Path(f'{new_molecule_directory}/am1.dat'),
                                Path(f'{new_molecule_directory}/{site}/am1.dat'))
                    # os.system('cp ../am1.dat ../mm.xyz .')
                    dispatcher_ts[radical](Path(f'{new_molecule_directory}/{site}/am1.dat'),
                                           Path(f'{new_molecule_directory}/{site}/mm.xyz'), site,
                                           Path(f'{new_molecule_directory}/{site}/'),
                                           pos_charged, False, True, new_mopac_2016)

                    if os.path.isfile(Path(f'{new_molecule_directory}/{site}/am1.arc')) is True:
                        ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/am1.out'),
                                             Path(f'{new_molecule_directory}/{site}/ts.xyz'),
                                             Path(f'{new_molecule_directory}/{site}'))
                        if ar.mopac_distance_checks(Path(f'{new_molecule_directory}/{site}/ts.xyz'), f'{radical}',
                                                    site) is False:

                            os.remove(Path(f'{new_molecule_directory}/{site}/am1.out'))
                            os.remove(Path(f'{new_molecule_directory}/{site}/am1.dat'))
                            check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/{site}/mm.xyz'), '-omopin',
                                        '-O', Path(f'{new_molecule_directory}/{site}/am1.dat')], stdout=DEVNULL,
                                       stderr=STDOUT)
                            dispatcher_ts[radical](Path(f'{new_molecule_directory}/{site}/am1.dat'),
                                                   Path(f'{new_molecule_directory}/{site}/mm.xyz'), site,
                                                   Path(f'{new_molecule_directory}/{site}/'),
                                                   pos_charged, True, False, new_mopac_2016)
                            if os.path.isfile(Path(f'{new_molecule_directory}/{site}/am1.arc')) is True:

                                ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/am1.out'),
                                                     Path(f'{new_molecule_directory}/{site}/constrained.xyz'),
                                                     Path(f'{new_molecule_directory}/{site}'))
                                if ar.mopac_distance_checks(Path(f'{new_molecule_directory}/{site}/constrained.xyz'),
                                                            f'{radical}', site) is True:
                                    check_call(['obabel', '-ixyz',
                                                Path(f'{new_molecule_directory}/{site}/constrained.xyz'), '-omopin',
                                                '-O', Path(f'{new_molecule_directory}/{site}/ts2.dat')],
                                               stdout=DEVNULL, stderr=STDOUT)
                                    with open(Path(f'{new_molecule_directory}/{site}/ts2.dat'), 'r+') as relaxed:
                                        lines = relaxed.readlines()
                                        if new_mopac_2016 is False:
                                            lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF CYCLES=10000\n'
                                        else:
                                            lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF CYCLES=10000\n'

                                    os.remove(Path(f'{new_molecule_directory}/{site}/ts2.dat'))
                                    with open(Path(f'{new_molecule_directory}/{site}/ts2.dat'), 'w+') as file:
                                        for line in lines:
                                            file.write(str(line))
                                    check_call([mopac_path, Path(f'{new_molecule_directory}/{site}/ts2')],
                                               stdout=DEVNULL,
                                               stderr=STDOUT)
                                    if os.path.isfile(Path(f'{new_molecule_directory}/{site}/ts2.arc')) is True:

                                        ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/ts2.out'),
                                                             Path(f'{new_molecule_directory}/{site}/ts.xyz'),
                                                             Path(f'{new_molecule_directory}/{site}'))
                                        ts_converged[site] = ar.mopac_freq_check(
                                            Path(f'{new_molecule_directory}/{site}/ts.xyz'), pos_charged, radical, site,
                                            Path(f'{new_molecule_directory}/{site}'),
                                            new_mopac_2016, multiple=False)

                                    else:
                                        ar.tweak_distance(Path(f'{new_molecule_directory}/{site}/ts2.dat'),
                                                          f'{radical}', Path(f'{new_molecule_directory}/{site}/'))
                                        if os.path.isfile(Path(f'{new_molecule_directory}/{site}/ts2.arc')):
                                            ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/ts2.out'),
                                                                 Path(f'{new_molecule_directory}/{site}/ts.xyz'),
                                                                 Path(f'{new_molecule_directory}/{site}'))
                                            ts_converged[site] = ar.mopac_freq_check(
                                                Path(f'{new_molecule_directory}/{site}/ts.xyz'), pos_charged, radical,
                                                site, Path(f'{new_molecule_directory}/{site}'),
                                                new_mopac_2016, multiple=False)

                                else:
                                    ts_converged[site] = False
                                    print('Did not converge')
                            else:
                                ts_converged[site] = False
                                print('Did not converge')
                        else:
                            ts_converged[site] = ar.mopac_freq_check(Path(f'{new_molecule_directory}/{site}/ts.xyz'),
                                                                     pos_charged, radical, site,
                                                                     Path(f'{new_molecule_directory}/{site}'),
                                                                     new_mopac_2016, multiple=False)

                    else:
                        os.remove(Path(f'{new_molecule_directory}/{site}/am1.out'))
                        os.remove(Path(f'{new_molecule_directory}/{site}/am1.dat'))
                        check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/{site}/mm.xyz'), '-omopin', '-O',
                                    Path(f'{new_molecule_directory}/{site}/am1.dat')], stdout=DEVNULL,
                                   stderr=STDOUT)
                        dispatcher_ts[radical](Path(f'{new_molecule_directory}/{site}/am1.dat'),
                                               Path(f'{new_molecule_directory}/{site}/mm.xyz'), site,
                                               Path(f'{new_molecule_directory}/{site}/'),
                                               pos_charged, True, False, new_mopac_2016)
                        if os.path.isfile(Path(f'{new_molecule_directory}/{site}/am1.arc')) is True:

                            ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/am1.out'),
                                                 Path(f'{new_molecule_directory}/{site}/constrained.xyz'),
                                                 Path(f'{new_molecule_directory}/{site}'))
                            if ar.mopac_distance_checks(Path(f'{new_molecule_directory}/{site}/constrained.xyz'),
                                                        f'{radical}', site) is True:
                                check_call(['obabel', '-ixyz', Path(f'{new_molecule_directory}/{site}/constrained.xyz'),
                                            '-omopin', '-O', Path(f'{new_molecule_directory}/{site}/ts2.dat')],
                                           stdout=DEVNULL, stderr=STDOUT)
                                with open(Path(f'{new_molecule_directory}/{site}/ts2.dat'), 'r+') as relaxed:
                                    lines = relaxed.readlines()
                                    if new_mopac_2016 is False:
                                        lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF CYCLES=10000\n'
                                    else:
                                        lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF CYCLES=10000\n'

                                os.remove(Path(f'{new_molecule_directory}/{site}/ts2.dat'))
                                with open(Path(f'{new_molecule_directory}/{site}/ts2.dat'), 'w+') as file:
                                    for line in lines:
                                        file.write(str(line))
                                check_call([mopac_path, Path(f'{new_molecule_directory}/{site}/ts2')], stdout=DEVNULL,
                                           stderr=STDOUT)
                                if os.path.isfile(Path(f'{new_molecule_directory}/{site}/ts2.arc')) is True:

                                    ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/ts2.out'),
                                                         Path(f'{new_molecule_directory}/{site}/ts.xyz'),
                                                         Path(f'{new_molecule_directory}/{site}'))
                                    ts_converged[site] = ar.mopac_freq_check(
                                        Path(f'{new_molecule_directory}/{site}/ts.xyz'), pos_charged, radical, site,
                                        Path(f'{new_molecule_directory}/{site}'),
                                        new_mopac_2016, multiple=False)

                                else:
                                    ar.tweak_distance(Path(f'{new_molecule_directory}/{site}/ts2.dat'), f'{radical}',
                                                      Path(f'{new_molecule_directory}/{site}/'))
                                    if os.path.isfile(Path(f'{new_molecule_directory}/{site}/ts2.arc')):
                                        ar.read_mopac_output(Path(f'{new_molecule_directory}/{site}/ts2.out'),
                                                             Path(f'{new_molecule_directory}/{site}/ts.xyz'),
                                                             Path(f'{new_molecule_directory}/{site}'))
                                        ts_converged[site] = ar.mopac_freq_check(
                                            Path(f'{new_molecule_directory}/{site}/ts.xyz'), pos_charged, radical, site,
                                            Path(f'{new_molecule_directory}/{site}'),
                                            new_mopac_2016, multiple=False)

                            else:
                                ts_converged[site] = False
                                print('Did not converge')

                        else:
                            ts_converged[site] = False
                            print('Did not converge')
            hf_jobs = str(Path(f'{str(Path("/")).join(new_molecule_directory.split(str(Path("/")))[:-1])}/hf_jobs.txt'))
            print(f'HPC CALCS ARE {hpc_calcs}')
            if hpc_calcs:
                hf.run_pyreact(f'-d a --nwchem --mult 2 --FC --NoEigen -F UHF --nProc 8 --TS --Rega --array {hf_jobs}')
            if compound_yes_rad_no is True:
                molecule_dict = sd.store_data(compound_name, canonical_smiles, aromatic_atom, radical,
                                              reagent_optimised, ts_converged, Path(f'{new_molecule_directory}/'),
                                              hpc_calcs, compounds_dict[compound_name])
            else:
                molecule_dict = sd.store_data(compound_name, canonical_smiles, aromatic_atom, radical,
                                              reagent_optimised, ts_converged, Path(f'{new_molecule_directory}/'),
                                              hpc_calcs)

            compounds_dict[str(compound_name)] = molecule_dict
            with open(f'{new_molecule_directory}_results.json', 'w') as json_file:
                json.dump(compounds_dict, json_file)
            with open(f'{new_molecule_directory}_results.csv', 'w') as csvfile:
                if hpc_calcs:
                    header = ['compound_name', 'smiles', 'radical', 'Site', 'hf_reactant_energy', 'hf_ts_energy',
                          'hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0']
                else:
                    header = ['compound_name', 'smiles', 'radical', 'Site', 'am1_reactant_energy', 'am1_ts_energy',
                          'am1_activation_energy', 'am1_ratio']
                writer = csv.writer(csvfile)
                writer.writerow(header)
                for com in compounds_dict.keys():
                    for radical in compounds_dict[com]['radicals'].values():
                        for key, value in radical.items():
                            if str.isdigit(key):
                                if hpc_calcs:
                                    row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                           radical['radical'], key,
                                           radical['hf_reactant_energy'], radical[key]['hf_ts_energy'],
                                           radical[key]['hf_act_energy'], radical[key]['hf_ratio'], radical[key]['fa+'],
                                           radical[key]['fa-'], radical[key]['fa0']]
                                else:
                                    row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                           radical['radical'], key,
                                           radical['am1_reactant_energy'], radical[key]['am1_ts_energy'],
                                           radical[key]['am1_act_energy'], radical[key]['am1_ratio']]
                                writer.writerow(row)
        else:
            with open(f'{new_molecule_directory}_results.json', 'r+') as completed_json:
                compounds_dict = json.load(completed_json)
            with open(f'{new_molecule_directory}_results.csv', 'w') as csvfile:
                header = ['compound_name', 'smiles', 'radical', 'Site', 'hf_reactant_energy', 'hf_ts_energy',
                          'hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0']
                writer = csv.writer(csvfile)
                writer.writerow(header)
                for com in compounds_dict.keys():
                    for radical in compounds_dict[com]['radicals'].values():
                        for key, value in radical.items():
                            if str.isdigit(key):
                                row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                       radical['radical'], key,
                                       radical['hf_reactant_energy'], radical[key]['hf_ts_energy'],
                                       radical[key]['hf_act_energy'], radical[key]['hf_ratio'], radical[key]['fa+'],
                                       radical[key]['fa-'], radical[key]['fa0']]
                                writer.writerow(row)

    elif previously_calculated == 3:
        new_molecule_directory = previously_calculated_compound_dir
        print('ADDING HPC CALCS')
        print(f'AROMATIC ATOM IS {aromatic_atom}')
        if hpc_calcs:
            hf_jobs = str(Path(f'{str(Path("/")).join(new_molecule_directory.split(str(Path("/")))[:-1])}/hf_jobs.txt'))
            if os.path.isfile(hf_jobs) is True:
                with open(hf_jobs, 'r') as fil:
                    lins = fil.readlines()
                for li in lins:
                    li = str(Path(li)).split(str(Path('/')))
                    print(f'LI[-4] IS {li[-4]}')
                    # removes the previous hf_jobs.txt file if other compounds calculated in the directory before.
                    if new_molecule_directory.split(str(Path("/")))[-1] != li[-3]:
                        print('REMOVING FILE')
                        os.remove(hf_jobs)
                        break
                with open(hf_jobs, 'a') as fi:
                    for site in aromatic_atom:
                        print(f'NEW MOLECULE DIRECTORY {new_molecule_directory}')
                        print(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                        if os.path.isfile(Path(f'{new_molecule_directory}/{site}/hf.sdf')):
                            fi.write(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                    if os.path.isfile(Path(f'{new_molecule_directory}/reagent/hf.sdf')):
                        fi.write(str(Path(f'{new_molecule_directory}/reagent/hf.sdf\n')))
                    fi.write(str(Path(f'{new_molecule_directory}/fukui/hf.sdf\n')))
            else:
                with open(hf_jobs, 'w') as fil:
                    for site in aromatic_atom:
                        print(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                        if os.path.isfile(Path(f'{new_molecule_directory}/{site}/hf.sdf')):
                            fil.write(str(Path(f'{new_molecule_directory}/{site}/hf.sdf\n')))
                    if os.path.isfile(Path(f'{new_molecule_directory}/reagent/hf.sdf')):
                        fil.write(str(Path(f'{new_molecule_directory}/reagent/hf.sdf\n')))
                    fil.write(str(Path(f'{new_molecule_directory}/fukui/hf.sdf\n')))
            hf.run_pyreact(f'-d a --nwchem --mult 2 --FC --NoEigen -F UHF --nProc 8 --TS --Rega --array {hf_jobs}')
            if compound_yes_rad_no is True:
                molecule_dict = sd.store_data(compound_name, canonical_smiles, aromatic_atom, radical,
                                              reagent_optimised, ts_converged, Path(f'{new_molecule_directory}/'),
                                              hpc_calcs, compounds_dict[compound_name])
            else:
                molecule_dict = sd.store_data(compound_name, canonical_smiles, aromatic_atom, radical,
                                              reagent_optimised, ts_converged, Path(f'{new_molecule_directory}/'),
                                              hpc_calcs)

            compounds_dict[str(compound_name)] = molecule_dict
            with open(f'{new_molecule_directory}_results.json', 'w') as json_file:
                json.dump(compounds_dict, json_file)
            with open(f'{new_molecule_directory}_results.csv', 'w') as csvfile:
                if hpc_calcs:
                    header = ['compound_name', 'smiles', 'radical', 'Site', 'hf_reactant_energy', 'hf_ts_energy',
                              'hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0']
                else:
                    header = ['compound_name', 'smiles', 'radical', 'Site', 'am1_reactant_energy', 'am1_ts_energy',
                              'am1_activation_energy', 'am1_ratio']
                writer = csv.writer(csvfile)
                writer.writerow(header)
                for com in compounds_dict.keys():
                    for radical in compounds_dict[com]['radicals'].values():
                        for key, value in radical.items():
                            if str.isdigit(key):
                                if hpc_calcs:
                                    row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                           radical['radical'], key,
                                           radical['hf_reactant_energy'], radical[key]['hf_ts_energy'],
                                           radical[key]['hf_act_energy'], radical[key]['hf_ratio'], radical[key]['fa+'],
                                           radical[key]['fa-'], radical[key]['fa0']]
                                else:
                                    row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                           radical['radical'], key,
                                           radical['am1_reactant_energy'], radical[key]['am1_ts_energy'],
                                           radical[key]['am1_act_energy'], radical[key]['am1_ratio']]
                                writer.writerow(row)
        else:
            with open(f'{new_molecule_directory}_results.json', 'w') as json_file:
                json.dump(compounds_dict, json_file)
            with open(f'{new_molecule_directory}_results.csv', 'w') as csvfile:
                header = ['compound_name', 'smiles', 'radical', 'Site', 'am1_reactant_energy', 'am1_ts_energy',
                          'am1_activation_energy', 'am1_ratio']
                writer = csv.writer(csvfile)
                writer.writerow(header)
                for com in compounds_dict.keys():
                    for radical in compounds_dict[com]['radicals'].values():
                        for key, value in radical.items():
                            if str.isdigit(key):
                                row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                       radical['radical'], key,
                                       radical['am1_reactant_energy'], radical[key]['am1_ts_energy'],
                                       radical[key]['am1_act_energy'], radical[key]['am1_ratio']]
                                writer.writerow(row)
    else:
        new_molecule_directory = previously_calculated_compound_dir
        with open(f'{new_molecule_directory}_results.json', 'w') as json_file:
            json.dump(compounds_dict, json_file)
        with open(f'{new_molecule_directory}_results.csv', 'w') as csvfile:
            header = ['compound_name', 'smiles', 'radical', 'Site', 'hf_reactant_energy', 'hf_ts_energy',
                      'hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0']
            writer = csv.writer(csvfile)
            writer.writerow(header)
            for com in compounds_dict.keys():
                for radical in compounds_dict[com]['radicals'].values():
                    for key, value in radical.items():
                        if str.isdigit(key):
                            row = [compounds_dict[com]['compound_name'], compounds_dict[com]['smiles'],
                                   radical['radical'], key,
                                   radical['hf_reactant_energy'], radical[key]['hf_ts_energy'],
                                   radical[key]['hf_act_energy'], radical[key]['hf_ratio'], radical[key]['fa+'],
                                   radical[key]['fa-'], radical[key]['fa0']]
                            writer.writerow(row)
    return compounds_dict, aromatic_atom
