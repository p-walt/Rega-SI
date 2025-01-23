from func_calc.store_data import store_data as sd
from pathlib import Path
import os
import json
# TODO RA: Make this a function. 
with open('all_calcs1.txt', 'r') as all_calcs:
    alines1 = all_calcs.readlines()
    aline2 = []
    dir_list = []
    dir_dict = {}
    results_dict = {}
    directory = Path.cwd()
    for fi in os.listdir(directory):
        d = os.path.join(directory, fi)
        if os.path.isdir(d):
            dir_list.append(d)
    for aline1 in alines1:
        st = aline1.strip('\n')
        aline2.append(st)
        # TODO RA: What is this doing?
    with open('sites_complete.txt', 'r') as comp_calcs:
        clines1 = comp_calcs.readlines()
        cline2 = []
        cline3 = []
        for cline1 in clines1:
            ct = cline1.strip('\n')
            ct_split = ct.split(str(Path('/')))
            cline3.append('/'.join(ct_split[1:]))
            cline2.append(ct) # TODO RA: Avoid using .append. It is incredibly slow and so should be avaoided at all costs. Its' less of a problem if you dont do it many times but it tends to be bad in for loops. Consider fstrings instead? 

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
    print(compounds_dict)
    for molecule in compounds_di.keys():
        filelist = os.listdir(molecule)
        previous_dict_path = Path(f'{str(Path("/")).join(molecule.split(str(Path("/")))[:-1])}/{molecule.split(str(Path("/")))[-1]}_results.json')
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
        print(f'PREV MOL DICT {previous_dict}')
        molecule_dict = sd.store_data(compound_name=mol_name, smiles=mol_smi, compound_sites=mol_sites, radical='cf3',
                                      reactant_optimised=True, ts_converged=compounds_di[molecule],
                                      directory=str(Path("/")).join(molecule.split(str(Path("/")))[:-1]), hpc_calcs=True,
                                      prev_molecule_dict=previous_dict)
        results_dict[str(mol_name)] = molecule_dict
print(results_dict)
# with open('test_results2.json', 'w') as results:
#     json.dump(results_dict, results)

    # print(compounds_di)
        # site_converged_dict = {}
        # site_list = []
        # for fi in os.listdir(molecule_dir):
        #     d = os.path.join(molecule_dir, fi)
        #     if os.path.isdir(d):
        #         if isinstance(d.split(str(Path('/')))[-1], int):
        #             site_list.append(d.split(str(Path('/')))[-1])
        # for site in site_list:
        #     if f'{molecule_dir}' in ts_converged_dict.keys():
        #         if Path(f'{molecule_dir}/{site}/hfoutp.in') in lines:
        #             ts_converged_dict[molecule_dir][site] = True
        #         else:
        #             ts_converged_dict[molecule_dir][site] = False
