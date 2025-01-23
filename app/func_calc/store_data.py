import os
import re
import numpy as np
from pathlib import Path
from math import log10, floor, exp
from subprocess import check_call, STDOUT, DEVNULL
import chemcoord as cc


# def round_it(x, sig):
#     return round(x, sig - int(floor(log10(abs(x)))) - 1)


def findlargest(arr):
    # TODO RA: docs?
    secondlargest = arr[0]
    largest = arr[0]
    for i in range(len(arr)):
        if arr[i] > largest:
            largest = arr[i]

    for i in range(len(arr)):
        if arr[i] > secondlargest and arr[i] != largest:
            secondLargest = arr[i] # TODO RA: This will break

    return secondlargest


def store_data(compound_name, smiles, compound_sites, radical, reactant_optimised,
                   ts_converged, directory, hpc_calcs, prev_molecule_dict=None, product_smiles_dict=None):
    """
    Function that grabs the energies from the HF and AM1 calcualtions and stores them in the data dictionary.
    Also calculates the ratio of products for each compound based on difference between AM1 or HF activation energies.
        @param compound_name: Name in molecules_list.csv.
        @param smiles: The compound of interest SMILES string.
        @param compound_sites: Site of reaction in the molecule in question.
        @param radical: The radical of interest.
        @param reactant_optimised: Whether the reactant successfully optimised or not.
        @param ts_converged: Whether the transition state was actually a true AM1 transition state (checked through freq
        calculation) for each site.
        @param directory: Compound Directory.
        @param prev_molecule_dict: If calculated previously the molecule dictionary is passed into the function and the
        new radical data is added.
        @return: Data dictionary."""# TODO RA: types?
    fplus_dict = {}
    fmin_dict = {}
    f0_dict = {}
    gross_pop_dict = {}
    bond_indices_dict = {}
    mulliken_dict = {}
    esp_dia_dict = {}
    elec_density_dict = {}
    # TODO RA: Hard coded values...

    if radical == 'cf3':
        radical_reactant_energy = -1583.81910
    else:
        pass
    if not reactant_optimised:
        if radical == 'cf2h':
            radical_reactant_energy = -1112.03737
        elif radical == 'ipr':
            radical_reactant_energy = -480.12463
    elif reactant_optimised and radical != 'cf3':
        radical_reactant_energy = 0.0
    if prev_molecule_dict is None:
        molecule_dict = {'compound_name': compound_name, 'smiles': smiles, 'radicals': {}}
    else:
        molecule_dict = prev_molecule_dict
    radical_dict = {'radical': radical}
    if hpc_calcs:
        try:
            with open(Path(f'{directory}/reagent/hf2outp.out'), encoding='utf-8') as reag:
                for line in reversed(reag.readlines()):
                    lines = line.strip('\n')
                    if 'Total SCF energy' in lines:
                        en = lines.split(' ')
                        hf_energy = float(en[-1])
                        break
                    else:
                        pass
            try:
                radical_dict.update({'hf_reactant_energy': hf_energy})
            except UnboundLocalError:
                energy = 'No Reagent Energy Found'
                radical_dict.update({'hf_reactant_energy': energy})
        except FileNotFoundError:
            try:
                with open(Path(f'{directory}/reagent/hfoutp.out'), encoding='utf-8') as reag:
                    for line in reversed(reag.readlines()):
                        lines = line.strip('\n')
                        if 'Total SCF energy' in lines:
                            en = lines.split(' ')
                            hf_energy = float(en[-1])
                            break
                        else:
                            pass
                try:
                    radical_dict.update({'hf_reactant_energy': hf_energy})
                except UnboundLocalError:
                    energy = 'No Reagent Energy Found'
                    radical_dict.update({'hf_reactant_energy': energy})
            except FileNotFoundError:
                energy = 'No Reagent Energy Found'
                radical_dict.update({'hf_reactant_energy': energy})
    try:
        with open(Path(f'{directory}/reagent/am1.arc')) as reag:
            for line in reag:
                lines = line.strip('\n')
                if 'TOTAL ENERGY' in lines:
                    energy = lines
                else:
                    pass
        reactant_e = float(re.split(' +', str(energy))[4]) + radical_reactant_energy
        print(reactant_e)
        radical_dict.update({'am1_reactant_energy': reactant_e})
    except FileNotFoundError:
        e = 'No TS Found'
        radical_dict.update({'am1_reactant_energy': e})
    if not compound_sites:
        molecule_dict['No Sites of reaction found'] = 'No sites of reaction found'
        return molecule_dict
    if hpc_calcs:
        try:
            with open(Path(f'{directory}/fukui/hf2outp.out'), encoding='utf-8') as rf, open(Path(f'{directory}/fukui/out_rev.txt'),
                                                                              'w', encoding='utf-8') as wf:
                for line in reversed(rf.readlines()):
                    wf.write(line)
        except FileNotFoundError:
            try:
                with open(Path(f'{directory}/fukui/hfoutp.out'), encoding='utf-8') as rf, open(
                        Path(f'{directory}/fukui/out_rev.txt'),
                        'w', encoding='utf-8') as wf:
                    for line in reversed(rf.readlines()):
                        wf.write(line)
            except FileNotFoundError:
                pass
            a = []
            try:
                with open(Path(f'{directory}/fukui/out_rev.txt'), encoding='utf-8') as file:
                    for line in file:
                        lines = line.strip('\n')
                        if 'Condensed Fukui function [fsn(+)]' in lines:
                            # collect block-related lines
                            while True:
                                try:
                                    lines = next(file)

                                except StopIteration:
                                    # there is no lines left
                                    break
                                if 'Condensed Fukui function [fnn(+)]' in lines:
                                    # we've reached the end of block
                                    break
                                a.append(lines)
                            # stop iterating over file
                            break

                alla = a[2:-4]
                i = 0
                for el in alla:
                    if el == ' -----------   ----------------\n':
                        fminus = list(reversed(alla[0:i]))
                        break
                    i += 1
                j = 0
                alla.reverse()
                try:
                    for atom in fminus:
                        ato = atom.strip('\n')
                        line = ato.split()
                        fmin_dict[line[0]] = line[3]
                except NameError:
                    pass

                for e in alla:
                    if e == '\n':
                        fplus = list(alla[0:j])
                        break
                    j += 1

                try:
                    for atom in fplus:
                        ato = atom.strip('\n')
                        line = ato.split()
                        fplus_dict[line[0]] = line[3]
                except NameError:
                    pass
                try:
                    for atom in fmin_dict.keys():
                        fmi = float(fmin_dict[atom])
                        fplu = float(fplus_dict[atom])
                        f0 = 0.5 * (fplu + fmi)
                        f0_dict[atom] = f0
                except KeyError:
                    pass
            except FileNotFoundError:
                for atomy in compound_sites:
                    fplus_dict[atomy] = 'NOT FOUND'
                    fmin_dict[atomy] = 'NOT FOUND'
                    f0_dict[atomy] = 'NOT FOUND'
        # print(fplus_dict)
        # print(fmin_dict)
        # print(f0_dict)
        try:
            with open(Path(f'{directory}/property/hfoutp.out'), encoding='utf-8') as rf, open(Path(f'{directory}/property/out_rev.txt'), 'w', encoding='utf-8') as wf:
                for line in reversed(rf.readlines()):
                    wf.write(line)
            a = []
            with open(Path(f'{directory}/property/out_rev.txt'), encoding='utf-8') as file:
                lines1 = iter(file.readlines())

                for line in lines1:
                    lines = line.strip('\n')
                    if 'Total spin density' in lines:
                        # collect block-related lines
                        while True:
                            try:
                                lines = next(lines1)

                            except StopIteration:
                                # there is no lines left
                                break
                            if '----------------------' in lines:
                                # we've reached the end of block
                                break
                            a.append(lines)
                        # stop iterating over file
                        break
                spin_density = list(reversed(a[1:-1]))
                elec_density_dict = {}
                for spind in spin_density:
                    spind.strip('\n')
                    line_elements1 = spind.split()
                    elec_density_dict[line_elements1[0]] = line_elements1[5]
                b = []
                esp_dia_dict = {}
                for line in lines1:
                    lines = line.strip('\n')
                    if 'scftyp:UHF' in lines:
                        while True:
                            try:
                                lines = next(lines1)
                            except StopIteration:
                                break
                            if 'Electrostatic potential/diamagnetic shielding' in lines:
                                break
                            b.append(lines)
                raw_data = list(reversed(b[0:-4]))
                for data in raw_data:
                    data.strip('\n')
                    data_elements = data.split()
                    esp_dia_dict[data_elements[0]] = [data_elements[5], data_elements[6]]
                # value index 0 is electrostatic potential and index 1 is diamagnetic shielding
            with open(Path(f'{directory}/property/out_rev.txt'), encoding='utf-8') as file:
                lines2 = iter(file.readlines())
                c = []
                mulliken_dict = {}

                for line in lines2:
                    lines = line.strip('\n')
                    if 'Electrostatic potential/diamagnetic shielding' in lines:
                        while True:
                            try:
                                lines = next(lines2)
                            except StopIteration:
                                break
                            if 'Valency' in lines:
                                break
                            c.append(lines)

                raw_data2 = list(reversed(c[2:]))
                for data2 in raw_data2:
                    data2.strip('\n')
                    data_elements2 = data2.split()

                    mulliken_dict[data_elements2[0]] = data_elements2[5]
            with open(Path(f'{directory}/property/out_rev.txt'), encoding='utf-8') as file:
                lines3 = iter(file.readlines())
                d = []
                bond_indices_dict = {}

                for line in lines3:
                    lines = line.strip('\n')
                    if 'Valency' in lines:
                        while True:
                            try:
                                lines = next(lines3)
                            except StopIteration:
                                break
                            if 'Large bond indices' in lines:
                                break
                            d.append(lines)

                raw_data3 = list(reversed(d[3:-1]))
                for data3 in raw_data3:
                    data3.strip('\n')
                    data_elements3 = data3.split()
                    # Get bond indices focusing on C-H bonds. Site number of aromatic C-H will then be used to get correct bond index.
                    if data_elements3[1] == 'C' and data_elements3[4] == 'H':
                        bond_indices_dict[data_elements3[0]] = data_elements3[5]

            with open(Path(f'{directory}/property/out_rev.txt'), encoding='utf-8') as file:
                lines4 = iter(file.readlines())
                e = []
                gross_pop_dict = {}

                for line in lines4:
                    lines = line.strip('\n')
                    if 'Atomic spin population' in lines:
                        while True:
                            try:
                                lines = next(lines4)
                            except StopIteration:
                                break
                            if 'Total      gross population on atoms' in lines:
                                break
                            e.append(lines)

                raw_data4 = list(reversed(e[1:-1]))
                for data4 in raw_data4:
                    data4.strip('\n')
                    data_elements4 = data4.split()
                    gross_pop_dict[data_elements4[0]] = data_elements4[3]

        except FileNotFoundError:
            gross_pop_dict = {}
            bond_indices_dict = {}
            mulliken_dict = {}
            esp_dia_dict = {}
            elec_density_dict = {}

    for site in compound_sites:
        site_dict = {'site': site}
        if product_smiles_dict is not None:
            try:
                site_dict['product_smiles'] = product_smiles_dict[site]
            except KeyError:
                site_dict['product_smiles'] = 'NO PRODUCT SMILE FOUND'
        if hpc_calcs:
            try:
                with open(Path(f'{directory}/property/hfoutp.out'), 'r', encoding='utf-8') as file:
                    rev_file = reversed(file.readlines())
                    lines2 = iter(list(rev_file))
                    f = []
                    property_geom = []
                    property_geometry_ch_length_dict = {}

                    for line in lines2:
                        lines = line.strip('\n')
                        if 'Atomic Mass' in lines:
                            while True:
                                try:
                                    lines = next(lines2)
                                except StopIteration:
                                    break
                                if '---- ----------------' in lines:
                                    break

                                f.append(lines)
                            break
                    raw_data5 = list(reversed(f[1:]))
                    for data5 in raw_data5:
                        data5.strip('\n')
                        data_elements5 = data5.split()
                        atom = f'{data_elements5[1]}\t{data_elements5[3]}\t{data_elements5[4]}\t{data_elements5[5]}\n'
                        property_geometry_ch_length_dict[data_elements5[0]] = data_elements5[3]
                        property_geom.append(atom)
                    atom_count = len(property_geom)
                    with open('test.xyz', 'w') as new:
                        new.write(f'{atom_count}\n\n')
                        for ato in property_geom:
                            new.write(ato)
                    compound = cc.Cartesian.read_xyz('test.xyz', start_index=1)
                    connection_table = compound.get_bonds()
                    print(connection_table[int(site)])
                    for connection in connection_table[int(site)]:
                        if 'H' in property_geom[connection - 1]:
                            site_dict['HF_C-H_Bond_Length'] = compound.get_bond_lengths([int(site), connection])[0]
            except FileNotFoundError:
                gross_pop_dict = {}
                bond_indices_dict = {}
                mulliken_dict = {}
                esp_dia_dict = {}
                elec_density_dict = {}

            try:
                with open(Path(f'{directory}/{site}/hf2outp.out'), encoding='utf-8') as reag:
                    rev_file = reversed(reag.readlines())
                    lines = list(rev_file)
                    complete = False
                    for com in lines[0:5]:
                        co = com.strip('\n')
                        if 'Total times' in co:
                            complete = True
                    if complete:
                        for line in lines:
                            lin = line.strip('\n')
                            if 'Total SCF energy' in lin:
                                en = lin.split(' ')
                                hf_ts_energy = float(en[-1])
                                break
                            else:
                                pass
                    else:
                        hf_ts_energy = 'No TS Found'
                with open(Path(f'{directory}/{site}/hf2outp.out'), encoding='utf-8') as reag2:
                    rev_file2 = reversed(reag2.readlines())
                    lines2 = iter(list(rev_file2))
                    geometry = []
                    g = []
                    for line2 in lines2:
                        li = line2.strip('\n')
                        if 'Atomic Mass' in li:
                            print('found')
                            while True:
                                try:
                                    li = next(lines2)
                                except StopIteration:
                                    break
                                if '---- ----------------' in li:
                                    print('ended')
                                    break

                                g.append(li)
                            break
                    raw_data6 = list(reversed(g[1:]))
                    for data6 in raw_data6:
                        data6.strip('\n')
                        data_elements6 = data6.split()
                        atom1 = f'{data_elements6[1]}\t{data_elements6[3]}\t{data_elements6[4]}\t{data_elements6[5]}\n'
                        geometry.append(atom1)
                    atom_count1 = len(geometry)
                    with open(Path(f'{directory}/{site}/hf_converged.xyz'), 'w', encoding='utf-8') as new_geo:
                        new_geo.write(f'{atom_count1}\n\n')
                        for atoms in geometry:
                            new_geo.write(f'{atoms}')
                    check_call(['obabel', '-ixyz', str(Path(f'{directory}/{site}/hf_converged.xyz')), '-omopin', '-O',
                                str(Path(f'{directory}/{site}/hf_converged.dat'))], stdout=DEVNULL, stderr=STDOUT)
                    with open(str(Path(f'{directory}/{site}/hf_converged.dat')), 'r') as dat:
                        geom = dat.readlines()
                        radical_c = geom[-4].split()
                        if radical_c[7] == str(site):
                            site_dict['HF_TS_C-C_Bond_Length'] = radical_c[1]
                        elif radical_c[7] == int(site):
                            site_dict['HF_TS_C-C_Bond_Length'] = radical_c[1]
                        else:
                            site_dict['HF_TS_C-C_Bond_Length'] = 'NOT FOUND'
            except FileNotFoundError:
                try:
                    with open(Path(f'{directory}/{site}/hfoutp.out'), encoding='utf-8') as reag:
                        rev_file = reversed(reag.readlines())
                        lines = list(rev_file)
                        complete = False
                        for com in lines[0:5]:
                            co = com.strip('\n')
                            if 'Total times' in co:
                                complete = True
                        if complete:
                            for line in lines:
                                lin = line.strip('\n')
                                if 'Total SCF energy' in lin:
                                    en = lin.split(' ')
                                    hf_ts_energy = float(en[-1])
                                    break
                                else:
                                    pass
                        else:
                            hf_ts_energy = 'No TS Found'
                    with open(Path(f'{directory}/{site}/hfoutp.out'), encoding='utf-8') as reag2:
                        rev_file2 = reversed(reag2.readlines())
                        lines2 = iter(list(rev_file2))
                        geometry = []
                        g = []
                        for line2 in lines2:
                            li = line2.strip('\n')
                            if 'Atomic Mass' in li:
                                print('found')
                                while True:
                                    try:
                                        li = next(lines2)
                                    except StopIteration:
                                        break
                                    if '---- ----------------' in li:
                                        print('ended')
                                        break

                                    g.append(li)
                                break
                        raw_data6 = list(reversed(g[1:]))
                        for data6 in raw_data6:
                            data6.strip('\n')
                            data_elements6 = data6.split()
                            atom1 = f'{data_elements6[1]}\t{data_elements6[3]}\t{data_elements6[4]}\t{data_elements6[5]}\n'
                            geometry.append(atom1)
                        atom_count1 = len(geometry)
                        with open(Path(f'{directory}/{site}/hf_converged.xyz'), 'w', encoding='utf-8') as new_geo:
                            new_geo.write(f'{atom_count1}\n\n')
                            for atoms in geometry:
                                new_geo.write(f'{atoms}')
                        check_call(
                            ['obabel', '-ixyz', str(Path(f'{directory}/{site}/hf_converged.xyz')), '-omopin', '-O',
                             str(Path(f'{directory}/{site}/hf_converged.dat'))], stdout=DEVNULL, stderr=STDOUT)
                        with open(str(Path(f'{directory}/{site}/hf_converged.dat')), 'r') as dat:
                            geom = dat.readlines()
                            radical_c = geom[-4].split()
                            print(f'SITE IS {site}')
                            print(f'RADICAL_C IS {radical_c[7]}')
                            if radical_c[7] == str(site):
                                site_dict['HF_TS_C-C_Bond_Length'] = radical_c[1]
                            elif radical_c[7] == int(site):
                                site_dict['HF_TS_C-C_Bond_Length'] = radical_c[1]
                            else:
                                site_dict['HF_TS_C-C_Bond_Length'] = 'NOT FOUND'

                except FileNotFoundError:
                    hf_ts_energy = 'No TS Found'
                    site_dict['HF_TS_C-C_Bond_Length'] = 'NOT FOUND'
        if ts_converged[site] is True:
            if os.path.isfile(Path(f'{directory}/{site}/ts2.arc')):
                with open(Path(f'{directory}/{site}/ts2.arc')) as reag:
                    for line in reag:
                        lines = line.strip('\n')
                        if 'TOTAL ENERGY' in lines:
                            am1_ts_energy = lines
                        else:
                            pass

            elif os.path.isfile(Path(f'{directory}/{site}/ts.arc')):
                with open(Path(f'{directory}/{site}/ts.arc')) as reag:
                    for line in reag:
                        lines = line.strip('\n')
                        if 'TOTAL ENERGY' in lines:
                            am1_ts_energy = lines
                        else:
                            pass

            elif os.path.isfile(Path(f'{directory}/{site}/am1.arc')):
                with open(Path(f'{directory}/{site}/am1.arc')) as reag:
                    for line in reag:
                        lines = line.strip('\n')
                        if 'TOTAL ENERGY' in lines:
                            am1_ts_energy = lines
                        else:
                            pass
        else:
            am1_ts_energy = '          TOTAL ENERGY            =      No TS Found EV'
        try:
            am1_ts_e = float(re.split(' +', str(am1_ts_energy))[4])
            # print(f'AM1 TS ENERGY {am1_ts_e}')
            site_dict['am1_ts_energy'] = am1_ts_e
            if reactant_e == 'No TS Found':
                site_dict['am1_act_energy'] = 'No TS Found'
            else:
                # print(round((am1_ts_e - reactant_e) * 23.061, 1))
                site_dict['am1_act_energy'] = round((float(am1_ts_e) - float(reactant_e)) * 23.061, 1)
                print(f'AM1 ACT ENERGY {site_dict["am1_act_energy"]}')
            # if os.path.isfile(Path(f'{directory}/{site}/ts.xyz')):
            #     with open(Path(f'{directory}/{site}/ts.xyz'), 'r+') as file:
            #         fi = file.readlines()
            #         cartesian_coordinates = fi[2:]
            #         system_size = fi[0].strip('\n')
            #         site_dict['cartesian_coordinates'] = cartesian_coordinates
            #         site_dict['system_size'] = system_size
            # else:
            #     pass
            radical_dict[str(site)] = site_dict
            radicals_dict = molecule_dict['radicals']
            radicals_dict.update({str(radical): radical_dict})
            molecule_dict.update({'radicals': radicals_dict})
        except ValueError:
            am1_ts_e = 'No TS Found'
            site_dict['am1_ts_energy'] = am1_ts_e

            site_dict['am1_act_energy'] = 'No TS Found'
            radical_dict[str(site)] = site_dict
            radicals_dict = molecule_dict['radicals']
            radicals_dict.update({str(radical): radical_dict})
            molecule_dict.update({'radicals': radicals_dict})
        except UnboundLocalError:
            am1_ts_e = 'No TS Found'
            site_dict['am1_ts_energy'] = am1_ts_e

            site_dict['am1_act_energy'] = 'No TS Found'
            radical_dict[str(site)] = site_dict
            radicals_dict = molecule_dict['radicals']
            radicals_dict.update({str(radical): radical_dict})
            molecule_dict.update({'radicals': radicals_dict})
        if hpc_calcs:
            try:
                hf_ts_e = float(hf_ts_energy)
                site_dict['hf_ts_energy'] = hf_ts_e
                try:
                    site_dict['fa+'] = fplus_dict[str(site)]
                    site_dict['fa-'] = fmin_dict[str(site)]
                    site_dict['fa0'] = f0_dict[str(site)]
                    site_dict['electron_density'] = elec_density_dict[str(site)]
                    site_dict['electrostatic_potential'] = esp_dia_dict[str(site)][0]
                    site_dict['diamagnetic_shielding'] = esp_dia_dict[str(site)][1]
                    site_dict['mulliken_charge'] = mulliken_dict[str(site)]
                    site_dict['bond_index'] = bond_indices_dict[str(site)]
                    site_dict['gross_population'] = gross_pop_dict[str(site)]

                except KeyError:
                    site_dict['fa+'] = 'N/A'
                    site_dict['fa-'] = 'N/A'
                    site_dict['fa0'] = 'N/A'
                    site_dict['electron_density'] = 'N/A'
                    site_dict['electrostatic_potential'] = 'N/A'
                    site_dict['diamagnetic_shielding'] = 'N/A'
                    site_dict['mulliken_charge'] = 'N/A'
                    site_dict['bond_index'] = 'N/A'
                    site_dict['gross_population'] = 'N/A'
                if hf_ts_energy == 'No TS Found':
                    site_dict['hf_act_energy'] = 'No TS Found'
                else:
                    site_dict['hf_act_energy'] = round((hf_ts_e - hf_energy) * 627.5, 1)
                print(Path(f'{directory}/{site}/ts.xyz'))
                if os.path.isfile(Path(f'{directory}/{site}/hf_ts.xyz')):
                    with open(Path(f'{directory}/{site}/hf_ts.xyz'), 'r+') as file:
                        fi = file.readlines()
                        cartesian_coordinates = fi[2:]
                        system_size = fi[0].strip('\n')
                        site_dict['hf_cartesian_coordinates'] = cartesian_coordinates
                        site_dict['system_size'] = system_size

                radical_dict[str(site)] = site_dict
                radicals_dict = molecule_dict['radicals']
                radicals_dict.update({str(radical): radical_dict})
                molecule_dict.update({'radicals': radicals_dict})
            except ValueError:
                hf_ts_e = 'No TS Found'
                site_dict['hf_ts_energy'] = hf_ts_e

                site_dict['hf_act_energy'] = 'No TS Found'
                try:
                    site_dict['fa+'] = fplus_dict[str(site)]
                    site_dict['fa-'] = fmin_dict[str(site)]
                    site_dict['fa0'] = f0_dict[str(site)]
                    site_dict['electron_density'] = elec_density_dict[str(site)]
                    site_dict['electrostatic_potential'] = esp_dia_dict[str(site)][0]
                    site_dict['diamagnetic_shielding'] = esp_dia_dict[str(site)][1]
                    site_dict['mulliken_charge'] = mulliken_dict[str(site)]
                    site_dict['bond_index'] = bond_indices_dict[str(site)]
                    site_dict['gross_population'] = gross_pop_dict[str(site)]
                except KeyError:
                    site_dict['fa+'] = 'N/A'
                    site_dict['fa-'] = 'N/A'
                    site_dict['fa0'] = 'N/A'
                    site_dict['electron_density'] = 'N/A'
                    site_dict['electrostatic_potential'] = 'N/A'
                    site_dict['diamagnetic_shielding'] = 'N/A'
                    site_dict['mulliken_charge'] = 'N/A'
                    site_dict['bond_index'] = 'N/A'
                    site_dict['gross_population'] = 'N/A'
                radical_dict[str(site)] = site_dict
                radicals_dict = molecule_dict['radicals']
                radicals_dict.update({str(radical): radical_dict})
                molecule_dict.update({'radicals': radicals_dict})
            except UnboundLocalError:
                hf_ts_e = 'No TS Found'
                site_dict['hf_ts_energy'] = hf_ts_e

                site_dict['hf_act_energy'] = 'No TS Found'
                try:
                    site_dict['fa+'] = fplus_dict[str(site)]
                    site_dict['fa-'] = fmin_dict[str(site)]
                    site_dict['fa0'] = f0_dict[str(site)]
                    site_dict['electron_density'] = elec_density_dict[str(site)]
                    site_dict['electrostatic_potential'] = esp_dia_dict[str(site)][0]
                    site_dict['diamagnetic_shielding'] = esp_dia_dict[str(site)][1]
                    site_dict['mulliken_charge'] = mulliken_dict[str(site)]
                    site_dict['bond_index'] = bond_indices_dict[str(site)]
                    site_dict['gross_population'] = gross_pop_dict[str(site)]
                except KeyError:
                    site_dict['fa+'] = 'N/A'
                    site_dict['fa-'] = 'N/A'
                    site_dict['fa0'] = 'N/A'
                    site_dict['electron_density'] = 'N/A'
                    site_dict['electrostatic_potential'] = 'N/A'
                    site_dict['diamagnetic_shielding'] = 'N/A'
                    site_dict['mulliken_charge'] = 'N/A'
                    site_dict['bond_index'] = 'N/A'
                    site_dict['gross_population'] = 'N/A'
                radical_dict[str(site)] = site_dict
                radicals_dict = molecule_dict['radicals']
                radicals_dict.update({str(radical): radical_dict})
                molecule_dict.update({'radicals': radicals_dict})

    hf_site_ea = {} # Could use a single nested dictionary? Will clean this up abit. 
    hf_diff_ea = {}
    hf_arr_ea = {}
    hf_ratio_ea = {}
    am1_site_ea = {}
    am1_diff_ea = {}
    am1_arr_ea = {}
    am1_ratio_ea = {}
    for site in compound_sites:
        am1ea = molecule_dict["radicals"][str(radical)][str(site)]["am1_act_energy"]
        if hpc_calcs:
            hfea = molecule_dict["radicals"][str(radical)][str(site)]["hf_act_energy"]
            print(f'Site is {site} and HF Ea is {hfea}')
            hf_site_ea[site] = hfea
        print(f'site is {site} and AM1 Ea is {am1ea}')

        am1_site_ea[site] = am1ea

    hf_ea_list = []
    am1_ea_list = []

    for key3, value3 in am1_site_ea.items():
        if isinstance(value3, float) is True:
            am1_ea_list.append(value3)
        else:
            continue
    if hpc_calcs:
        for key3, value3 in hf_site_ea.items():
            if isinstance(value3, float) is True:
                hf_ea_list.append(value3)
            else:
                continue
        try:
            mini2 = min(hf_ea_list)
            lowest_ea_site = list(hf_site_ea.keys())[list(hf_site_ea.values()).index(mini2)]
            for key, value in hf_site_ea.items():
                if isinstance(value, str) is True:
                    hf_diff_ea[key] = 'No TS Found'
                else:
                    diff = value - hf_site_ea[lowest_ea_site]
                    if diff > 4:
                        hf_diff_ea[key] = 'Likely Not Observed'
                    else:
                        hf_diff_ea[key] = diff
            for differ_key, differ_value in hf_diff_ea.items():
                if differ_value == 'No TS Found':
                    hf_arr_ea[differ_key] = 'No TS Found'
                elif isinstance(differ_value, float) is True:
                    arr_value = exp(-differ_value / (0.001987 * 348.15))
                    hf_arr_ea[differ_key] = arr_value
                else:
                    hf_arr_ea[differ_key] = 'Likely Not Observed'
            arrh_list = []
            for key2, value2 in hf_arr_ea.items():
                if isinstance(value2, float) is True:
                    arrh_list.append(value2)
                else:
                    continue
            mini = min(arrh_list)
            highest_ea_site = list(hf_arr_ea.keys())[list(hf_arr_ea.values()).index(mini)]
            for arr_key, arr_value in hf_arr_ea.items():
                if arr_value == 'No TS Found':
                    hf_ratio_ea[arr_key] = 'No TS Found'
                elif isinstance(arr_value, float) is True:
                    ratio = round(arr_value * (1 / hf_arr_ea[highest_ea_site]), 1)
                    hf_ratio_ea[arr_key] = ratio
                else:
                    hf_ratio_ea[arr_key] = 'Likely Not Observed'
            # print(compound_sites)
            # print(hf_ratio_ea)
            for site in compound_sites:
                print(hf_ratio_ea[site])
                molecule_dict["radicals"][str(radical)][str(site)]["hf_ratio"] = hf_ratio_ea[site]
        except ValueError:
            for site in compound_sites:
                hf_ratio_ea[site] = "No Ratio"
                molecule_dict["radicals"][str(radical)][str(site)]["hf_ratio"] = hf_ratio_ea[site]

    try:
        mini2 = min(am1_ea_list)
        lowest_ea_site = list(am1_site_ea.keys())[list(am1_site_ea.values()).index(mini2)]
        print(lowest_ea_site)
        for key, value in am1_site_ea.items():
            if isinstance(value, str) is True:
                am1_diff_ea[key] = 'No TS Found'
            else:
                diff = value - am1_site_ea[lowest_ea_site]
                if diff > 4:
                    am1_diff_ea[key] = 'Likely Not Observed'
                else:
                    am1_diff_ea[key] = diff
        for differ_key, differ_value in am1_diff_ea.items():
            if differ_value == 'No TS Found':
                am1_arr_ea[differ_key] = 'No TS Found'
            elif isinstance(differ_value, float) is True:
                arr_value = exp(-differ_value / (0.001987 * 348.15))
                am1_arr_ea[differ_key] = arr_value
            else:
                am1_arr_ea[differ_key] = 'Likely Not Observed'
        arrh_list = []
        for key2, value2 in am1_arr_ea.items():
            if isinstance(value2, float) is True:
                arrh_list.append(value2)
            else:
                continue
        mini = min(arrh_list)
        highest_ea_site = list(am1_arr_ea.keys())[list(am1_arr_ea.values()).index(mini)]
        for arr_key, arr_value in am1_arr_ea.items():
            if arr_value == 'No TS Found':
                am1_ratio_ea[arr_key] = 'No TS Found'
            elif isinstance(arr_value, float) is True:
                ratio = round(arr_value * (1 / am1_arr_ea[highest_ea_site]), 1)
                am1_ratio_ea[arr_key] = ratio
            else:
                am1_ratio_ea[arr_key] = 'Likely Not Observed'
        print(compound_sites)
        print(am1_ratio_ea)
        for site in compound_sites:
            print(am1_ratio_ea[site])
            molecule_dict["radicals"][str(radical)][str(site)]["am1_ratio"] = am1_ratio_ea[site]
    except ValueError:
        for site in compound_sites:
            am1_ratio_ea[site] = "No Ratio"
            molecule_dict["radicals"][str(radical)][str(site)]["am1_ratio"] = am1_ratio_ea[site]
    
    return molecule_dict


# TODO RA: Would be good to have somewere a description of the data dictionary. (What its structure is)