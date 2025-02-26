import os
import chemcoord as cc
import re
from subprocess import DEVNULL, STDOUT, check_call
import glob, shutil
import func_calc.hf_calcs as hf
from app import pycharm
from pathlib import Path

if pycharm == 1: # TODO RA: Is this hard coded? Also is pycharm a variable to identify windows?
    mopac_path = "C:\\Program Files\\MOPAC\\MOPAC2016.exe"
else:
    mopac_path = 'mopac'


def read_mopac_output(in_file, out_file, directory):
    """
    Function that reads the .out file from a MOPAC calculation and decides whether it completed successfully or not and
    writes the TS geometry to a file if yes.
    @param in_file: File to be read.
    @param out_file: File name to write to.
    @param directory: Directory where files can be found.
    @return: Code of either 1 (completed successfully) or 2 (Failed calculation)
    """ # TODO RA: docstrings should really contain types. 
    # Read am1 output line by line but reverse order and put into new temporary file out_rev.txt.
    # TODO RA: This feels like making a temp file just to have do some file analysis is abit inefficient... 

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
        coordinates.append(coordinate1) # TODO RA: .append is bad... incredibly slow. try: Nat = len(trimmed); coordinates = [None]*Nat; coordinates[i] = coordinate1 (I know this doesnt actually work but something along these lines would be better)
    system_size = len(trimmed)
    # Write AM1 optimised coordinates to an xyz file
    with open(out_file, 'w+') as xyz:
        xyz.write(str(system_size) + '\n\n')
        for i in coordinates:
            xyz.write(str(i) + '\n')
    os.remove(Path(f'{directory}/out_rev.txt')) # TODO RA: I dont think this being a file is useful... just read it into a list and then delete the list
    return 1


def mopac_distance_checks(geometry, carbon_radical, site):
    """
    Function that checks whether the transition state geometry has a C-C bond length within a typical range.
    @param geometry: Cartesian coordinates file of transition state.
    @param carbon_radical: Which Functional group is being investigated.
    @param site: Which site within the compound is being investigated.
    @return: Boolean that decides whether the distance is within the typical range.
    """ # TODO RA: again should have types
    try:
        with open(f'{geometry}', 'r+')as f:
            lin = f.readlines()
            coords = lin[2:]
            cart = []
            for atom in coords:
                strip = atom.strip('\n')
                line_elements = re.split(' +', str(strip))
                line = []
                # print(line_elements)
                for element in line_elements:
                    line.append(element.split('\t'))
                flat_list = [item for sublist in line for item in sublist]
                cart.append(flat_list) # TODO RA: Again... avoid .append
            if carbon_radical == 'cf3':
                carbon = cart[-4]
            elif carbon_radical == 'cf2h':
                carbon = cart[-4]
            elif carbon_radical == 'ipr':
                carbon = cart[-10]
            site_of_interest = cart[site-1]

            [a, x1, y1, z1] = carbon
            [b, x2, y2, z2] = site_of_interest
            dist = float((((float(x2) - float(x1)) ** 2) + ((float(y2) - float(y1)) ** 2) + (
                        (float(z2) - float(z1)) ** 2)) ** (1 / 2))
            # print(dist)
            if 1.8 <= dist <= 2.2:
                return True
            else:
                return False
    except FileNotFoundError:
        return False


def read_mopac_freq(input_file, directory):
    """return codes: 1 = true transition state
                     2 = close to transition state with a second small imaginary frequency
                     3 = 3 or more negative frequencies or large second imaginary frequency, needs further intervention
                     4 = no imaginary frequency in correct range, needs other intervention""" # TODO RA: This should really have docstring formatting. 
    # Read am1 output line by line but reverse order and put into new temporary file out_rev.txt.
    with open(input_file) as rf, open(Path(f'{directory}/out_rev.txt'), 'w') as wf: # TODO RA: Is this just copying the file??? if not os.system(f'cp {input_file} {Path(f'{directory}/out_rev.txt')}') would be better??
        for line in rf.readlines():
            wf.write(line)

    # Create empty list to put cartesian coordinates from output file into
    a = []
    # Scan the  file for the keyword MASS-WEIGHTED COORDINATE ANALYSIS (NORMAL COORDINATES) and collect the following
    # lines until Root No.       9        10        11        12        13        14        15        16 is reached and
    # add them to the list a.
    with open(Path(f'{directory}/out_rev.txt')) as file:
        for line in file:
            lines = line.strip('\n')
            if 'MASS-WEIGHTED COORDINATE ANALYSIS (NORMAL COORDINATES)' in lines:
                # collect block-related lines
                while True:
                    try:
                        lines = next(file)
                    except StopIteration:
                        # there is no lines left
                        break
                    if 'Root No.       9        10        11        12        13        14        15        16' in lines:
                        # we've reached the end of block
                        break
                    a.append(lines)
                # stop iterating over file
                break
    # Remove the lines CARTESIAN COORDINATES and EIGENVALUES and the blank lines before and after the coordinates of
    # interest.
    os.remove(Path(f'{directory}/out_rev.txt')) # TODO RA: Again... why is this a file? as apposed to just an object. Temp files should be for large data sets that may cause memory issues...
    try:
        b = re.split(' +', a[6])[1:]
        if - 800 <= float(b[0]) <= - 300:
            if 0.00 <= float(b[1]):
                return 1
            elif - 100 <= float(b[1]) <= - 0.01:
                if 0.00 <= float(b[2]):
                    return 2
                elif - 100 <= float(b[2]) <= - 0.01:
                    return 3
            elif - 600 <= float(b[1]) <= - 101.1:
                return 3
        else:
            return 4
    except IndexError:
        return 4


def mopac_freq_check(inpu, positively_charged, radical, site, directory, new_mopac_2016=False, multiple=False):
    """
    Function that Decides whether further TS searches are needed by running a frequency calculation on the resultant
    geometry.
    @param inpu: File name to work with.
    @param positively_charged: Whether the compound is a charged species or not.
    @param radical: Name of the functional group being investigated.
    @param site: Site that is being investigated.
    @param directory: Folder where the files are found.
    @param new_mopac_2016: Whether the new MOPAC 2016 is being used or not.
    @param multiple: Whether this is a calculation for multiple compounds at once or not.
    @return: Boolean 
    """ # TODO RA: types? Also what does the boolean mean?
    check_call(['obabel', '-ixyz', str(inpu), '-omopin', '-O', Path(f'{directory}/tmp.dat')],
               stdout=DEVNULL, stderr=STDOUT)
    # os.system(f'obabel -ixyz {inpu} -omopin -O tmp.dat > /dev/null 2>&1')
    with open(Path(f'{directory}/tmp.dat'), 'r+') as file:
        lines = file.readlines()
        if positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 FORCE LET MMOK GEO-OK UHF\n'
            else:
                lines[0] = 'AM1 FORCE LET MMOK GEO-OK DISP UHF\n'
        else:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 FORCE LET MMOK GEO-OK UHF CHARGE=+1\n'
            else:
                lines[0] = 'AM1 FORCE LET MMOK GEO-OK DISP UHF CHARGE=+1\n'
    os.remove(Path(f'{directory}/tmp.dat'))
    with open(Path(f'{directory}/tmp.dat'), 'w+') as tmp:
        for string in lines:
            tmp.write(str(string))
    check_call([mopac_path, Path(f'{directory}/tmp')], stdout=DEVNULL, stderr=STDOUT)
    actual_ts = read_mopac_freq(Path(f'{directory}/tmp.out'), directory)
    if actual_ts == 1:
        print('Transition state')
        # os.system('ls')
        hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=False, multiple=multiple)
        with open(Path(f'{directory}/am1_converged.txt'), 'w+') as am1:
            am1.write('semi-empirircal converged')
        files = Path(f'{directory}/').glob('tmp.*')
        for file in files:
            os.remove(file)
        return True
    elif actual_ts == 2:
        print('Close to Transition State')
        hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=False, multiple=multiple)
        with open(Path(f'{directory}/am1_converged.txt'), 'w+') as am1:
            am1.write('semi-empirircal converged')
        files = Path(f'{directory}/').glob('tmp.*')
        for file in files:
            os.remove(file)
        return True
    elif actual_ts == 3:
        # print('2nd step')
        files = Path(f'{directory}/').glob('tmp.*')
        for file in files:
            os.remove(file)
        if positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF\n'
        else:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF CHARGE=+1\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF CHARGE=+1\n'
        with open(Path(f'{directory}/tmp.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))
        check_call([mopac_path, Path(f'{directory}/tmp')], stdout=DEVNULL, stderr=STDOUT)
        if read_mopac_output(Path(f'{directory}/tmp.out'), Path(f'{directory}/ts2.xyz'), Path(f'{directory}')) == 1:

            if mopac_distance_checks(Path(f'{directory}/ts2.xyz'), radical, site) is True:
                check_call(['obabel', '-ixyz', Path(f'{directory}/ts2.xyz'), '-omopin', '-O',
                            Path(f'{directory}/tmp2.dat')], stdout=DEVNULL, stderr=STDOUT)
                with open(Path(f'{directory}/tmp2.dat'), 'r+') as file:
                    lines = file.readlines()
                    if positively_charged is False:
                        if new_mopac_2016 is False:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK UHF\n'
                        else:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK DISP UHF\n'
                    else:
                        if new_mopac_2016 is False:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK UHF CHARGE=+1\n'
                        else:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK DISP UHF CHARGE=+1\n'
                os.remove(Path(f'{directory}/tmp2.dat'))
                with open(Path(f'{directory}/tmp2.dat'), 'w+') as tmp:
                    for string in lines:
                        tmp.write(str(string))
                check_call([mopac_path, Path(f'{directory}/tmp2')], stdout=DEVNULL, stderr=STDOUT)
                actual_ts = read_mopac_freq(Path(f'{directory}/tmp2.out'), Path(f'{directory}'))
                if actual_ts == 1:
                    print('Now TS')
                    shutil.move(Path(f'{directory}/ts2.xyz'), Path(f'{directory}/ts.xyz'))
                    hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=False, multiple=multiple)
                    with open(Path(f'{directory}/am1_converged.txt'), 'w+') as am1:
                        am1.write('semi-empirircal converged')
                    files = Path(f'{directory}/').glob('tmp.*')
                    for file in files:
                        os.remove(file)
                    # os.system('mv ts2.xyz ts.xyz')
                    # os.system('rm ./tmp.*')
                    return True
                elif actual_ts == 2:
                    print('Close to Transition State')
                    shutil.move(Path(f'{directory}/ts2.xyz'), Path(f'{directory}/ts.xyz'))
                    hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=False, multiple=multiple)
                    with open(Path(f'{directory}/am1_converged.txt'), 'w+') as am1:
                        am1.write('semi-empirircal converged')
                    files = Path(f'{directory}/').glob('tmp.*')
                    for file in files:
                        os.remove(file)
                    # os.system('mv ts2.xyz ts.xyz')
                    # os.system('rm ./tmp.*')
                    return True
                elif actual_ts == 3 or actual_ts == 4:
                    print('Fail')
                    hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=True, multiple=multiple)
                    files = Path(f'{directory}/').glob('tmp.*')
                    for file in files:
                        os.remove(file)
                    # os.system('rm ./tmp.*')
                    return False
            else:
                print('Fail')
                hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=True, multiple=multiple)
                files = Path(f'{directory}/').glob('tmp.*')
                for file in files:
                    os.remove(file)
                # os.system('rm ./tmp.*')
                return False
        else:
            print('Fail')
            hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=True, multiple=multiple)
            files = Path(f'{directory}/').glob('tmp.*')
            for file in files:
                os.remove(file)
            # os.system('rm ./tmp.*')
            return False

    elif actual_ts == 4:
        files = Path(f'{directory}/').glob('tmp.*')
        for file in files:
            os.remove(file)
        if positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF\n'
        else:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF CHARGE=+1\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF CHARGE=+1\n'
        with open(Path(f'{directory}/tmp.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))
        check_call([mopac_path, Path(f'{directory}/tmp')], stdout=DEVNULL, stderr=STDOUT)
        if read_mopac_output(Path(f'{directory}/tmp.out'), Path(f'{directory}/ts2.xyz'), Path(f'{directory}')) == 1:
            if mopac_distance_checks(Path(f'{directory}/ts2.xyz'), radical, site) is True:
                check_call(['obabel', '-ixyz', Path(f'{directory}/ts2.xyz'), '-omopin', '-O',
                            Path(f'{directory}/tmp2.dat')], stdout=DEVNULL, stderr=STDOUT)
                # os.system(f'obabel -ixyz ts2.xyz -omopin -O tmp2.dat > /dev/null 2>&1')
                with open(Path(f'{directory}/tmp2.dat'), 'r+') as file:
                    lines = file.readlines()
                    if positively_charged is False:
                        if new_mopac_2016 is False:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK UHF\n'
                        else:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK DISP UHF\n'
                    else:
                        if new_mopac_2016 is False:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK UHF CHARGE=+1\n'
                        else:
                            lines[0] = 'AM1 FORCE LET MMOK GEO-OK DISP UHF CHARGE=+1\n'
                os.remove(Path(f'{directory}/tmp2.dat'))
                with open(Path(f'{directory}/tmp2.dat'), 'w+') as tmp:
                    for string in lines:
                        tmp.write(str(string))
                check_call([mopac_path, Path(f'{directory}/tmp2')], stdout=DEVNULL, stderr=STDOUT)
                actual_ts = read_mopac_freq(Path(f'{directory}/tmp2.out'), Path(f'{directory}'))
                path2 = os.getcwd()
                if actual_ts == 1:
                    print('Now TS')
                    shutil.move(Path(f'{directory}/ts2.xyz'), Path(f'{directory}/ts.xyz'))
                    hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=False, multiple=multiple)
                    with open(Path(f'{directory}/am1_converged.txt'), 'w+') as am1:
                        am1.write('semi-empirircal converged')
                    files = Path(f'{directory}/').glob('tmp.*')
                    for file in files:
                        os.remove(file)
                    # os.system('mv ts2.xyz ts.xyz')
                    # os.system('rm ./tmp.*')
                    return True
                elif actual_ts == 2:
                    print('Close to Transition State')
                    shutil.move(Path(f'{directory}/ts2.xyz'), Path(f'{directory}/ts.xyz'))
                    hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=False, multiple=multiple)
                    with open(Path(f'{directory}/am1_converged.txt'), 'w+') as am1:
                        am1.write('semi-empirircal converged')
                    files = Path(f'{directory}/').glob('tmp.*')
                    for file in files:
                        os.remove(file)
                    # os.system('mv ts2.xyz ts.xyz')
                    # os.system('rm ./tmp.*')
                    return True
                elif actual_ts == 3 or actual_ts == 4:
                    print('Fail')
                    hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=True, multiple=multiple)
                    files = Path(f'{directory}/').glob('tmp.*')
                    for file in files:
                        os.remove(file)
                    # os.system('rm ./tmp.*')
                    return False
            else:
                print('Fail')
                hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=True, multiple=multiple)
                files = Path(f'{directory}/').glob('tmp.*')
                for file in files:
                    os.remove(file)
                # os.system('rm ./tmp.*')
                return False
        else:
            print('Fail')
            hf.generate_hf_geometry(Path(f'{directory}'), radical, failed=True, multiple=multiple)
            files = Path(f'{directory}/').glob('tmp.*')
            for file in files:
                os.remove(file)
            # os.system('rm ./tmp.*')
            return False


def check_for_clashes(geometry, query_atom, directory, clash_threshold=1): # TODO RA: directory is not used. 
    """
    Function that checks the built transition state geometry from the template to see whether there are any clashes
    between the compound of interest and the radical being added to the system. It checks the distance between each
     radical atom and each atom in the compound structure is no closer than the threshold distance.
    @param geometry: The input file for the function.
    @param query_atom: The atom number in system that distance checks need ot performed on.
    @param directory: The firectory the calculations are being executed in
    @param clash_threshold: The minimum distance between the atom of interest and other atoms in the system.
    @return: """ # TODO RA: types?
    clashes = False

    with open(f'{geometry}', 'r+') as inp:
        lines = inp.readlines()
        cartesian = lines[2:]
        cart = []
        for atom in cartesian:
            strip = atom.strip('\n')
            line_elements = re.split(' +', str(strip))
            line = []
            # print(line_elements)
            for element in line_elements:
                line.append(element.split('\t'))
            flat_list = [item for sublist in line for item in sublist]
            cart.append(flat_list) # TODO RA: .append... 
    distances = []
    atom_of_concern = cart[query_atom-1]
    for atom_no in cart:
        if atom_no == atom_of_concern:
            pass
        else:
            [a, x1, y1, z1] = atom_no
            [b, x2, y2, z2] = atom_of_concern
            dist = float((((float(x2) - float(x1)) ** 2) + ((float(y2) - float(y1)) ** 2) + ((float(z2) - float(z1)) ** 2)) ** (1 / 2))
            distances.append(dist)
    for length in distances:
        if length < clash_threshold:
            clashes = True
        else:
            pass

    return clashes


def cf3(input_file, xyz, site, directory, positively_charged=False, constrained=False, precise=False, new_mopac_2016=False):
    """
    Function that sets up calculation for adding the CF3 group to the compound of interest in each site from a template
    AM1 transition state geometry.
    @param input_file: MOPAC .dat file containing the optimised geometry of the isolated compound of interest with no
     radical added.
    @param xyz: Cartesian coordinates of the optimised structure.
    @param site: Site within the compound to add the radical to, building the template.
    @param directory: Directory we are reading/writing files.
    @param positively_charged: Boolean on whether the compound is positively charged or not.
    @param constrained: Boolean on whether the calculation should be constrained. Failed initial TS searches will be set
     off again with an initial constrained optimisation and subsequent relaxed TS search.
    @param precise: Whether the precise keyword should be added to the MOPAC calculation.
    @param new_mopac_2016: Whether we are running on the new mopac_2016 version or not (alters how the output
    information is formatted) 
    """# TODO RA: types?
    with open(input_file) as file:

        lines = file.readlines()

        if precise is False and positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK DISP CYCLES=10000 UHF\n'
        elif precise is False and positively_charged is True:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK CHARGE=+1 UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK CHARGE=+1 DISP UHF\n'
        elif precise is True and positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP CYCLES=10000 UHF\n'
        elif precise is True and positively_charged is True:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE CHARGE=+1 UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE CHARGE=+1 DISP UHF\n'

        if not constrained:
            const = 1
        else:
            const = 0
            ts = lines[0].split()
            ts.remove('TS')
            ts[-1] = 'UHF\n'
            t = " "
            t = t.join(ts)
            lines[0] = t
        # Generate a connection table relating each atom's position to the previous in the carbon chain. A carbon's
        # position is not related to the position of the hydrogen.
        compound = cc.Cartesian.read_xyz(str(xyz), start_index=1)
        with open(str(xyz)) as coord:
            coords = coord.readlines()
            atom_count = int(coords[0])
        connection_table = compound.get_bonds()
        site_connection = connection_table[site]
        for bonded in site_connection:
            if len(connection_table[bonded]) != 1 and bonded != site:
                angle_connection = connection_table[bonded]
                angle_atom = bonded

            else:
                pass
        dihedral_atoms = []
        for dihed in angle_connection:
            if len(connection_table[dihed]) != 1 and dihed != angle_atom and dihed != site:
                # dihedral_connection = connection_table[dihed]
                dihedral_atoms.append(dihed) 
            else:
                pass
        dihedral_atom = min(dihedral_atoms)
        # print(f'site is: {site}, angle atom is: {angle_atom}, dihedral atom is: {dihedral_atom}')
        c_dihed = 102.600000
        f1_dihed = -179.300000
        f2_dihed = 63.800000
        f3_dihed = -58.000000

        lines.append(f'C    2.000000  {const}  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    1.350000  1  107.800000  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    1.345000  1  114.600000  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]
        while any(check_atoms) is True:
            # print('here')
            rotation_spacing = 360/72
            c_dihed += rotation_spacing
            f1_dihed += rotation_spacing
            f2_dihed += rotation_spacing
            f3_dihed += rotation_spacing

            lines[-4] = f'C    2.000000  {const}  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[-3] = f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-2] = f'F    1.350000  1  107.800000  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-1] = f'F    1.345000  1  114.600000  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))
            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]

    os.remove(Path(f'{directory}/temp_ts.dat'))
    os.remove(Path(f'{directory}/temp_ts.xyz'))
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)


def cf2h(input_file, xyz, site, directory, positively_charged=False, constrained=False, precise=False, new_mopac_2016=False):
    """
    Function that sets up calculation for adding the CF3 group to the compound of interest in each site from a template
    AM1 transition state geometry.
    @param input_file: MOPAC .dat file containing the optimised geometry of the isolated compound of interest with no
     radical added.
    @param xyz: Cartesian coordinates of the optimised structure.
    @param site: Site within the compound to add the radical to, building the template.
    @param directory: Directory we are reading/writing files.
    @param positively_charged: Boolean on whether the compound is positively charged or not.
    @param constrained: Boolean on whether the calculation should be constrained. Failed initial TS searches will be set
     off again with an initial constrained optimisation and subsequent relaxed TS search.
    @param precise: Whether the precise keyword should be added to the MOPAC calculation.
    @param new_mopac_2016: Whether we are running on the new mopac_2016 version or not (alters how the output
    information is formatted)
    """ # TODO RA: types?
    with open(input_file) as file:

        lines = file.readlines()

        if precise is False and positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK DISP UHF\n'
        elif precise is False and positively_charged is True:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK CHARGE=+1 UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK CHARGE=+1 DISP UHF\n'
        elif precise is True and positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF\n'
        elif precise is True and positively_charged is True:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE CHARGE=+1 UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE CHARGE=+1 DISP UHF\n'
        # TODO RA: This could be a fstring ^^^

        if not constrained:
            const = 1
        else:
            const = 0
            ts = lines[0].split()
            ts.remove('TS')
            ts[-1] = 'UHF\n'
            t = " "
            t = t.join(ts)
            lines[0] = t

        compound = cc.Cartesian.read_xyz(str(xyz), start_index=1)
        with open(str(xyz)) as coord:
            coords = coord.readlines()
            atom_count = int(coords[0])
        connection_table = compound.get_bonds()
        site_connection = connection_table[site]
        for bonded in site_connection:
            if len(connection_table[bonded]) != 1 and bonded != site:
                angle_connection = connection_table[bonded]
                angle_atom = bonded

            else:
                pass
        dihedral_atoms = []
        for dihed in angle_connection:
            if len(connection_table[dihed]) != 1 and dihed != angle_atom and dihed != site:
                # dihedral_connection = connection_table[dihed]
                dihedral_atoms.append(dihed)
            else:
                pass
        dihedral_atom = min(dihedral_atoms)
        # print(f'site is: {site}, angle atom is: {angle_atom}, dihedral atom is: {dihedral_atom}')
        c_dihed = 102.600000
        f1_dihed = 57.2481322
        h_dihed = -66.4380492
        f2_dihed = 170.4046283

        lines.append(f'C    2.000000  {const}  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'H    1.350000  1  107.800000  1  {h_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    1.345000  1  114.600000  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]
        while any(check_atoms) is True:
            # print('here')
            rotation_spacing = 360/72
            c_dihed += rotation_spacing
            f1_dihed += rotation_spacing
            h_dihed += rotation_spacing
            f2_dihed += rotation_spacing

            lines[-4] = f'C    2.000000  {const}  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[-3] = f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-2] = f'H    1.350000  1  107.800000  1  {h_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-1] = f'F    1.345000  1  114.600000  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]

    os.remove(Path(f'{directory}/temp_ts.dat'))
    os.remove(Path(f'{directory}/temp_ts.xyz'))
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)


def ipr(input_file, xyz, site, directory, positively_charged=False, constrained=False, precise=False, new_mopac_2016=False):
    """
    Function that sets up calculation for adding the CF3 group to the compound of interest in each site from a template
    AM1 transition state geometry.
    @param input_file: MOPAC .dat file containing the optimised geometry of the isolated compound of interest with no
     radical added.
    @param xyz: Cartesian coordinates of the optimised structure.
    @param site: Site within the compound to add the radical to, building the template.
    @param directory: Directory we are reading/writing files.
    @param positively_charged: Boolean on whether the compound is positively charged or not.
    @param constrained: Boolean on whether the calculation should be constrained. Failed initial TS searches will be set
     off again with an initial constrained optimisation and subsequent relaxed TS search.
    @param precise: Whether the precise keyword should be added to the MOPAC calculation.
    @param new_mopac_2016: Whether we are running on the new mopac_2016 version or not (alters how the output
    information is formatted)
    """
    with open(input_file) as file:

        lines = file.readlines()

        if precise is False and positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK DISP UHF\n'
        elif precise is False and positively_charged is True:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK CHARGE=+1 UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK CHARGE=+1 DISP UHF\n'
        elif precise is True and positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE DISP UHF\n'
        elif precise is True and positively_charged is True:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE CHARGE=+1 UHF\n'
            else:
                lines[0] = 'AM1 TS LET MMOK GEO-OK PRECISE CHARGE=+1 DISP UHF\n'

        if not constrained:
            const = 1
        else:
            const = 0
            ts = lines[0].split()
            ts.remove('TS')
            ts[-1] = 'UHF\n'
            t = " "
            t = t.join(ts)
            lines[0] = t

        compound = cc.Cartesian.read_xyz(str(xyz), start_index=1)
        with open(str(xyz)) as coord:
            coords = coord.readlines()
            atom_count = int(coords[0])
        connection_table = compound.get_bonds()
        site_connection = connection_table[site]
        for bonded in site_connection:
            if len(connection_table[bonded]) != 1 and bonded != site:
                angle_connection = connection_table[bonded]
                angle_atom = bonded

            else:
                pass
        dihedral_atoms = []
        for dihed in angle_connection:
            if len(connection_table[dihed]) != 1 and dihed != angle_atom and dihed != site:
                # dihedral_connection = connection_table[dihed]
                dihedral_atoms.append(dihed)
            else:
                pass
        dihedral_atom = min(dihedral_atoms) # TODO RA: Should these be hard coded? If so consider having a dictionary of standard values and store at the top of the file. Makes editing much easier. 
        # print(f'site is: {site}, angle atom is: {angle_atom}, dihedral atom is: {dihedral_atom}')
        c1_dihed = 254.264918
        c2_dihed = 301.336871
        c3_dihed = 179.509647
        h1_dihed = 60.272731
        h2_dihed = 291.898120
        h3_dihed = 172.085792
        h4_dihed = 52.709817
        h5_dihed = 69.148928
        h6_dihed = 308.139965
        h7_dihed = 188.805597
        # TODO RA: Again is all this meant to be hard coded?
        lines.append(f'C    2.021139  {const}  100.078262  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(f'C    1.482068  1  104.004890  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'C    1.481423  1  104.119367  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'H    1.103400  1   97.795054  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'H    1.117064  1  111.484379  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n')
        lines.append(f'H    1.120711  1  109.826665  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
        lines.append(f'H    1.118015  1  110.943566  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
        lines.append(f'H    1.117322  1  111.323645  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
        lines.append(f'H    1.118059  1  111.047346  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
        lines.append(f'H    1.120507  1  109.884160  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-9, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-8, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-7, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-6, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-5, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-4, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]
        while any(check_atoms) is True:
            # print('here')
            rotation_spacing = 360/72
            c1_dihed += rotation_spacing
            c2_dihed += rotation_spacing
            c3_dihed += rotation_spacing
            h1_dihed += rotation_spacing
            h2_dihed += rotation_spacing
            h3_dihed += rotation_spacing
            h4_dihed += rotation_spacing
            h5_dihed += rotation_spacing
            h6_dihed += rotation_spacing
            h7_dihed += rotation_spacing

            lines[-10] = f'C    2.021139  {const}  100.078262  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[-9] = f'C    1.482068  1  104.004890  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-8] = f'C    1.481423  1  104.119367  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-7] = f'H    1.103400  1   97.795054  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-6] = f'H    1.117064  1  111.484379  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n'
            lines[-5] = f'H    1.120711  1  109.826665  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
            lines[-4] = f'H    1.118015  1  110.943566  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
            lines[-3] = f'H    1.117322  1  111.323645  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
            lines[-2] = f'H    1.118059  1  111.047346  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
            lines[-1] = f'H    1.120507  1  109.884160  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-9, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-8, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-7, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-6, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-5, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-4, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]

    os.remove(Path(f'{directory}/temp_ts.dat'))
    os.remove(Path(f'{directory}/temp_ts.xyz'))
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)


def cf3_reactant(input_file, xyz, site, directory, positively_charged=False, reagent_optimised=True, new_mopac_2016=False, multiple=False):

    with open(input_file) as file:

        lines = file.readlines()

        if positively_charged is False:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 LET MMOK GEO-OK PRECISE\n'
            else:
                lines[0] = 'AM1 LET MMOK GEO-OK DISP PRECISE\n'
        else:
            if new_mopac_2016 is False:
                lines[0] = 'AM1 LET MMOK GEO-OK PRECISE CHARGE=+1\n'
            else:
                lines[0] = 'AM1 LET MMOK GEO-OK PRECISE DISP CHARGE=+1\n'
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))
    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
    with open(f'{input_file}', 'w+') as file:
        compound = cc.Cartesian.read_xyz(str(xyz), start_index=1)
        with open(str(xyz)) as coord:
            coords = coord.readlines()
            atom_count = int(coords[0])
        connection_table = compound.get_bonds()
        site_connection = connection_table[site]
        for bonded in site_connection:
            if len(connection_table[bonded]) != 1 and bonded != site:
                angle_connection = connection_table[bonded]
                angle_atom = bonded

            else:
                pass
        dihedral_atoms = []
        for dihed in angle_connection:
            if len(connection_table[dihed]) != 1 and dihed != angle_atom and dihed != site:
                # dihedral_connection = connection_table[dihed]
                dihedral_atoms.append(dihed)
            else:
                pass
        dihedral_atom = min(dihedral_atoms)
        # print(f'site is: {site}, angle atom is: {angle_atom}, dihedral atom is: {dihedral_atom}')
        c_dihed = 102.600000
        f1_dihed = -179.300000
        f2_dihed = 63.800000
        f3_dihed = -58.000000
        fukui = lines.copy()
        fukui_dir = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-1])}{str(Path("/fukui"))}'
        lines.append(f'C    12.000000  1  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    1.350000  1  107.800000  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    1.345000  1  114.600000  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]
        while any(check_atoms) is True:
            # print('here')
            rotation_spacing = 360/72
            c_dihed += rotation_spacing
            f1_dihed += rotation_spacing
            f2_dihed += rotation_spacing
            f3_dihed += rotation_spacing

            lines[-4] = f'C    12.000000  1  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[-3] = f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-2] = f'F    1.350000  1  107.800000  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-1] = f'F    1.345000  1  114.600000  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]

        os.remove(Path(f'{directory}/temp_ts.dat'))
        os.remove(Path(f'{directory}/temp_ts.xyz'))
        lines[-4] = f'C    50.000000  1  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
    with open(Path(f'{directory}/HF_reagent.dat'), 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))
    with open(Path(f'{fukui_dir}/hf.dat'), 'w+') as f:
        for st in fukui:
            f.write(str(st))

    check_call(['obabel', '-imopin', Path(f'{directory}/HF_reagent.dat'), '-oxyz', '-O',
                Path(f'{directory}/HF_reagent.xyz')], stdout=DEVNULL, stderr=STDOUT)
    check_call(['obabel', '-imopin', Path(f'{directory}/HF_reagent.dat'), '-osdf', '-O',
                Path(f'{directory}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
    check_call(['obabel', '-imopin', Path(f'{fukui_dir}/hf.dat'), '-osdf', '-O',
                Path(f'{fukui_dir}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
    hf_calcs_file = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-2])}{str(Path("/hf_jobs.txt"))}'
    if os.path.isfile(hf_calcs_file) is True:
        with open(hf_calcs_file, 'r') as fil:
            lins = fil.readlines()
        for li in lins:
            li = str(Path(li)).split(str(Path('/')))
            # removes the previous hf_jobs.txt file if other compounds calculated in the directory before. Also works
            # for multiple compounds calculated at once where the folder structure is slightly different.
            if not multiple:
                print('in single')
                if str(directory).split(str(Path('/')))[-2] != li[-3]:
                    print('REMOVING FILE')
                    os.remove(hf_calcs_file)
                    break
            else:
                print(directory, li)
                print(str(directory).split(str(Path('/')))[-3], li[-4])
                if str(directory).split(str(Path('/')))[-3] != li[-4]:
                    print('REMOVING FILE')
                    os.remove(hf_calcs_file)
                    break
        with open(hf_calcs_file, 'a') as fi:
            fi.write(str(Path(f'{directory}/hf.sdf\n')))
            fi.write(str(Path(f'{fukui_dir}/hf.sdf\n')))
    else:
        with open(hf_calcs_file, 'w') as fil:
            fil.write(str(Path(f'{directory}/hf.sdf\n')))
            fil.write(str(Path(f'{fukui_dir}/hf.sdf\n')))
    # os.system('obabel -imopin HF_reagent.dat -oxyz -O HF_reagent.xyz > /dev/null 2>&1')


def cf2h_reactant(input_file, xyz, site, directory, positively_charged=False, reagent_optimised=True, new_mopac_2016=False, multiple=False):

    if not reagent_optimised:
        with open(input_file) as file:

            lines = file.readlines()

            if positively_charged is False: # TODO RA: You should make a function that does this, rather than re-making the same code over and over.
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK DISP PRECISE\n'
            else:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE CHARGE=+1\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE DISP CHARGE=+1\n'
        os.remove(f'{input_file}')
        with open(f'{input_file}', 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))
        check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
    else:
        with open(input_file) as file:

            lines = file.readlines()

            if positively_charged is False:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK DISP PRECISE\n'
            else:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE CHARGE=+1\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE DISP CHARGE=+1\n'

            compound = cc.Cartesian.read_xyz(str(xyz), start_index=1)
            with open(str(xyz)) as coord:
                coords = coord.readlines()
                atom_count = int(coords[0])
            connection_table = compound.get_bonds()
            site_connection = connection_table[site]
            for bonded in site_connection:
                if len(connection_table[bonded]) != 1 and bonded != site:
                    angle_connection = connection_table[bonded]
                    angle_atom = bonded

                else:
                    pass
            dihedral_atoms = []
            for dihed in angle_connection:
                if len(connection_table[dihed]) != 1 and dihed != angle_atom and dihed != site:
                    # dihedral_connection = connection_table[dihed]
                    dihedral_atoms.append(dihed)
                else:
                    pass
            dihedral_atom = min(dihedral_atoms)
            # print(f'site is: {site}, angle atom is: {angle_atom}, dihedral atom is: {dihedral_atom}') 
            # TODO RA: Should these be hard coded? If so consider having a dictionary of standard values and store at the top of the file. Makes editing much easier.
            c_dihed = 102.600000
            f1_dihed = 57.2481322
            h_dihed = -66.4380492
            f2_dihed = 170.4046283
            fukui = lines.copy()
            fukui_dir = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-1])}{str(Path("/fukui"))}'
            lines.append(f'C    12.000000  1  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
            lines.append(f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'H    1.350000  1  107.800000  1  {h_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'F    1.345000  1  114.600000  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]

            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]
            while any(check_atoms) is True:
                # print('here')
                rotation_spacing = 360/72
                c_dihed += rotation_spacing
                f1_dihed += rotation_spacing
                h_dihed += rotation_spacing
                f2_dihed += rotation_spacing

                lines[-4] = f'C    12.000000  1  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
                lines[-3] = f'F    1.351000  1  107.100000  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-2] = f'H    1.350000  1  107.800000  1  {h_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-1] = f'F    1.345000  1  114.600000  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

                os.remove(Path(f'{directory}/temp_ts.dat'))
                os.remove(Path(f'{directory}/temp_ts.xyz'))
                with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                    for string in lines:
                        tmp.write(str(string))

                check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz',
                            '-O', Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
                with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                    lin = inp.readlines()
                    cartesian = lin[2:]
                check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]

        os.remove(Path(f'{directory}/temp_ts.dat'))
        os.remove(Path(f'{directory}/temp_ts.xyz'))
        os.remove(f'{input_file}')
        with open(f'{input_file}', 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))

        check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
        lines[-4] = f'C    50.000000  1  101.000000  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
        with open(Path(f'{directory}/HF_reagent.dat'), 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))
        with open(Path(f'{fukui_dir}/hf.dat'), 'w+') as f:
            for st in fukui:
                f.write(str(st))
        check_call(['obabel', '-imopin', Path(f'{directory}/HF_reagent.dat'), '-oxyz', '-O',
                    Path(f'{directory}/HF_reagent.xyz')], stdout=DEVNULL, stderr=STDOUT)
        check_call(['obabel', '-imopin', Path(f'{directory}/HF_reagent.dat'), '-osdf', '-O',
                    Path(f'{directory}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
        check_call(['obabel', '-imopin', Path(f'{fukui_dir}/hf.dat'), '-osdf', '-O',
                    Path(f'{fukui_dir}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
        hf_calcs_file = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-2])}{str(Path("/hf_jobs.txt"))}'
        if os.path.isfile(hf_calcs_file) is True:
            with open(hf_calcs_file, 'r') as fil:
                lins = fil.readlines()
            for li in lins:
                li = str(Path(li)).split(str(Path('/')))
                if not multiple:
                    if str(directory).split(str(Path('/')))[-2] != li[-3]:
                        print('REMOVING FILE')
                        os.remove(hf_calcs_file)
                        break
                else:
                    if str(directory).split(str(Path('/')))[-3] != li[-4]:
                        print('REMOVING FILE')
                        os.remove(hf_calcs_file)
                        break
            with open(hf_calcs_file, 'a') as fi:
                fi.write(str(Path(f'{directory}/hf.sdf\n')))
                fi.write(str(Path(f'{fukui_dir}/hf.sdf\n')))
        else:
            with open(hf_calcs_file, 'w') as fil:
                fil.write(str(Path(f'{directory}/hf.sdf\n')))
                fil.write(str(Path(f'{fukui_dir}/hf.sdf\n')))
        # os.system('obabel -imopin HF_reagent.dat -oxyz -O HF_reagent.xyz > /dev/null 2>&1')


def ipr_reactant(input_file, xyz, site, directory, positively_charged=False, reagent_optimised=True, new_mopac_2016=False, multiple=False):
    # TODO RA: docstring? 
    if not reagent_optimised:
        with open(input_file) as file:

            lines = file.readlines()

            if positively_charged is False:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK DISP PRECISE\n'
            else:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE CHARGE=+1\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE DISP CHARGE=+1\n'
        os.remove(f'{input_file}')
        with open(f'{input_file}', 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))
        check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
    else:
        with open(input_file) as file:

            lines = file.readlines()

            if positively_charged is False:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK DISP PRECISE\n'
            else:
                if new_mopac_2016 is False:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE CHARGE=+1\n'
                else:
                    lines[0] = 'AM1 LET MMOK GEO-OK PRECISE DISP CHARGE=+1\n'

            compound = cc.Cartesian.read_xyz(str(xyz), start_index=1)
            with open(str(xyz)) as coord:
                coords = coord.readlines()
                atom_count = int(coords[0])
            connection_table = compound.get_bonds()
            site_connection = connection_table[site]
            for bonded in site_connection:
                if len(connection_table[bonded]) != 1 and bonded != site:
                    angle_connection = connection_table[bonded]
                    angle_atom = bonded

                else:
                    pass
            dihedral_atoms = []
            for dihed in angle_connection:
                if len(connection_table[dihed]) != 1 and dihed != angle_atom and dihed != site:
                    # dihedral_connection = connection_table[dihed]
                    dihedral_atoms.append(dihed)
                else:
                    pass
            dihedral_atom = min(dihedral_atoms)
            # print(f'site is: {site}, angle atom is: {angle_atom}, dihedral atom is: {dihedral_atom}')
            c1_dihed = 254.264918
            c2_dihed = 301.336871
            c3_dihed = 179.509647
            h1_dihed = 60.272731
            h2_dihed = 291.898120
            h3_dihed = 172.085792
            h4_dihed = 52.709817
            h5_dihed = 69.148928
            h6_dihed = 308.139965
            h7_dihed = 188.805597
            fukui = lines.copy()
            fukui_dir = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-1])}{str(Path("/fukui"))}'
            lines.append(f'C    12.021139  1  100.078262  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
            lines.append(f'C    1.482068  1  104.004890  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'C    1.481423  1  104.119367  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'H    1.103400  1   97.795054  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'H    1.117064  1  111.484379  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n')
            lines.append(f'H    1.120711  1  109.826665  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
            lines.append(f'H    1.118015  1  110.943566  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
            lines.append(f'H    1.117322  1  111.323645  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
            lines.append(f'H    1.118059  1  111.047346  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
            lines.append(f'H    1.120507  1  109.884160  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')

            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]

            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-9, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-8, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-7, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-6, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-5, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-4, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)]
            while any(check_atoms) is True:
                # print('here')
                rotation_spacing = 360/72
                c1_dihed += rotation_spacing
                c2_dihed += rotation_spacing
                c3_dihed += rotation_spacing
                h1_dihed += rotation_spacing
                h2_dihed += rotation_spacing
                h3_dihed += rotation_spacing
                h4_dihed += rotation_spacing
                h5_dihed += rotation_spacing
                h6_dihed += rotation_spacing
                h7_dihed += rotation_spacing
                lines[-10] = f'C    12.021139  1  100.078262  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
                lines[-9] = f'C    1.482068  1  104.004890  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-8] = f'C    1.481423  1  104.119367  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-7] = f'H    1.103400  1   97.795054  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-6] = f'H    1.117064  1  111.484379  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n'
                lines[-5] = f'H    1.120711  1  109.826665  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
                lines[-4] = f'H    1.118015  1  110.943566  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
                lines[-3] = f'H    1.117322  1  111.323645  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
                lines[-2] = f'H    1.118059  1  111.047346  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
                lines[-1] = f'H    1.120507  1  109.884160  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'

                os.remove(Path(f'{directory}/temp_ts.dat'))
                os.remove(Path(f'{directory}/temp_ts.xyz'))
                with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                    for string in lines:
                        tmp.write(str(string))

                check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz',
                            '-O', Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
                with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                    lin = inp.readlines()
                    cartesian = lin[2:]
                check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-9, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-8, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-7, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-6, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-5, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-4, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-3, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-2, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian)-1, directory),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian), directory)] # TODO RA: For loop?

        os.remove(Path(f'{directory}/temp_ts.dat'))
        os.remove(Path(f'{directory}/temp_ts.xyz'))
        os.remove(f'{input_file}')
        with open(f'{input_file}', 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))

        check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
        lines[-10] = f'C    50.000000  1  100.078262  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
        with open(Path(f'{directory}/HF_reagent.dat'), 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))
        with open(Path(f'{fukui_dir}/hf.dat'), 'w+') as f:
            for st in fukui:
                f.write(str(st))
        check_call(['obabel', '-imopin', Path(f'{directory}/HF_reagent.dat'), '-oxyz', '-O',
                    Path(f'{directory}/HF_reagent.xyz')], stdout=DEVNULL, stderr=STDOUT)
        check_call(['obabel', '-imopin', Path(f'{directory}/HF_reagent.dat'), '-osdf', '-O',
                    Path(f'{directory}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
        check_call(['obabel', '-imopin', Path(f'{fukui_dir}/hf.dat'), '-osdf', '-O',
                    Path(f'{fukui_dir}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
        hf_calcs_file = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-2])}{str(Path("/hf_jobs.txt"))}'
        if os.path.isfile(hf_calcs_file) is True:
            with open(hf_calcs_file, 'r') as fil:
                lins = fil.readlines()
            for li in lins:
                li = str(Path(li)).split(str(Path('/')))
                if not multiple:
                    if str(directory).split(str(Path('/')))[-2] != li[-3]:
                        print('REMOVING FILE')
                        os.remove(hf_calcs_file)
                        break
                else:
                    if str(directory).split(str(Path('/')))[-3] != li[-4]:
                        print('REMOVING FILE')
                        os.remove(hf_calcs_file)
                        break
            with open(hf_calcs_file, 'a') as fi:
                fi.write(str(Path(f'{directory}/hf.sdf\n')))
                fi.write(str(Path(f'{fukui_dir}/hf.sdf')))
        else:
            with open(hf_calcs_file, 'w') as fil:
                fil.write(str(Path(f'{directory}/hf.sdf\n')))
                fil.write(str(Path(f'{fukui_dir}/hf.sdf\n')))


def tweak_distance(input_file, radical, directory):
     # TODO RA: Docstring?
    os.remove(Path(f'{directory}/ts2.out'))
    with open(f'{input_file}', 'r+') as relaxed:
        lines = relaxed.readlines()
        if radical == 'cf3':
            tweak = lines[-4].split()
            tweak[1] = '1.900000'
            tweak.append('\n')
            t = "\t"
            t = t.join(tweak)
            lines[-4] = t

        elif radical == 'cf2h':
            tweak = lines[-4].split()
            tweak[1] = '1.900000'
            tweak.append('\n')
            t = "\t"
            t = t.join(tweak)
            lines[-4] = t

        elif radical == 'ipr':
            tweak = lines[-10].split()
            tweak[1] = '1.900000'
            tweak.append('\n')
            t = "\t"
            t = t.join(tweak)
            lines[-10] = t
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))
    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
    if os.path.isfile(Path(f'{directory}/ts2.arc')) is True:
        return
    else:
        with open(f'{input_file}', 'r+') as relaxed:
            lines = relaxed.readlines()
            if radical == 'cf3':
                tweak = lines[-4].split()
                tweak[1] = '2.100000'
                tweak.append('\n')
                t = "\t"
                t = t.join(tweak)
                lines[-4] = t

            elif radical == 'cf2h':
                tweak = lines[-4].split()
                tweak[1] = '2.100000'
                tweak.append('\n')
                t = "\t"
                t = t.join(tweak)
                lines[-4] = t

            elif radical == 'ipr':
                tweak = lines[-10].split()
                tweak[1] = '2.100000'
                tweak.append('\n')
                t = "\t"
                t = t.join(tweak)
                lines[-10] = t

    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))
    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)
    if os.path.isfile(Path(f'{directory}/ts2.arc')) is True:
        return
    else:
        print('Did not converge')
