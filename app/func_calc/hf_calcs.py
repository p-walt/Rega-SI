import os
import time
from pathlib import Path
from subprocess import DEVNULL, STDOUT, check_call
import subprocess
import re
import sys
from typing import Optional

# Define PyReact path dynamically or fall back to default if the environment variable is not set
PyReact_path = Path(os.getenv("PYREACT_PATH",
                              r'C:\\Users\\pcypw1\\Documents\\PhD\\GitHub\\ML-for-CH\\app\\PyReact_Scripts\\PyReact.py'))


def generate_hf_geometry(directory: Path, radical: str, failed: bool = False, multiple: bool = False) -> Optional[int]:
    """
    Generate Hartree-Fock optimized geometry files (HF) based on input data and directory structure.
    The function processes output files from prior calculations or optimization steps,
    extracts the required Cartesian coordinates, modifies them based on provided
    parameters, and creates various formatted geometry files for downstream tools.

    :param directory: The directory containing result files for Hartree-Fock geometry processing.
    :type directory: Path
    :param radical: A string defining the radical type (e.g., 'cf3', 'cf2h', 'ipr') to adjust
        specific atom coordinates in the generated file.
    :type radical: str
    :param failed: Flag indicating whether the optimization runs had failed or not.
        If True, alternative file processing paths will be used.
    :type failed: bool
    :param multiple: Flag to determine whether multiple file entries in the hf_jobs.txt
        should correspond to the processing directory.
    :type multiple: bool
    :return: An integer indicating the termination status or None for successful completion.
        Returns 2 if specific termination messages or errors are encountered within processed files.
    :rtype: Optional[int]
    """
    if failed is False:

        if os.path.isfile(Path(f'{directory}/ts2.out')):
            infile = Path(f'{directory}/ts2.out')
        elif os.path.isfile(Path(f'{directory}/ts.out')):
            infile = Path(f'{directory}/ts.out')
        elif os.path.isfile(Path(f'{directory}/am1.out')):
            infile = Path(f'{directory}/am1.out')
            # Read am1 output line by line but reverse order and put into new temporary file out_rev.txt.
        with open(infile) as rf, open(Path(f'{directory}/out_rev.txt'), 'w') as wf:
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
                elif 'CARTESIAN COORDINATES' in lines:
                    # collect block-related lines
                    while True:
                        try:
                            lines = next(file)
                        except StopIteration:
                            # there is no lines left
                            break
                        if '(I)' in lines:
                            # we've reached the end of block
                            break
                        a.append(lines)
                    # stop iterating over file
                    break
        # Reverse the list back to the right way round.
        a.reverse()
        # Remove the lines CARTESIAN COORDINATES and EIGENVALUES and the blank lines before and after the coordinates of
        # interest.
        trimmed = a[0:-1]
        coordinates = []
        # Remove the first column from the coordinates output of MOPAC which gives the atom number which is not needed.
        # This leaves ATOM LABEL, X, Y Z
        for atom in trimmed:
            strip = atom.strip('\n')
            line_elements = re.split(' +', str(strip))
            coordinate = line_elements[2:]
            co2 = []
            for co in coordinate:
                if co == '*':
                    co = '1'
                co2.append(co)
            if not coordinate:
                break
            else:
                del coordinate[0]

            coordinate1 = '\t'.join(map(str, co2))
            coordinates.append(coordinate1)
        # Write AM1 optimised coordinates to an xyz file
        with open(Path(f'{directory}/hf.dat'), 'w+') as dat:
            dat.write('PUT KEYWORDS HERE\n\n\n')
            for i in coordinates:
                dat.write(str(i) + '\n')
        os.remove(Path(f'{directory}/out_rev.txt'))
    else:

        if os.path.isfile(Path(f'{directory}/tmp2.out')):
            infile = Path(f'{directory}/tmp2.out')
        elif os.path.isfile(Path(f'{directory}/tmp.out')):
            infile = Path(f'{directory}/tmp.out')

        a = []
        # Scan the reversed file for the keyword EIGENVALUES and collect the following lines until CARTESIAN COORDINATES
        # is reached and add them to the list a.
        with open(infile) as file:
            for line in file:
                lines = line.strip('\n')
                if 'CALCULATION IS TERMINATED' in lines:
                    return 2
                elif 'Too many variables. By definition, at least one force constant is exactly zero' in lines:
                    return 2
                elif '(I)' in lines:
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

        # Remove the lines CARTESIAN COORDINATES and EIGENVALUES and the blank lines before and after the coordinates of
        # interest.
        trimmed = a[0:-1]
        coordinates = []
        # Remove the first column from the coordinates output of MOPAC which gives the atom number which is not needed.
        # This leaves ATOM LABEL, X, Y Z
        for atom in trimmed:
            strip = atom.strip('\n')
            line_elements = re.split(' +', str(strip))
            coordinate = line_elements[2:]
            co2 = []
            for co in coordinate:
                if co == '*':
                    co = '1'
                co2.append(co)
            if not coordinate:
                break
            else:
                del coordinate[0]

            coordinate1 = '\t'.join(map(str, co2))
            coordinates.append(coordinate1)
        # Write AM1 optimised coordinates to an xyz file
        with open(Path(f'{directory}/hf.dat'), 'w+') as dat:
            dat.write('PUT KEYWORDS HERE\n\n\n')
            for i in coordinates:
                dat.write(str(i) + '\n')


    # dihedral_atom = min(dihedral_atoms)
    with open(Path(f'{directory}/hf.dat')) as f:
        lines = f.readlines()
        if radical == 'cf3':
            carbon = lines[-4].split()
            carbon[1] = '2.120000\t'
            carbon.append('\n')
            t = " "
            t = t.join(carbon)
            lines[-4] = t
        elif radical == 'cf2h':
            carbon = lines[-4].split()
            carbon[1] = '2.120000\t'
            carbon.append('\n')
            t = "   "
            t = t.join(carbon)
            lines[-4] = t
        elif radical == 'ipr':
            carbon = lines[-10].split()
            carbon[1] = '2.120000\t'
            carbon.append('\n')
            t = " "
            t = t.join(carbon)
            lines[-10] = t
    os.remove(Path(f'{directory}/hf.dat'))
    with open(Path(f'{directory}/hf.dat'), 'w+') as new:
        for line in lines:
            new.write(str(line))
    check_call(['obabel', '-imopin', Path(f'{directory}/hf.dat'), '-oxyz',
               '-O', Path(f'{directory}/hf.xyz')], stdout=DEVNULL, stderr=STDOUT)
    check_call(['obabel', '-imopin', Path(f'{directory}/hf.dat'), '-osdf',
                '-O', Path(f'{directory}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
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
            print(f'DIRECTORY IS {directory}')
            fi.write(str(Path(f'{directory}/hf.sdf\n')))
    else:
        with open(hf_calcs_file, 'w') as fil:
            fil.write(str(Path(f'{directory}/hf.sdf\n')))

    # check_call('C:\\Users\\pcypw1\\Documents\\PhD\\PyReact_Scripts\\PyReact.py' )


def run_pyreact(command): # Could you not make pyReact a function and then just call the function... 
    """ Function that calls the Pyreact scripts using the subprocess module.
    @param command: The string of commands given to PyReact.
    @param files_to_run: The list of files for each compound/site to be calcualted on the HPC. Full path names required.
    @param first_run: True/False boolean that tells the function whether this is the first or second round of HF TS
    search.
    @return: The list of files to run and the dictionary containing the list of frequencies as the value."""
    print("Before Popen")

    # proc = subprocess.Popen(["pwd"],
    #                         stdout=subprocess.PIPE,
    #                         stderr=subprocess.PIPE, shell=True
    #                         )

    command2 = str(str(PyReact_path) + ' ' + command)
    print(command2)
    with open('test.log', 'wb') as f:
        # Subprocess Debugging
        # try:
        #     proc = subprocess.run([sys.executable] + command2.split(" "),
        #                           check=True,
        #                           capture_output=True,
        #                           text=True
        #                           )
        # except subprocess.CalledProcessError as error:
        #     print('Exception')
        #     print('output : ' + error.output)
        #     print('stderr : ' + error.stderr)
        proc = subprocess.Popen([sys.executable] + command2.split(" "), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        for c in iter(lambda: proc.stdout.read(1), b""):
            sys.stdout.buffer.write(c)
            f.write(c)
    print("After Popen")

    while True:
        poll = proc.poll()  # returns the exit code or None if the process is still running

        if poll is not None:
            break

    print(f"Final: {poll}")

    if poll is None:
        raise Exception("Timeout")


def check_hf_geometry(geometry: str, carbon_radical: str, site: int) -> bool:
    """
    Evaluates whether a given radical in a geometry file maintains a valid bond distance 
    with a specified site of interest. The function reads atomic coordinates from a file, 
    extracts the specified radical and site, and calculates the distance between them. 
    The calculation determines if the distance resides within an acceptable range 
    of bond lengths.

    :param geometry: Path to the file containing the atomic geometry information.
    :type geometry: str
    :param carbon_radical: Type of the carbon radical in the molecule. It can be cf3, 
        cf2h, or ipr.
    :type carbon_radical: str
    :param site: Index of the site of interest (1-based index) whose distance from 
        the carbon radical is to be calculated.
    :type site: int
    :return: True if the bond distance between the specified radical and site 
        lies between 1.8 and 2.2 units, otherwise False.
    :rtype: bool
    """
    with open(f'{geometry}', 'r+')as f:
        lin = f.readlines()
        coords = lin[2:]
        cart = []
        for atom in coords:
            strip = atom.strip('\n')
            line_elements = re.split(' +', str(strip))
            line = []

            for element in line_elements:
                line.append(element.split('\t'))
            flat_list = [item for sublist in line for item in sublist]
            cart.append(flat_list)
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


def hf_freq_check(freq_dict: dict, nwchem: bool) -> None:
    """
    Analyzes frequency data from Hartree-Fock calculations and determines necessary 
    actions such as geometry regeneration or pyreact simulation execution. 

    This function inspects a dictionary of frequency data and determines whether 
    to regenerate geometries or execute pyreact simulations based on a set of 
    conditions. It handles transition states, near-transition states, and failures, 
    and builds a list of files that need to be processed further. For certain cases, 
    this function also invokes external geometry regeneration routines or pyreact 
    simulations based on the computational suite used.

    :param freq_dict: A dictionary where keys are file paths (str) and values are lists 
        of frequency and energy values per respective file. The analysis is performed 
        based on these numerical values and file path properties.
    :type freq_dict: dict
    :param nwchem: Indicates whether the computations utilize the NWChem computational suite. 
        This boolean flag adapts the pyreact call based on the specific software being used.
    :type nwchem: bool
    :return: This function does not return any result. Instead, it indirectly operates by printing 
        log/diagnostic information, modifying file paths, and invoking geometry regeneration 
        or simulation commands as required.
    :rtype: None
    """
    redo_list = []
    print(freq_dict)
    for key, value in freq_dict.items():
        if type(value) == list:
            if key.split(str(Path("/")))[-2] == 'reagent':
                continue
            else:
                if - 800 <= float(value[0]) <= - 300:
                    if 0.00 <= float(value[1]):
                        print('Transition State')
                        redo_geometry_generator(key, nwchem, transition_state=True)
                        continue
                    elif - 100 <= float(value[1]) <= - 0.01:
                        if 0.00 <= float(value[2]):
                            if key.split(str(Path("/")))[-1] == 'hf2.sdf':
                                print('Close to TS')
                                redo_geometry_generator(key, nwchem, transition_state=True)
                                continue
                            else:
                                redo_list.append(key)
                                continue
                        elif - 100 <= float(value[2]) <= - 0.01:
                            if key.split(str(Path("/")))[-1] == 'hf2.sdf':
                                print('Fail')
                                continue
                            else:
                                redo_list.append(key)
                                continue
                    elif - 600 <= float(value[1]) <= - 101.1:
                        if key.split(str(Path("/")))[-1] == 'hf2outp.out':
                            print('Close to TS')
                            continue
                        else:
                            redo_list.append(key)
                            continue
                else:
                    print('Fail')
                    continue
    print(redo_list)
    if not redo_list:
        return
    else:
        file_list = []
        for file in redo_list:
            redo_geometry_generator(file, nwchem, transition_state=False)
            file_list.append(str(Path("/")).join(file.split(str(Path("/")))[:-1]))
        pyreact_files = " ".join([str(Path(f'{i}')) + str(Path('/hf2.sdf')) for i in file_list])
        if nwchem is False:
            run_pyreact(f'-d a --mult 2 --FC --NoEigen -F UHF --nProc 8 --TS --Rega --array {pyreact_files}',
                        pyreact_files.split(' '), first_run=False)
        else:
            run_pyreact(f'-d a --nwchem --mult 2 --FC --NoEigen -F UHF --nProc 8 --TS --Rega --array '
                        f'{pyreact_files}', pyreact_files.split(' '), first_run=False)


def redo_geometry_generator(out_file: str, nwchem: bool, transition_state: bool = False) -> None:
    """
    Function that processes the results of a Hartree-Fock (HF) computational chemistry
    calculation and generates geometry files for either the finalized transition state
    geometry or the next iteration of HF transition state search. The source file is
    analyzed based on the specified computational chemistry tool, and output files
    are created in specific formats.

    :param out_file: The full path to the file containing the HF calculation output data.
    :type out_file: str
    :param nwchem: Boolean value indicating whether the NWChem software package was used for the calculations.
    :type nwchem: bool
    :param transition_state: Specifies if the geometry output should represent the final Hartree-Fock
        (true) transition state or the starting point for the next stage of HF TS searching.
        Default value is False.
    :type transition_state: bool
    :return: This function does not return any value.
    :rtype: None
    """
    if nwchem is False:
        atom_dict = {'1': 'H', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl',
                     '35': 'Br',
                     '53': 'I'}
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hfoutp.out"))}', 'r') as rf, open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}', 'w') as wf:
            for line in reversed(rf.readlines()):
                wf.write(line)
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}', 'r') as file:
            a = []
            coordinates = []
            for line in file:
                lines = line.strip('\n')
                if 'Leave Link  202' in lines:
                    # collect block-related lines
                    while True:
                        try:
                            lines = next(file)
                        except StopIteration:
                            # there is no lines left
                            break
                        if 'Number' in lines:
                            # we've reached the end of block
                            break
                        a.append(lines)
                    # stop iterating over file
                    break
                # Remove the lines CARTESIAN COORDINATES and EIGENVALUES and the blank lines before and after the coordinates of
                # interest.
            a = a[1:-2]
            a.reverse()
            for atom in a:
                strip = atom.strip('\n')
                line_elements = re.split(' +', str(strip))
                element = line_elements[2]
                element = atom_dict[element]
                coord = line_elements[4:]
                co = []
                co.append(element)
                for coor in coord:
                    co.append(coor)
                coordinate = '\t'.join(map(str, co))
                coordinates.append(coordinate)
            system_size = len(a)
            if transition_state is False:
                with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.xyz"))}',
                          'w+') as xyz:
                    xyz.write(str(system_size) + '\n\n')
                    for i in coordinates:
                        xyz.write(str(i) + '\n')
                os.remove(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}')
                check_call(['obabel', '-ixyz',
                            f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.xyz"))}',
                            '-omopin',
                            '-O', f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.dat"))}'],
                           stdout=DEVNULL, stderr=STDOUT)
                check_call(['obabel', '-imopin',
                            f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.dat"))}',
                            '-osdf',
                            '-O', f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.sdf"))}'],
                           stdout=DEVNULL, stderr=STDOUT)
            else:
                with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf_ts.xyz"))}',
                          'w+') as xyz:
                    xyz.write(str(system_size) + '\n\n')
                    for i in coordinates:
                        xyz.write(str(i) + '\n')
                os.remove(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}')
    else:
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hfoutp.out"))}', 'r') as rf, open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}', 'w') as wf:
            for line in reversed(rf.readlines()):
                wf.write(line)

        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}', 'r') as file:
            a = []
            coordinates = []
            for line in file:
                lines = line.strip('\n')
                if 'Atomic Mass' in lines:
                    # collect block-related lines
                    while True:
                        try:
                            lines = next(file)
                        except StopIteration:
                            # there is no lines left
                            break
                        if 'No.       Tag' in lines:
                            # we've reached the end of block
                            break
                        a.append(lines)
                    # stop iterating over file
                    break
                # Remove the lines CARTESIAN COORDINATES and EIGENVALUES and the blank lines before and after the coordinates of
                # interest.
            a = a[1:-1]
            a.reverse()
            for atom in a:
                strip = atom.strip('\n')
                line_elements = re.split(' +', str(strip))
                atom_letter = str(line_elements[2])
                coordinate = line_elements[4:]
                coordinate.insert(0, atom_letter)

                if not coordinate:
                    break
                else:
                    del coordinate[0]

                coordinate1 = "\t".join(coordinate)
                coordinate2 = []
                coordinate2.append(coordinate1)
                coordinate2.insert(0, atom_letter)
                coordinate3 = "\t".join(coordinate2)
                coordinates.append(coordinate3)
        system_size = len(a)
        if transition_state is False:
            with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.xyz"))}', 'w+') as xyz:
                xyz.write(str(system_size) + '\n\n')
                for i in coordinates:
                    xyz.write(str(i) + '\n')
            os.remove(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}')
            check_call(['obabel', '-ixyz',
                        f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.xyz"))}', '-omopin',
                        '-O', f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.dat"))}'],
                       stdout=DEVNULL, stderr=STDOUT)
            check_call(['obabel', '-imopin',
                        f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.dat"))}', '-osdf',
                        '-O', f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.sdf"))}'],
                       stdout=DEVNULL, stderr=STDOUT)
        else:
            with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf_ts.xyz"))}', 'w+') as xyz:
                xyz.write(str(system_size) + '\n\n')
                for i in coordinates:
                    xyz.write(str(i) + '\n')
            os.remove(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}')
#             # frequencies = [item for sublist in frequencies for item in sublist]
#             # print(frequencies)
#             # print('First Frequency is: ' + str(frequencies[0]))
#             # freq_dict[out_file] = frequencies
