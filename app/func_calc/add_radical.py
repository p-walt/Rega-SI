import os
from typing import Union

import chemcoord as cc
import re
from subprocess import DEVNULL, STDOUT, check_call
import glob, shutil
import func_calc.hf_calcs as hf
from app import pycharm
from pathlib import Path

RADICAL_PARAMS = {
    # CF3 parameters
    'cf3': {
        # Bond lengths
        'c_c_dist': 2.000000,  # C-C bond length
        'c_f1_dist': 1.351000,  # C-F1 bond length
        'c_f2_dist': 1.350000,  # C-F2 bond length
        'c_f3_dist': 1.345000,  # C-F3 bond length

        # Bond angles
        'c_c_angle': 101.000000,  # C-C-X angle
        'c_f1_angle': 107.100000,  # C-F1 angle
        'c_f2_angle': 107.800000,  # C-F2 angle
        'c_f3_angle': 114.600000,  # C-F3 angle

        # Dihedral angles
        'c_dihed': 180.0,  # Central carbon dihedral angle
        'f1_dihed': 180.0,  # First fluorine dihedral angle
        'f2_dihed': 60.0,  # Second fluorine dihedral angle
        'f3_dihed': -60.0,  # Third fluorine dihedral angle
    },

    # CF3 reactant parameters
    'cf3_reactant': {
        # Bond lengths
        'c_c_dist': 12.000000,  # C-C bond length (extended for reactant)
        'c_f1_dist': 1.351000,  # C-F1 bond length
        'c_f2_dist': 1.350000,  # C-F2 bond length
        'c_f3_dist': 1.345000,  # C-F3 bond length

        # Bond angles
        'c_c_angle': 101.000000,  # C-C-X angle
        'c_f1_angle': 107.100000,  # C-F1 angle
        'c_f2_angle': 107.800000,  # C-F2 angle
        'c_f3_angle': 114.600000,  # C-F3 angle

        # Dihedral angles
        'c_dihed': 180.0,  # Central carbon dihedral angle
        'f1_dihed': 180.0,  # First fluorine dihedral angle
        'f2_dihed': 60.0,  # Second fluorine dihedral angle
        'f3_dihed': -60.0,  # Third fluorine dihedral angle
    },

    # CF2H parameters
    'cf2h': {
        # Bond lengths
        'c_c_dist': 2.000000,  # C-C bond length
        'c_f1_dist': 1.351000,  # C-F1 bond length
        'c_f2_dist': 1.350000,  # C-F2 bond length
        'c_h_dist': 1.090000,  # C-H bond length

        # Bond angles
        'c_c_angle': 101.000000,  # C-C-X angle
        'c_f1_angle': 107.100000,  # C-F1 angle
        'c_f2_angle': 107.800000,  # C-F2 angle
        'c_h_angle': 111.200000,  # C-H angle

        # Dihedral angles
        'c_dihed': 180.0,  # Central carbon dihedral angle
        'f1_dihed': 120.0,  # First fluorine dihedral angle
        'f2_dihed': -120.0,  # Second fluorine dihedral angle
        'h_dihed': 0.0,  # Hydrogen dihedral angle
    },

    # CF2H reactant parameters
    'cf2h_reactant': {
        # Bond lengths
        'c_c_dist': 12.000000,  # C-C bond length (extended for reactant)
        'c_f1_dist': 1.351000,  # C-F1 bond length
        'c_f2_dist': 1.350000,  # C-F2 bond length
        'c_h_dist': 1.090000,  # C-H bond length

        # Bond angles
        'c_c_angle': 101.000000,  # C-C-X angle
        'c_f1_angle': 107.100000,  # C-F1 angle
        'c_f2_angle': 107.800000,  # C-F2 angle
        'c_h_angle': 111.200000,  # C-H angle

        # Dihedral angles
        'c_dihed': 180.0,  # Central carbon dihedral angle
        'f1_dihed': 120.0,  # First fluorine dihedral angle
        'f2_dihed': -120.0,  # Second fluorine dihedral angle
        'h_dihed': 0.0,  # Hydrogen dihedral angle
    },

    # Isopropyl (IPR) parameters
    'ipr': {
        # Bond lengths
        'c_c_dist': 2.021139,  # C-C bond length
        'c1_c_dist': 1.482068,  # C1-C bond length (first methyl)
        'c2_c_dist': 1.481423,  # C2-C bond length (second methyl)
        'h_c_dist': 1.103400,  # H-C bond length (central carbon hydrogen)
        'h1_c1_dist': 1.117064,  # H1-C1 bond length (first methyl)
        'h2_c1_dist': 1.120711,  # H2-C1 bond length (first methyl)
        'h3_c1_dist': 1.118015,  # H3-C1 bond length (first methyl)
        'h4_c2_dist': 1.117322,  # H4-C2 bond length (second methyl)
        'h5_c2_dist': 1.118059,  # H5-C2 bond length (second methyl)
        'h6_c2_dist': 1.120507,  # H6-C2 bond length (second methyl)

        # Bond angles
        'c_c_angle': 100.078262,  # Central C-C-X angle
        'c1_c_angle': 104.004890,  # C1-C-X angle (first methyl)
        'c2_c_angle': 104.119367,  # C2-C-X angle (second methyl)
        'h_c_angle': 97.795054,  # H-C-X angle (central carbon hydrogen)
        'h1_c1_angle': 111.484379,  # H1-C1-X angle (first methyl)
        'h2_c1_angle': 109.826665,  # H2-C1-X angle (first methyl)
        'h3_c1_angle': 110.943566,  # H3-C1-X angle (first methyl)
        'h4_c2_angle': 111.323645,  # H4-C2-X angle (second methyl)
        'h5_c2_angle': 111.047346,  # H5-C2-X angle (second methyl)
        'h6_c2_angle': 109.884160,  # H6-C2-X angle (second methyl)

        # Dihedral angles
        'c1_dihed': 180.0,  # Central carbon dihedral angle
        'c2_dihed': 60.0,  # First methyl carbon dihedral angle
        'c3_dihed': -60.0,  # Second methyl carbon dihedral angle
        'h1_dihed': 180.0,  # Central carbon hydrogen dihedral angle
        'h2_dihed': 60.0,  # First methyl first H dihedral angle
        'h3_dihed': -60.0,  # First methyl second H dihedral angle
        'h4_dihed': 180.0,  # First methyl third H dihedral angle
        'h5_dihed': 60.0,  # Second methyl first H dihedral angle
        'h6_dihed': -60.0,  # Second methyl second H dihedral angle
        'h7_dihed': 180.0,  # Second methyl third H dihedral angle
    },

    # Isopropyl (IPR) reactant parameters
    'ipr_reactant': {
        # Bond lengths
        'c_c_dist': 12.021139,  # C-C bond length (extended for reactant)
        'c1_c_dist': 1.482068,  # C1-C bond length (first methyl)
        'c2_c_dist': 1.481423,  # C2-C bond length (second methyl)
        'h_c_dist': 1.103400,  # H-C bond length (central carbon hydrogen)
        'h1_c1_dist': 1.117064,  # H1-C1 bond length (first methyl)
        'h2_c1_dist': 1.120711,  # H2-C1 bond length (first methyl)
        'h3_c1_dist': 1.118015,  # H3-C1 bond length (first methyl)
        'h4_c2_dist': 1.117322,  # H4-C2 bond length (second methyl)
        'h5_c2_dist': 1.118059,  # H5-C2 bond length (second methyl)
        'h6_c2_dist': 1.120507,  # H6-C2 bond length (second methyl)

        # Bond angles
        'c_c_angle': 100.078262,  # Central C-C-X angle
        'c1_c_angle': 104.004890,  # C1-C-X angle (first methyl)
        'c2_c_angle': 104.119367,  # C2-C-X angle (second methyl)
        'h_c_angle': 97.795054,  # H-C-X angle (central carbon hydrogen)
        'h1_c1_angle': 111.484379,  # H1-C1-X angle (first methyl)
        'h2_c1_angle': 109.826665,  # H2-C1-X angle (first methyl)
        'h3_c1_angle': 110.943566,  # H3-C1-X angle (first methyl)
        'h4_c2_angle': 111.323645,  # H4-C2-X angle (second methyl)
        'h5_c2_angle': 111.047346,  # H5-C2-X angle (second methyl)
        'h6_c2_angle': 109.884160,  # H6-C2-X angle (second methyl)

        # Dihedral angles
        'c1_dihed': 180.0,  # Central carbon dihedral angle
        'c2_dihed': 60.0,  # First methyl carbon dihedral angle
        'c3_dihed': -60.0,  # Second methyl carbon dihedral angle
        'h1_dihed': 180.0,  # Central carbon hydrogen dihedral angle
        'h2_dihed': 60.0,  # First methyl first H dihedral angle
        'h3_dihed': -60.0,  # First methyl second H dihedral angle
        'h4_dihed': 180.0,  # First methyl third H dihedral angle
        'h5_dihed': 60.0,  # Second methyl first H dihedral angle
        'h6_dihed': -60.0,  # Second methyl second H dihedral angle
        'h7_dihed': 180.0,  # Second methyl third H dihedral angle
    }
}




if pycharm == 1:
    mopac_path = "C:\\Program Files\\MOPAC\\MOPAC2016.exe"
else:
    mopac_path = 'mopac'


def read_mopac_output(in_file: Path, out_file: Path, directory: Path) -> int:
    """
    Reads a MOPAC AM1 calculation output file, processes the data to extract Cartesian
    coordinates, and writes the processed coordinates to an .xyz file.

    The function first reverses the MOPAC output file contents and scans for specific
    keywords to identify important blocks of data. It processes the data block to extract
    and format Cartesian coordinates and writes them to the specified output file. Temporary
    intermediary files are created and removed during this process. The function also handles
    specific error messages in the input file and provides appropriate return codes
    accordingly.

    :param in_file: Path to the input MOPAC output file.
    :type in_file: str
    :param out_file: Path to the output .xyz file containing formatted Cartesian coordinates.
    :type out_file: str
    :param directory: Directory where intermediary file(s) will be created and processed.
    :type directory: str
    :return: Status code indicating the success or error during processing. Returns ``1``
        for successful processing, ``2`` if specific termination errors are encountered.
    :rtype: int
    """

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


def mopac_distance_checks(geometry: Path, carbon_radical: str, site: int) -> bool:
    """
    Checks whether the distance between a specified carbon atom and another site of interest falls
    within a defined range, for a given geometry file and specific carbon radical type. The function
    reads atomic coordinates from a geometry file, identifies the atoms based on the given radical
    type and site index, computes their Euclidean distance, and evaluates whether it lies within
    the valid range (1.8 to 2.2 units). The function also handles missing geometry files gracefully
    and ensures compatibility with specific radical types.

    :param geometry: Path to the geometry file containing atomic coordinates. The file is expected
        to contain atomic data starting from the third line. Each line after the second represents
        an atom with its properties, including Cartesian coordinates.
    :type geometry: Path
    :param carbon_radical: Specifies the type of carbon radical. Valid options include 'cf3', 'cf2h',
        and 'ipr', each of which determines the specific carbon atom of interest for the distance
        calculation based on its unique position in the atom list.
    :type carbon_radical: str
    :param site: Index of the site of interest in the atomic list of the geometry file. The provided
        index determines which atom's coordinates are considered for the distance calculation relative
        to the selected carbon atom.
    :type site: int
    :return: `True` if the calculated distance between the specified carbon atom and the target site
        falls within the range of 1.8 and 2.2 units. Returns `False` if the distance is outside this
        range or the geometry file cannot be located.
    :rtype: bool
    """
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
    except FileNotFoundError:
        return False


def read_mopac_freq(input_file: Path, directory: Union[str, Path]) -> int:
    """
    Analyzes the frequency data from MOPAC-generated AM1 output files to determine the state of
    the system based on the number and magnitude of imaginary frequencies. The function reads
    the input file in reverse order, processes relevant frequency blocks, and identifies if the
    system is a true transition state, near transition state, or requires further intervention.

    :param input_file: Path to the MOPAC frequency output file.
    :type input_file: Path
    :param directory: Directory path where temporary files are created and processed.
    :type directory: Union[str, Path]
    :return: Status code indicating the interpretation of the frequency analysis:
        - 1: True transition state.
        - 2: Close to transition state with a small second imaginary frequency.
        - 3: Multiple negative frequencies or a large second imaginary frequency,
             requiring further intervention.
        - 4: No imaginary frequency in the correct range, needing other interventions.
    :rtype: int
    """

    # Read am1 output line by line but reverse order and put into new temporary file out_rev.txt.
    with open(input_file) as rf, open(Path(f'{directory}/out_rev.txt'), 'w') as wf:
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
    os.remove(Path(f'{directory}/out_rev.txt'))
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
        else:
            return 4
    except IndexError:
        return 4


def mopac_freq_check(inpu: Path, positively_charged: bool, radical: str, site: str, directory: Path,
                     new_mopac_2016: bool = False, multiple: bool = False) -> bool | None:
    """
    Determines if additional transition state searches are required by executing a frequency
    calculation on the resulting geometry. This function performs various steps including
    file manipulations, executing external commands, and interpreting the results. The decision
    is based on the nature of the compound, its charge state, and the calculation configurations.

    :param inpu: Path to the input geometry file in XYZ format.
    :type inpu: Path
    :param positively_charged: Boolean indicating whether the compound is positively charged.
    :type positively_charged: bool
    :param radical: Name of the radical functional group analyzed in the transition state calculation.
    :type radical: str
    :param site: Specific site under investigation for the transition state analysis.
    :type site: str
    :param directory: Parent directory where calculation and temporary files are stored.
    :type directory: Path
    :param new_mopac_2016: Boolean flag to use an updated version of MOPAC 2016 if available.
    :type new_mopac_2016: bool, optional
    :param multiple: Boolean indicating if multiple compounds are processed simultaneously.
    :type multiple: bool, optional
    :return: True if a transition state geometry or condition is determined; False otherwise.
    :rtype: bool
    """
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
                return None
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
                return None
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
    return None


def check_for_clashes(geometry: str | Path, query_atom: int, clash_threshold: float = 1) -> bool:
    """
    Checks for spatial clashes between atoms within a structure described in a geometry file. The
    function reads atomic coordinates from the given file, computes the inter-atomic distances
    between a specified query atom and all other atoms in the structure, and determines whether
    the distances fall below a specified clash threshold, indicating a potential clash.

    :param geometry: File path to the text geometry file containing atomic coordinates.
    :type geometry: str
    :param query_atom: Index of the atom of interest (starting from 1) whose distances to all other
        atoms will be calculated.
    :type query_atom: int
    :param clash_threshold: Threshold distance below which a clash is considered. Default is 1.
    :type clash_threshold: float, optional
    :return: Boolean value indicating whether a clash is detected. Returns True if any inter-atomic
        distance with the query atom is below the given threshold, otherwise False.
    :rtype: bool
    """
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
            cart.append(flat_list)
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


def cf3(input_file: str, xyz: str, site: int, directory: str, positively_charged: bool = False,
        constrained: bool = False, precise: bool = False, new_mopac_2016: bool = False) -> None:
    """
    This function modifies a MOPAC input file based on the specified parameters, generates temporary files for
    transition state geometry, and utilizes Open Babel and MOPAC to handle file conversions and computations. It
    uses geometrical constraints and connectivity information derived from the provided molecular structure file
    to ensure proper atom arrangements in the resulting output.

    :param input_file: Path to the MOPAC input file to be modified.
    :type input_file: str
    :param xyz: Path to the XYZ format file containing molecular structure information.
    :type xyz: str
    :param site: Index of the atom in the molecule serving as the central site for bond connections.
    :type site: int
    :param directory: Path to the directory where temporary files and output results will be stored.
    :type directory: str
    :param positively_charged: Boolean flag to indicate if the system involves a positively charged molecule or not.
    :type positively_charged: bool
    :param constrained: Boolean flag to determine if geometrical constraints should be applied or not.
    :type constrained: bool
    :param precise: Boolean flag to enable or disable the use of precise optimizations for calculations.
    :type precise: bool
    :param new_mopac_2016: Boolean flag to specify if MOPAC 2016-specific features/settings should be used.
    :type new_mopac_2016: bool
    :return: None
    :rtype: None
    """
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

        # Get parameters from the RADICAL_PARAMS dictionary
        params = RADICAL_PARAMS['cf3']

        # Extract necessary values
        c_c_dist = params['c_c_dist']
        c_f1_dist = params['c_f1_dist']
        c_f2_dist = params['c_f2_dist']
        c_f3_dist = params['c_f3_dist']

        c_c_angle = params['c_c_angle']
        c_f1_angle = params['c_f1_angle']
        c_f2_angle = params['c_f2_angle']
        c_f3_angle = params['c_f3_angle']

        c_dihed = params['c_dihed']
        f1_dihed = params['f1_dihed']
        f2_dihed = params['f2_dihed']
        f3_dihed = params['f3_dihed']

        lines.append(
            f'C    {c_c_dist}  {const}  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(
            f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(
            f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(
            f'F    {c_f3_dist}  1  {c_f3_angle}  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        # os.system('obabel -imopin temp_ts.dat -oxyz -O temp_ts.xyz > /dev/null 2>&1')
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]
        while any(check_atoms) is True:
            # print('here')
            rotation_spacing = 360 / 72
            c_dihed += rotation_spacing
            f1_dihed += rotation_spacing
            f2_dihed += rotation_spacing
            f3_dihed += rotation_spacing

            lines[
                -4] = f'C    {c_c_dist}  {const}  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[
                -3] = f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -2] = f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -1] = f'F    {c_f3_dist}  1  {c_f3_angle}  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

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
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]

    os.remove(Path(f'{directory}/temp_ts.dat'))
    os.remove(Path(f'{directory}/temp_ts.xyz'))
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)


def cf2h(input_file: str, xyz: str, site: int, directory: str, positively_charged: bool = False,
         constrained: bool = False, precise: bool = False, new_mopac_2016: bool = False) -> None:
    """
    This function modifies a MOPAC input file for CF2H radical, generates transition state geometry files,
    and handles file conversions and MOPAC computations with proper atom placement.

    :param input_file: Path to the MOPAC input file to be modified.
    :type input_file: str
    :param xyz: Path to the XYZ format file containing molecular structure.
    :type xyz: str
    :param site: Index of the atom serving as central site for bond connections.
    :type site: int
    :param directory: Path to the directory for temporary files and output.
    :type directory: str
    :param positively_charged: Flag for positively charged molecule.
    :type positively_charged: bool
    :param constrained: Flag for applying geometrical constraints.
    :type constrained: bool
    :param precise: Flag for precise optimizations.
    :type precise: bool
    :param new_mopac_2016: Flag for MOPAC 2016-specific features.
    :type new_mopac_2016: bool
    :return: None
    """
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

        # Get parameters from the RADICAL_PARAMS dictionary
        params = RADICAL_PARAMS['cf2h']

        # Extract necessary values
        c_c_dist = params['c_c_dist']
        c_f1_dist = params['c_f1_dist']
        c_f2_dist = params['c_f2_dist']
        c_h_dist = params['c_h_dist']

        c_c_angle = params['c_c_angle']
        c_f1_angle = params['c_f1_angle']
        c_f2_angle = params['c_f2_angle']
        c_h_angle = params['c_h_angle']

        c_dihed = params['c_dihed']
        f1_dihed = params['f1_dihed']
        f2_dihed = params['f2_dihed']
        h_dihed = params['h_dihed']

        lines.append(
            f'C    {c_c_dist}  {const}  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(
            f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(
            f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'H    {c_h_dist}  1  {c_h_angle}  1  {h_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]
        while any(check_atoms) is True:
            rotation_spacing = 360 / 72
            c_dihed += rotation_spacing
            f1_dihed += rotation_spacing
            f2_dihed += rotation_spacing
            h_dihed += rotation_spacing

            lines[
                -4] = f'C    {c_c_dist}  {const}  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[
                -3] = f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -2] = f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -1] = f'H    {c_h_dist}  1  {c_h_angle}  1  {h_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))
            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]

    os.remove(Path(f'{directory}/temp_ts.dat'))
    os.remove(Path(f'{directory}/temp_ts.xyz'))
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)


def ipr(input_file: str, xyz: str, site: int, directory: str, positively_charged: bool = False,
        constrained: bool = False, precise: bool = False, new_mopac_2016: bool = False) -> None:
    """
    This function modifies a MOPAC input file for isopropyl radical, generates transition state geometry files,
    and handles file conversions and MOPAC computations with proper atom placement.

    :param input_file: Path to the MOPAC input file to be modified.
    :type input_file: str
    :param xyz: Path to the XYZ format file containing molecular structure.
    :type xyz: str
    :param site: Index of the atom serving as central site for bond connections.
    :type site: int
    :param directory: Path to the directory for temporary files and output.
    :type directory: str
    :param positively_charged: Flag for positively charged molecule.
    :type positively_charged: bool
    :param constrained: Flag for applying geometrical constraints.
    :type constrained: bool
    :param precise: Flag for precise optimizations.
    :type precise: bool
    :param new_mopac_2016: Flag for MOPAC 2016-specific features.
    :type new_mopac_2016: bool
    :return: None
    """
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

        # Get parameters from the RADICAL_PARAMS dictionary
        params = RADICAL_PARAMS['ipr']

        # Extract necessary values
        c_c_dist = params['c_c_dist']
        c1_c_dist = params['c1_c_dist']
        c2_c_dist = params['c2_c_dist']
        h_c_dist = params['h_c_dist']
        h1_c1_dist = params['h1_c1_dist']
        h2_c1_dist = params['h2_c1_dist']
        h3_c1_dist = params['h3_c1_dist']
        h4_c2_dist = params['h4_c2_dist']
        h5_c2_dist = params['h5_c2_dist']
        h6_c2_dist = params['h6_c2_dist']

        c_c_angle = params['c_c_angle']
        c1_c_angle = params['c1_c_angle']
        c2_c_angle = params['c2_c_angle']
        h_c_angle = params['h_c_angle']
        h1_c1_angle = params['h1_c1_angle']
        h2_c1_angle = params['h2_c1_angle']
        h3_c1_angle = params['h3_c1_angle']
        h4_c2_angle = params['h4_c2_angle']
        h5_c2_angle = params['h5_c2_angle']
        h6_c2_angle = params['h6_c2_angle']

        c1_dihed = params['c1_dihed']
        c2_dihed = params['c2_dihed']
        c3_dihed = params['c3_dihed']
        h1_dihed = params['h1_dihed']
        h2_dihed = params['h2_dihed']
        h3_dihed = params['h3_dihed']
        h4_dihed = params['h4_dihed']
        h5_dihed = params['h5_dihed']
        h6_dihed = params['h6_dihed']
        h7_dihed = params['h7_dihed']

        lines.append(
            f'C    {c_c_dist}  {const}  {c_c_angle}  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(
            f'C    {c1_c_dist}  1  {c1_c_angle}  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(
            f'C    {c2_c_dist}  1  {c2_c_angle}  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'H    {h_c_dist}  1  {h_c_angle}  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(
            f'H    {h1_c1_dist}  1  {h1_c1_angle}  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n')
        lines.append(
            f'H    {h2_c1_dist}  1  {h2_c1_angle}  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
        lines.append(
            f'H    {h3_c1_dist}  1  {h3_c1_angle}  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
        lines.append(
            f'H    {h4_c2_dist}  1  {h4_c2_angle}  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
        lines.append(
            f'H    {h5_c2_dist}  1  {h5_c2_angle}  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
        lines.append(
            f'H    {h6_c2_dist}  1  {h6_c2_angle}  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 9),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 8),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 7),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 6),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 5),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 4),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]
        while any(check_atoms) is True:
            rotation_spacing = 360 / 72
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

            lines[
                -10] = f'C    {c_c_dist}  {const}  {c_c_angle}  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[
                -9] = f'C    {c1_c_dist}  1  {c1_c_angle}  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -8] = f'C    {c2_c_dist}  1  {c2_c_angle}  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -7] = f'H    {h_c_dist}  1  {h_c_angle}  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[
                -6] = f'H    {h1_c1_dist}  1  {h1_c1_angle}  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n'
            lines[
                -5] = f'H    {h2_c1_dist}  1  {h2_c1_angle}  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
            lines[
                -4] = f'H    {h3_c1_dist}  1  {h3_c1_angle}  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
            lines[
                -3] = f'H    {h4_c2_dist}  1  {h4_c2_angle}  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
            lines[
                -2] = f'H    {h5_c2_dist}  1  {h5_c2_angle}  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
            lines[
                -1] = f'H    {h6_c2_dist}  1  {h6_c2_angle}  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))
            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 9),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 8),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 7),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 6),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 5),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 4),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]

    os.remove(Path(f'{directory}/temp_ts.dat'))
    os.remove(Path(f'{directory}/temp_ts.xyz'))
    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)


def cf3_reactant(input_file: str, xyz: str, site: int, directory: str, positively_charged: bool = False,
                 reagent_optimised: bool = True, new_mopac_2016: bool = False, multiple: bool = False) -> None:
    """
    Generates the CF3 reactant structure, modifies atomic coordinates, and optimizes molecule files
    to handle geometric clashes while preparing reaction intermediates. It utilizes input structural
    data and parameters to create temporary and final data files required for chemical modeling.

    :param input_file: The path to the input file containing initial molecule data.
    :type input_file: str
    :param xyz: The xyz file path containing atomic coordinates of the molecule.
    :type xyz: str
    :param site: The index (0-based) of the atom site where the CF3 group is to be added.
    :type site: int
    :param directory: The directory path for output files and intermediate results.
    :type directory: str
    :param positively_charged: Specifies whether the molecule is positively charged. Defaults to False.
    :type positively_charged: bool
    :param reagent_optimised: Specifies whether the reactant has been already optimized. Defaults to True.
    :type reagent_optimised: bool
    :param new_mopac_2016: Specifies whether to use the new MOPAC 2016 parameterizations. Defaults to False.
    :type new_mopac_2016: bool
    :param multiple: Specifies whether multiple conformers or variations are processed. Defaults to False.
    :type multiple: bool
    :return: None
    :rtype: None
    """
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

        # Get parameters from the RADICAL_PARAMS dictionary
        params = RADICAL_PARAMS['cf3_reactant']

        # Extract necessary values
        c_c_dist = params['c_c_dist']
        c_f1_dist = params['c_f1_dist']
        c_f2_dist = params['c_f2_dist']
        c_f3_dist = params['c_f3_dist']

        c_c_angle = params['c_c_angle']
        c_f1_angle = params['c_f1_angle']
        c_f2_angle = params['c_f2_angle']
        c_f3_angle = params['c_f3_angle']

        c_dihed = params['c_dihed']
        f1_dihed = params['f1_dihed']
        f2_dihed = params['f2_dihed']
        f3_dihed = params['f3_dihed']

        fukui = lines.copy()
        fukui_dir = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-1])}{str(Path("/fukui"))}'

        lines.append(f'C    {c_c_dist}  1  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
        lines.append(f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
        lines.append(f'F    {c_f3_dist}  1  {c_f3_angle}  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

        with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
            for string in lines:
                tmp.write(str(string))

        check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                    Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
        with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
            lin = inp.readlines()
            cartesian = lin[2:]

        check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                       check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]
        while any(check_atoms) is True:
            rotation_spacing = 360/72
            c_dihed += rotation_spacing
            f1_dihed += rotation_spacing
            f2_dihed += rotation_spacing
            f3_dihed += rotation_spacing

            lines[-4] = f'C    {c_c_dist}  1  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
            lines[-3] = f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-2] = f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
            lines[-1] = f'F    {c_f3_dist}  1  {c_f3_angle}  1  {f3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

            os.remove(Path(f'{directory}/temp_ts.dat'))
            os.remove(Path(f'{directory}/temp_ts.xyz'))
            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]
            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]

        os.remove(Path(f'{directory}/temp_ts.dat'))
        os.remove(Path(f'{directory}/temp_ts.xyz'))
        # Set C-C distance to 50.0 for final structure, as in original code
        lines[-4] = f'C    50.000000  1  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
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


def cf2h_reactant(input_file: str, xyz: str, site: int, directory: str, positively_charged: bool = False,
                  reagent_optimised: bool = True, new_mopac_2016: bool = False, multiple: bool = False) -> None:
    """
    This function modifies and prepares a MOPAC input file to simulate a CF2H reactant. Depending on the provided parameters,
    it applies necessary optimizations, charges, and dihedral atom configurations while preparing and validating the molecular
    geometry. If `multiple` is True, it executes a reactive scan to resolve steric clashes. Throughout the process, intermediate
    files are written to a specified directory and cleaned up after the computation.

    :param input_file: Path to the MOPAC input file for the CF2H reactant.
    :param xyz: Path to the XYZ file containing the Cartesian coordinates of the initial molecular structure.
    :param site: Index of the specific atom in the molecule serving as the reaction site.
    :param directory: Directory to store intermediate and result files during computation.
    :param positively_charged: Flag indicating whether the reactant is positively charged. Default is False.
    :param reagent_optimised: Flag determining whether optimization steps are enabled. If True, skips MOPAC and uses existing geometry. Default is True.
    :param new_mopac_2016: Flag indicating whether to include MOPAC 2016 dispersion parameters in the calculation. Default is False.
    :param multiple: Flag to indicate whether to perform multiple rotation scans to prevent steric overlaps. Default is False.
    :return: None
    :rtype: None
    """
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

            # Get parameters from the RADICAL_PARAMS dictionary
            params = RADICAL_PARAMS['cf2h_reactant']

            # Extract necessary values
            c_c_angle = params['c_c_angle']
            c_f1_angle = params['c_f1_angle']
            c_f2_angle = params['c_f2_angle']
            c_h_angle = params['c_h_angle']

            c_dihed = params['c_dihed']
            f1_dihed = params['f1_dihed']
            f2_dihed = params['f2_dihed']
            h_dihed = params['h_dihed']

            # Use parameters for bond distances
            c_c_dist = params['c_c_dist']
            c_f1_dist = params['c_f1_dist']
            c_f2_dist = params['c_f2_dist']
            c_h_dist = params['c_h_dist']

            fukui = lines.copy()
            fukui_dir = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-1])}{str(Path("/fukui"))}'

            lines.append(f'C    {c_c_dist}  1  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
            lines.append(f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'H    {c_h_dist}  1  {c_h_angle}  1  {h_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')

            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]

            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]
            while any(check_atoms) is True:
                rotation_spacing = 360/72
                c_dihed += rotation_spacing
                f1_dihed += rotation_spacing
                h_dihed += rotation_spacing
                f2_dihed += rotation_spacing

                lines[-4] = f'C    {c_c_dist}  1  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
                lines[-3] = f'F    {c_f1_dist}  1  {c_f1_angle}  1  {f1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-2] = f'H    {c_h_dist}  1  {c_h_angle}  1  {h_dihed}   1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-1] = f'F    {c_f2_dist}  1  {c_f2_angle}  1  {f2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'

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
                check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]

        os.remove(Path(f'{directory}/temp_ts.dat'))
        os.remove(Path(f'{directory}/temp_ts.xyz'))
        os.remove(f'{input_file}')
        with open(f'{input_file}', 'w+') as file:
            for string in lines:
                # print(string)
                file.write(str(string))

        check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)

        # Set C-C distance to 50.0 for final structure as in original code
        lines[-4] = f'C    50.000000  1  {c_c_angle}  1  {c_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
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


def ipr_reactant(input_file: str, xyz: str, site: int, directory: str, positively_charged: bool = False,
                 reagent_optimised: bool = True, new_mopac_2016: bool = False, multiple: bool = False) -> None:
    """
    Processes the reaction information given input parameters, modifies files, computes geometric
    data, appends new molecular structures, checks for atom clashes during molecular rotations,
    and generates output files.

    :param input_file: Path to the input file containing initial molecular data.
    :type input_file: str
    :param xyz: Path to the .xyz file containing Cartesian coordinates of the molecule.
    :type xyz: str
    :param site: The index of the target site in the molecule for modification.
    :type site: int
    :param directory: Directory path where temporary and output files are to be stored.
    :type directory: str
    :param positively_charged: Indicates if the molecule has a positive charge.
    :type positively_charged: bool
    :param reagent_optimised: Indicates if the reagent geometry is already optimised.
    :type reagent_optimised: bool
    :param new_mopac_2016: Specifies whether to use MOPAC 2016 features such as dispersion.
    :type new_mopac_2016: bool
    :param multiple: Indicates if multiple molecules are being processed.
    :type multiple: bool
    :return: None. The function processes files and computes new structures without directly
        returning any value.
    :rtype: None
    """

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
                    dihedral_atoms.append(dihed)
                else:
                    pass
            dihedral_atom = min(dihedral_atoms)

            # Get parameters from the RADICAL_PARAMS dictionary
            params = RADICAL_PARAMS['ipr_reactant']

            # Extract necessary values
            c_c_dist = params['c_c_dist']  # This will be overridden with 12.0 and later 50.0
            c1_c_dist = params['c1_c_dist']
            c2_c_dist = params['c2_c_dist']
            h_c_dist = params['h_c_dist']
            h1_c1_dist = params['h1_c1_dist']
            h2_c1_dist = params['h2_c1_dist']
            h3_c1_dist = params['h3_c1_dist']
            h4_c2_dist = params['h4_c2_dist']
            h5_c2_dist = params['h5_c2_dist']
            h6_c2_dist = params['h6_c2_dist']

            c_c_angle = params['c_c_angle']
            c1_c_angle = params['c1_c_angle']
            c2_c_angle = params['c2_c_angle']
            h_c_angle = params['h_c_angle']
            h1_c1_angle = params['h1_c1_angle']
            h2_c1_angle = params['h2_c1_angle']
            h3_c1_angle = params['h3_c1_angle']
            h4_c2_angle = params['h4_c2_angle']
            h5_c2_angle = params['h5_c2_angle']
            h6_c2_angle = params['h6_c2_angle']

            c1_dihed = params['c1_dihed']
            c2_dihed = params['c2_dihed']
            c3_dihed = params['c3_dihed']
            h1_dihed = params['h1_dihed']
            h2_dihed = params['h2_dihed']
            h3_dihed = params['h3_dihed']
            h4_dihed = params['h4_dihed']
            h5_dihed = params['h5_dihed']
            h6_dihed = params['h6_dihed']
            h7_dihed = params['h7_dihed']

            fukui = lines.copy()
            fukui_dir = f'{str(Path("/")).join(str(directory).split(str(Path("/")))[:-1])}{str(Path("/fukui"))}'

            # Using 12.0 for initial C-C distance as in original code (with the exact values from the original code)
            lines.append(f'C    {c_c_dist}  1  {c_c_angle}  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n')
            lines.append(f'C    {c1_c_dist}  1  {c1_c_angle}  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'C    {c2_c_dist}  1  {c2_c_angle}  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'H    {h_c_dist}  1  {h_c_angle}  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n')
            lines.append(f'H    {h1_c1_dist}  1  {h1_c1_angle}  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n')
            lines.append(f'H    {h2_c1_dist}  1  {h2_c1_angle}  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
            lines.append(f'H    {h3_c1_dist}  1  {h3_c1_angle}  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n')
            lines.append(f'H    {h4_c2_dist}  1  {h4_c2_angle}  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
            lines.append(f'H    {h5_c2_dist}  1  {h5_c2_angle}  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')
            lines.append(f'H    {h6_c2_dist}  1  {h6_c2_angle}  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n')

            with open(Path(f'{directory}/temp_ts.dat'), 'w+') as tmp:
                for string in lines:
                    tmp.write(str(string))

            check_call(['obabel', '-imopin', Path(f'{directory}/temp_ts.dat'), '-oxyz', '-O',
                        Path(f'{directory}/temp_ts.xyz')], stdout=DEVNULL, stderr=STDOUT)
            with open(Path(f'{directory}/temp_ts.xyz'), 'r+') as inp:
                lin = inp.readlines()
                cartesian = lin[2:]

            check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 9),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 8),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 7),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 6),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 5),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 4),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                           check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]
            while any(check_atoms) is True:
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

                lines[-10] = f'C    {c_c_dist}  1  {c_c_angle}  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
                lines[-9] = f'C    {c1_c_dist}  1  {c1_c_angle}  1  {c2_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-8] = f'C    {c2_c_dist}  1  {c2_c_angle}  1  {c3_dihed}  1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-7] = f'H    {h_c_dist}  1  {h_c_angle}  1  {h1_dihed} 1     {atom_count + 1}  {site}   {angle_atom}\n'
                lines[-6] = f'H    {h1_c1_dist}  1  {h1_c1_angle}  1  {h2_dihed}   1     {atom_count + 3}  {atom_count + 1}   {site}\n'
                lines[-5] = f'H    {h2_c1_dist}  1  {h2_c1_angle}  1  {h3_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
                lines[-4] = f'H    {h3_c1_dist}  1  {h3_c1_angle}  1  {h4_dihed}  1     {atom_count + 3}  {atom_count + 1}   {site}\n'
                lines[-3] = f'H    {h4_c2_dist}  1  {h4_c2_angle}  1  {h5_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
                lines[-2] = f'H    {h5_c2_dist}  1  {h5_c2_angle}  1  {h6_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'
                lines[-1] = f'H    {h6_c2_dist}  1  {h6_c2_angle}  1  {h7_dihed}  1     {atom_count + 2}  {atom_count + 1}   {site}\n'

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
                check_atoms = [check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 9),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 8),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 7),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 6),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 5),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 4),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 3),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 2),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian) - 1),
                               check_for_clashes(Path(f'{directory}/temp_ts.xyz'), len(cartesian))]

        os.remove(Path(f'{directory}/temp_ts.dat'))
        os.remove(Path(f'{directory}/temp_ts.xyz'))
        os.remove(f'{input_file}')
        with open(f'{input_file}', 'w+') as file:
            for string in lines:
                file.write(str(string))

        check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)

        # Set C-C distance to 50.0 for final structure as in original code
        lines[-10] = f'C    50.000000  1  {c_c_angle}  1  {c1_dihed}  1     {site}   {angle_atom}   {dihedral_atom}\n'
        with open(Path(f'{directory}/HF_reagent.dat'), 'w+') as file:
            for string in lines:
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


def tweak_distance(input_file: Path, radical: str, directory: Path) -> None:
    """
    Updates bond length in a given molecular structure file, executes external computation, and verifies
    the output. Handles specific adjustments based on the given radical type and retries with modified
    parameters if initial computation does not succeed.

    The function performs the following steps:
    - Deletes a specified output file if it exists.
    - Reads the input file and updates bond length based on the specified radical.
    - Writes the modified input back to the file.
    - Executes an external tool via subprocess.
    - Verifies if the desired output file is generated, retrying with altered parameters if needed.
    - Prints an error message if the computation does not converge.

    :param input_file: The file path to the molecular structure input file.
    :type input_file: Path
    :param radical: The type of radical (e.g., 'cf3', 'cf2h', 'ipr') used to determine the bond adjustment.
    :type radical: str
    :param directory: The working directory containing output and temporary computation files.
    :type directory: Path
    :return: None
    :rtype: None
    """

    if radical not in RADICAL_PARAMS:
        print(f"Unsupported radical type: {radical}")
        return

    os.remove(Path(f'{directory}/ts2.out'))

    # First attempt with initial bond length
    with open(f'{input_file}', 'r+') as relaxed:
        lines = relaxed.readlines()
        line_index = RADICAL_PARAMS[radical]['line_offset']
        tweak = lines[line_index].split()
        tweak[1] = RADICAL_PARAMS[radical]['bond_lengths']['initial']
        tweak.append('\n')
        t = "\t"
        t = t.join(tweak)
        lines[line_index] = t

    os.remove(f'{input_file}')
    with open(f'{input_file}', 'w+') as file:
        for string in lines:
            # print(string)
            file.write(str(string))

    check_call([mopac_path, f'{input_file}'], stdout=DEVNULL, stderr=STDOUT)

    if os.path.isfile(Path(f'{directory}/ts2.arc')) is True:
        return
    else:
        # Retry with alternative bond length
        with open(f'{input_file}', 'r+') as relaxed:
            lines = relaxed.readlines()
            line_index = RADICAL_PARAMS[radical]['line_offset']
            tweak = lines[line_index].split()
            tweak[1] = RADICAL_PARAMS[radical]['bond_lengths']['retry']
            tweak.append('\n')
            t = "\t"
            t = t.join(tweak)
            lines[line_index] = t

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

