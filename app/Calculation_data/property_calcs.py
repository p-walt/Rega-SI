import os
import shutil
from pathlib import Path
import subprocess

with open('property_calculations.txt', 'r', encoding='utf-8') as file:
    lines = file.readlines()
    cwd = Path.cwd()
    for line in lines:
        stripped = line.strip('\n')
        splitted = stripped.split('/')
        print(splitted)

        with open(str(Path(f'{cwd}/{splitted[0]}/{splitted[1]}/{splitted[2]}/hfoutp.in')), 'r', encoding='utf-8')as old_file:
            rl = old_file.readlines()
            del rl[len(rl) - 11:]
            nl = ['dft\n', ' direct\n', ' xc b3lyp\n', ' fukui\n', ' print "Fukui information"\n', ' odft\n', 'end\n',
                  'property\n', ' nbofile\n', ' dipole\n', ' quadrupole\n', ' octupole\n', ' mulliken\n', ' esp\n',
                  ' efield\n', ' efieldgrad\n', ' efieldgradz4\n', ' gshift\n', ' electrondensity\n', 'end\n',
                  'set fock:mirrmat f\n', 'task dft property']

            rl.extend(nl)
        os.remove(str(Path(f'{cwd}/{splitted[0]}/{splitted[1]}/{splitted[2]}/hfoutp.in')))
        with open(str(Path(f'{cwd}/{splitted[0]}/{splitted[1]}/{splitted[2]}/hfoutp.in')), 'w') as new_file:
            for r in rl:
                new_file.write(r)


# with open('property_jobs.tmp', 'w') as jobs:
#     for line in lines:
#         stripped = line.strip('\n')
#         splitted = stripped.split('/')
#         jobs.write(f'cd /gpfs01/home/pcypw1/property_calcs/{splitted[0]}/{splitted[1]}/{splitted[2]} ; echo $PWD/hfoutp.in ; mpirun nwchem hfoutp.in > hfoutp.out\n')