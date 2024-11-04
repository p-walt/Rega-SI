import os
from pathlib import Path
import subprocess

with open('all_calcs.txt', 'r', encoding='utf-8') as file:
    directory = Path.cwd()
    dir_list = []
    dir_dict = {}
    for fi in os.listdir(directory):
        d = os.path.join(directory, fi)
        if os.path.isdir(d):
            dir_list.append(d)
    for line in file:
        line2 = line.strip('\n')
        splitted = line2.split('/')

        for dire in dir_list:
            if os.path.isfile(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/mm.xyz'))):
                if dire.split(str(Path('/')))[-1] in dir_dict:
                    dir_dict[dire.split(str(Path('/')))[-1]].append(line2)
                else:
                    dir_dict[dire.split(str(Path('/')))[-1]] = [line2]
    for folder in dir_dict.keys():
        with open(f'all_calcs1.txt', 'w') as rem:
            for calc in dir_dict[folder]:
                rem.write(f'{("/").join(calc.split("/")[:-1])}/hfoutp.in\n')
            # print('rsync',  '-a', f'--files-from={str(Path("./6a0b59a5-2c0e-46d9-9671-1ccb6c8bad3b/augusta_filelist"))}', f'./{dire}/', 'pcypw1@login001.augusta.nottingham.ac.uk:/gpfs01/home/pcypw1/22May/')
            #
            # outp = subprocess.Popen(['rsync',  '-a', f'--files-from={str(Path("./6a0b59a5-2c0e-46d9-9671-1ccb6c8bad3b/augusta_filelist"))}', f'./{dire}/', 'pcypw1@login001.augusta.nottingham.ac.uk:/gpfs01/home/pcypw1/22May/'],
            #                         stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            # print('command attempted')
            # status = outp.decode()
            # print(status)