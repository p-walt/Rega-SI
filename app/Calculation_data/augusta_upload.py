from pathlib import Path
import os
import subprocess


# TODO RA: I would make this a function. And docstring this. 
with open('analysis_files.txt', 'r', encoding='utf-8') as file:
    directory = Path.cwd()

    dir_list = []
    for fi in os.listdir(directory):
        d = os.path.join(directory, fi)
        if os.path.isdir(d):
            dir_list.append(d.split(str(Path('/')))[-1])

    for dire in dir_list:# TODO RA: Hard-coded, needs fixing. 
        print('rsync',  '-a', f'--files-from={str(Path("./6a0b59a5-2c0e-46d9-9671-1ccb6c8bad3b/augusta_filelist"))}', f'./{dire}/', 'pcypw1@login001.augusta.nottingham.ac.uk:/gpfs01/home/pcypw1/25May/')

        outp = subprocess.Popen(['rsync',  '-a', f'--files-from={str(Path("./6a0b59a5-2c0e-46d9-9671-1ccb6c8bad3b/augusta_filelist"))}', f'./{dire}/', 'pcypw1@login001.augusta.nottingham.ac.uk:/gpfs01/home/pcypw1/25May/'],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print('command attempted')
        status = outp.decode()
        print(status)