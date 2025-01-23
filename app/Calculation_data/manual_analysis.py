import os
from pathlib import Path
import subprocess


# with open('sites_complete.txt', 'r', encoding='utf-8') as comp_ol:
#     old_comp_sites = comp_ol.readlines()
#
# with open('sites_work.txt', 'r', encoding='utf-8') as wo_ol:
#     old_work_sites = wo_ol.readlines()
# TODO RA: Make a function? 
with open('remote_filelist1.txt', 'r', encoding='utf-8') as file:
    directory = Path.cwd()

    dir_list = []
    for fi in os.listdir(directory):
        d = os.path.join(directory, fi)
        if os.path.isdir(d):
            dir_list.append(d)
    print(dir_list)
    sites_needing_work = []
    sites_complete = []
    for line in file:
        line2 = line.strip('\n')
        splitted = line2.split('/')
        print(splitted)

        for dire in dir_list:
            print(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/hfoutp.out')))
            if os.path.isfile(str(Path(f'{dire}/{str(Path("/")).join(splitted)}'))):
                print('FOUND IT')
                if os.path.isfile(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/hfoutp.out'))):
                    #     with open(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/hfoutp.out')), 'r', encoding='utf-8') as ou:
                    #         rl = ou.readlines()
                    #         fi = []
                    #         for ine in rl:
                    #             fi.append("".join(ine.split()))
                    #         print(fi[-100:-1])
                    #         if 'Pleasecitethefollowingreferencewhenpublishing' not in fi[-100:-1]:
                    #             print('didnt finish')
                    #             print('rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
                    #                                      'pcypw1@login.sulis.ac.uk:' + '/home/p/pcypw1/10May1622' + '/' + '/'.join(
                    #                                          splitted[:-1]) + '/hfoutp.out', str(Path('/')).join(['.', str(dire.split(str(Path("/")))[-1]), str(Path('/')).join(splitted[-5:-1])]) + '/hfoutp.out')
                    #             outp = subprocess.Popen(['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
                    #                                      'pcypw1@login.sulis.ac.uk:' + '/home/p/pcypw1/10May1622' + '/' + '/'.join(
                    #                                          splitted[:-1]) + '/hfoutp.out', str(Path('/')).join(['.', str(dire.split(str(Path("/")))[-1]), str(Path('/')).join(splitted[-5:-1])]) + '/hfoutp.out'],
                    #                                     stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    #             print('command attempted')
                    #             status = outp.decode()
                    #             print(status)
                    # if not os.path.isfile(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/hfoutp.out'))):
                    #     print('rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
                    #           'pcypw1@login.sulis.ac.uk:' + '/home/p/pcypw1/10May1622' + '/' + '/'.join(
                    #               splitted[:-1]) + '/hfoutp.out', str(Path('/')).join(['.', str(dire.split(str(Path("/")))[-1]),
                    #                                                                    str(Path('/')).join(
                    #                                                                        splitted[-5:-1])]) + '/hfoutp.out')
                    #     outp = subprocess.Popen(['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
                    #                              'pcypw1@login.sulis.ac.uk:' + '/home/p/pcypw1/10May1622' + '/' + '/'.join(
                    #                                  splitted[:-1]) + '/hfoutp.out', str(Path('/')).join(
                    #             ['.', str(dire.split(str(Path("/")))[-1]),
                    #              str(Path('/')).join(splitted[-5:-1])]) + '/hfoutp.out'],
                    #                             stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    #     print('command attempted')
                    #     status = outp.decode()
                    #     print(status)
                    # if os.path.isfile(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/hfoutp.out'))):
                    if splitted[-2] == 'reagent':
                        sites_complete.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                        continue
                    elif splitted[-2] == 'fukui':
                        sites_complete.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                        continue
                    else:
                        try:
                            with open(str(Path(f'{dire}/{str(Path("/")).join(splitted[:-1])}/hfoutp.out')), 'r',
                                      encoding='utf-8') as out:
                                lines = out.readlines()
                                frequencies = []
                                for li in lines:
                                    if 'P.Frequency' in li:
                                        tmp = ' '.join(li.split())
                                        new = tmp.split(' ')[1:]
                                        frequencies.append(new)
                                    else:
                                        continue
                                print(frequencies)
                                if frequencies:
                                    frequencies = [item for sublist in frequencies for item in sublist]
                                    # TODO RA: Maybe explain some of these if statements? the chemical backing maybe? 
                                    print('LOOKING AT TS FREQUENCIES') # TODO RA: What if freq less than -800?
                                    if - 800 <= float(frequencies[0]) <= - 300:
                                        if 0.00 <= float(frequencies[1]):
                                            print('Transition State')
                                            # if line2 not in old_work_sites:
                                            #     if line2 not in old_comp_sites:
                                            sites_complete.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                                            continue
                                        elif - 100 <= float(frequencies[1]) <= - 0.01: # TODO RA: Isnt this also a TS? 
                                            if 0.00 <= float(frequencies[2]):
                                                print('Close to transition state')
                                                # if line2 not in old_work_sites:
                                                #     if line2 not in old_comp_sites:
                                                sites_needing_work.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                                                continue
                                            elif - 100 <= float(frequencies[2]) <= - 0.01:
                                                print('Fail')
                                                # if line2 not in old_work_sites:
                                                #     if line2 not in old_comp_sites:
                                                sites_needing_work.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                                                continue
                                        elif - 600 <= float(frequencies[1]) <= - 101.1:
                                            print('Fail')
                                            # if line2 not in old_work_sites:
                                            #     if line2 not in old_comp_sites:
                                            sites_needing_work.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                                            continue
                                    else:
                                        print('Fail')
                                        # if line2 not in old_work_sites:
                                        #     if line2 not in old_comp_sites:
                                        sites_needing_work.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                                        continue
                                else:
                                    sites_needing_work.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))
                                    continue
                        except FileNotFoundError:
                            sites_needing_work.append(str(Path(f'{dire.split(str(Path("/")))[-1]}/{line2}')))

with open('sites_complete2.txt', 'w') as comp_new:
    for co in sites_complete:
        comp_new.write(f'{co}\n')

with open('sites_work2.txt', 'w') as wo_new:
    for wo in sites_needing_work:
        wo_new.write(f'{wo}\n')
