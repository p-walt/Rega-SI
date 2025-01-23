# import csv
# import os
# from pathlib import Path
#
# directory = Path.cwd()
# with open('Reagents_redo.csv') as csv_file:
#     dir_list = []
#     dir_list2 = []
#     location_dict = {}
#     csv_reader = list(csv.reader(csv_file, delimiter=','))[1:]
#     # print(csv_reader)
#     for fi in os.listdir(directory):
#         d = os.path.join(directory, fi)
#         if os.path.isdir(d):
#             dir_list.append(d)
#     for dire in dir_list:
#         for com in os.listdir(dire):
#             c = os.path.join(dire, com)
#             if os.path.isdir(c):
#                 dir_list2.append(c)
#     for line in csv_reader:
#         compound_name = list(csv.reader(line, delimiter=','))[0][0]
#         # site = list(csv.reader(line, delimiter=','))[3][0]
#         smiles = list(csv.reader(line, delimiter=','))[1][0]
#         # print(compound_name)
#         for folder in dir_list2:
#             if os.path.isfile(Path(f'{folder}/{compound_name}.txt')):
#                 if os.path.isdir(Path(f'{folder}/reagent')):
#                     location = str(Path(f'{folder}/reagent'))
#                     print(location)
#                     loc = '/'.join(location.split('\\')[-3:])
#                     print(loc)
#                     if compound_name not in location_dict.keys():
#                         location_dict[compound_name] = [loc]
#                     else:
#                         location_dict[compound_name].append(str(loc))
#             elif os.path.isfile(Path(Path(f'{folder}/smiles.smi'))):
#                 with open(Path(Path(f'{folder}/smiles.smi')), 'r') as smile:
#                     test_smiles = smile.readlines()[0]
#                     if smiles == test_smiles:
#                         print(f'SMILES ARE {smiles}')
#                         print(f'TEST SMILES ARE {test_smiles}')
#                         if os.path.isdir(Path(f'{folder}/reagent')):
#                             location = str(Path(f'{folder}/reagent'))
#                             print(location)
#                             loc = '/'.join(location.split('\\')[-3:])
#                             print(loc)
#                             if compound_name not in location_dict.keys():
#                                 location_dict[compound_name] = [loc]
#                             else:
#                                 location_dict[compound_name].append(str(loc))
#     print(location_dict)
#     with open('Reagents_redo.txt', 'w') as file:
#         for key, value in location_dict.items():
#             if len(value) > 1:
#                 for va in value:
#                     file.write(f'{va}      {key}\n')
#             else:
#                 file.write(f'{value}      {key}\n')

# TODO RA: This needs cleaning up. Also consider making a function?

locations = []
with open('Reagents_redo.txt', 'r') as file:
    lines = file.readlines()
    for line in lines:
        split_line = line.split(' ')
        location = split_line[0]
        location1 = location.strip('[]')
        location2 = location1.replace("'", "")
        print(location2)
        locations.append(location2)
with open('Reagents_redo1.txt', 'w') as file1:
    for location in locations:
        file1.write(f'{location}\n')

# complete = 0
# incomplete_list = []
# complete_list = []
# with open('2nd_calc_location1.txt', 'r') as inp:
#     lines = inp.readlines()
#     for line in lines:
#         completed = False
#         location = line.split()[0]
#         if '[' in location:
#             folder = location[2:-2]
#         else:
#             folder = str(location)
#         if os.path.isfile(str(Path(f'{directory}/{folder}/hfoutp.out'))):
#             print('here')
#             with open(str(Path(f'{directory}/{folder}/hfoutp.out')), 'r', encoding='utf-8') as output:
#                 rev_file = reversed(output.readlines())
#                 lines = list(rev_file)
#                 for com in lines[0:5]:
#                     co = com.strip('\n')
#                     if 'Total times' in co:
#                         complete += 1
#                         completed = True
#                         complete_list.append(folder)
#                 if not completed:
#                     incomplete_list.append(folder)
#
#     with open('calc_redo1.txt', 'w') as file:
#         for li in incomplete_list:
#             file.write(f'{li}\n')

#     print(complete)
#     print(complete_list)
#     print(len(incomplete_list))
#     print(incomplete_list)
#
# with open('calc_redo.txt', 'r') as file:
#     lines = file.readlines()
#     li = []
#     for line in lines:
#         location = line.split()[0]
#         if '[' in location:
#             folder = location[2:-2]
#         else:
#             folder = location
#         li.append(folder)
# with open('calc_redo2.txt', 'w') as fi:
#     for l in li:
#         fi.write(f'{l}\n')
#
# import subprocess
# from pathlib import Path
# import os
#
# filelist = []
# sites_complete = []
# sites_failed = []
# sites_needing_work = []
#
# # with open('calc_redo1.txt', 'r') as inp:
# #     lines = inp.readlines()
# #     for line in lines:
# #         compound = '/'.join(line.split('/')[-2:])
# #         location = str(subprocess.check_output(f'find */{compound}/hfoutp.in', shell=True))
# #         filelist.append('/'.join(location.split('/')[:-1]))
# #     print(filelist)
#
# with open('32_jobs4.tmp', 'r') as inp:
#     lines = inp.readlines()
#     for line in lines:
#         location = '/'.join(line.split()[1].split('/')[4:])
#         filelist.append(location)
# complete = 0
# incomplete_list = []
# complete_list = []
# for line in filelist:
#     finished = False
#     if os.path.isfile(str(Path(f'{line}/hf2outp.out'))):
#         with open(str(Path(f'{line}/hf2outp.out')), 'r', encoding='utf-8') as out:
#             rev_file = reversed(out.readlines())
#             lines1 = list(rev_file)
#             for com in lines1[0:5]:
#                 co = com.strip('\n')
#                 if 'Total times' in co:
#                     print('finished')
#                     finished = True
#             if finished:
#                 lines = reversed(lines1)
#                 frequencies = []
#                 for li in lines:
#                     if 'P.Frequency' in li:
#                         tmp = ' '.join(li.split())
#                         new = tmp.split(' ')[1:]
#                         frequencies.append(new)
#                     else:
#                         continue
#                 print(frequencies)
#                 if frequencies:
#                     frequencies = [item for sublist in frequencies for item in sublist]
#
#                     print('LOOKING AT TS FREQUENCIES')
#                     if - 800 <= float(frequencies[0]) <= - 300:
#                         if 0.00 <= float(frequencies[1]):
#                             print('Transition State')
#                             # if line2 not in old_work_sites:
#                             #     if line2 not in old_comp_sites:
#                             sites_complete.append(line)
#                             continue
#                         elif - 100 <= float(frequencies[1]) <= - 0.01:
#                             if 0.00 <= float(frequencies[2]):
#                                 print('Close to transition state')
#                                 # if line2 not in old_work_sites:
#                                 #     if line2 not in old_comp_sites:
#                                 sites_needing_work.append(line)
#                                 continue
#                             elif - 100 <= float(frequencies[2]) <= - 0.01:
#                                 print('Fail')
#                                 # if line2 not in old_work_sites:
#                                 #     if line2 not in old_comp_sites:
#                                 sites_needing_work.append(line)
#                                 continue
#                         elif - 600 <= float(frequencies[1]) <= - 101.1:
#                             print('Fail')
#                             # if line2 not in old_work_sites:
#                             #     if line2 not in old_comp_sites:
#                             sites_needing_work.append(line)
#                             continue
#                     else:
#                         print('Fail')
#                         # if line2 not in old_work_sites:
#                         #     if line2 not in old_comp_sites:
#                         sites_needing_work.append(line)
#                         continue
#                 else:
#                     continue
#             else:
#                 sites_failed.append(line)
# with open('sites_complete2.txt', 'w') as comp_new:
#     for co in sites_complete:
#         comp_new.write(f'{co}\n')
#
# with open('sites_work2.txt', 'w') as wo_new:
#     for wo in sites_needing_work:
#         wo_new.write(f'{wo}\n')
#
# with open('sites_failed2.txt', 'w') as com:
#     for co in sites_failed:
#         com.write(f'{co}\n')
#
#         # with open(str(Path(f'{line}/hfoutp.out')), 'r') as output:
#         #    rev_file = reversed(output.readlines())
#         #    lines = list(rev_file)
#         #    for com in lines[0:5]:
#         #        co = com.strip('\n')
#         #        if 'Total times' in co:
#         #            complete += 1
#         #            complete_list.append(line)
#     # else:
#     #    incomplete_list.append(line)
# print(complete)
# print(complete_list)
# print(incomplete_list)
