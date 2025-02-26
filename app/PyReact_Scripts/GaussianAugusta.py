#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:56:54 2017


@author: ke291

Contains all of the code specific to running Gaussian on Augusta.
Calls a lot of code from Gaussian.py  for input generation and calculation
execution. Called by PyReact.py.
"""

import os
import datetime
import re
import sys
from pathlib import Path
from glob import glob

try:
    import Gaussian
except ModuleNotFoundError:
    import PyReact_Scripts.Gaussian as Gaussian
import math
import subprocess
from subprocess import DEVNULL, STDOUT, check_call
from time import sleep

# Max concurrent jobs on Augusta
user = 'pcypw1'
CoresPerNode = 40
MaxCoresPerUser = 600
MaxTimeLimit = 168
MemPerCore = 2
sulis = False
if sulis:
    hostname = 'pcypw1@login.sulis.ac.uk'
    ssh_key = r"C:\Users\pcypw1\.ssh\sulis_key -p 22"
    rsync_key = r"-e 'ssh -i \Users\pcypw1\.ssh\sulis_key'"
    home = f'/home/p/{user}/'
else:
    hostname = 'pcypw1@hpclogin01.ada.nottingham.ac.uk'
    ssh_key = "-p 22"
    rsync_key = ''
    home = f'/gpfs01/home/{user}/'

rsync_path = "C:\\Program Files\\Git\\usr\\bin\\rsync.exe"
start = Path.cwd()

SetupGaussian = Gaussian.SetupGaussian

InitJobLog = Gaussian.InitJobLog

WriteJobLogFile = Gaussian.WriteJobLogFile

ReadJobLogFile = Gaussian.ReadJobLogFile

WriteGausFileRestart = Gaussian.WriteGausFileRestart

WriteGausFileOpt = Gaussian.WriteGausFileOpt

GetFiles2Run = Gaussian.GetFiles2Run

IsGausCompleted = Gaussian.IsGausCompleted

CheckOutputFiles = Gaussian.CheckOutputFiles

CheckConvergence = Gaussian.CheckConvergence

ReadSDFGeometry = Gaussian.ReadSDFGeometry

ReadGeometry = Gaussian.ReadGeometry

CopyGoutpFiles = Gaussian.CopyGoutpFiles


def RunCalcs(Files2Run, settings, molecules_dict, resubmission=False, resub_compound_and_sites=None,
             distance_tweak=False):

    print('\nRunning Gaussian on ' + hostname + '...')
    # Run Gaussian jobs on Nottingham Augusta cluster in folder named after date
    # and title and wait until the last file is completed
    now = datetime.datetime.now()
    if settings.nProc > CoresPerNode:
        print("Invalid number of cores requested, maximum is " + str(CoresPerNode))
        quit()

    if settings.TimeLimit > MaxTimeLimit:
        print("Too large job time limit requested, maximum is " + str(MaxTimeLimit))
        quit()

    MaxCon = math.floor(800 * (MaxCoresPerUser / settings.nProc))

    print('Maximum concurrent job limit set to ' + str(int(MaxCon)))

    if len(Files2Run) < MaxCon:
        print('---> ' + settings.Title)
        settings.folders.append(settings.now + settings.Title)
        if resubmission is False:

            RunOnAugusta(settings.now + settings.Title,
                         Files2Run, settings, molecules_dict)
        else:
            print('resubmitting')
            print(settings.now)
            print(settings.Title)
            print(Files2Run)
            print(settings)
            print(molecules_dict)
            print(resub_compound_and_sites)
            print(distance_tweak)
            QRes = ResubmitOnAugusta(settings.now + settings.Title, Files2Run, settings, molecules_dict,
                                     resub_compound_and_sites=resub_compound_and_sites, distance_tweak=distance_tweak)
            return QRes
    else:
        print("Job size/quantity exceeds limit, quitting... ")
        quit()


def RunOnAugusta(folder, GausJobs, settings, molecules_dict, wholecompound=False):
    """
    Function that executes calculations on the Augusta cluster.
    @param folder: The name of the remote folder where files will be uploaded/downloaded from and where calcualtions are
     executed.
    @param GausJobs: List of calculations to be run on Augusta
    @param settings: Settings Class
    @param molecules_dict: Dictionary with the hashed SMILES string of each compound as the key and the list of sites as
     the value.
    @param wholecompound: Whether the full compound is to be calculated or just a single calculation.
    @return: Error if folder already found on the cluster.
    """
    # Check that folder does not exist, create job folder on ziggy
    print("Checking if folder exists on " + hostname)
    print(f'ssh', f'{ssh_key}', f'{hostname}', f'ls')
    if sulis:
        outp = subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'ls'],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    else:
        outp = subprocess.Popen(['ssh', f'{hostname}', 'ls'],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    if folder in outp.decode():
        print("Folder exists on cluster, choose another folder name.")
        return
    fullfolder = home + folder
    print('ssh', f'{ssh_key}', f'{hostname}', 'mkdir', folder)
    if sulis:
        outp = subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'mkdir', folder],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    else:
        outp = subprocess.Popen(['ssh', f'{hostname}', 'mkdir', folder],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    if settings.Rega != '':
        compound_and_sites = {}
        for Gaus in GausJobs:
            compound = Gaus.split(str(Path("/")))[-3]
            if compound in compound_and_sites:
                site = Gaus.split(str(Path("/")))[-2]
                compound_and_sites[compound] += [site]
            else:
                compound_and_sites[compound] = []
                site = Gaus.split(str(Path("/")))[-2]
                compound_and_sites[compound] += [site]
        wholecompound = True
        print(compound_and_sites)
        for mol in compound_and_sites.keys():
            print('ssh', f'{ssh_key}', f'{hostname}', 'cd', fullfolder, ' ', 'mkdir', mol)
            if sulis:
                subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd', fullfolder, ' ', ';mkdir', mol], stderr=subprocess.STDOUT,
                                 stdout=subprocess.PIPE).communicate()[0]
            else:
                subprocess.Popen(
                    ['ssh', f'{ssh_key}', f'{hostname}', 'cd', fullfolder, ' ',
                     ';mkdir', mol], stderr=subprocess.STDOUT,
                    stdout=subprocess.PIPE).communicate()[0]
            fullfolder2 = fullfolder + '/' + mol
            for site in compound_and_sites[mol]:
                if sulis:
                    subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd', fullfolder2, ' ', ';mkdir', site], stderr=subprocess.STDOUT,
                                     stdout=subprocess.PIPE).communicate()[0]
                else:
                    subprocess.Popen(
                        ['ssh', f'{hostname}', 'cd', fullfolder2, ' ',
                         ';mkdir', site], stderr=subprocess.STDOUT,
                        stdout=subprocess.PIPE).communicate()[0]

    JobLogs = InitJobLog('/home/' + settings.user + '/' + folder, GausJobs, hostname, settings.nwchem)
    print(JobLogs)
    print(f'SETTINGS.INP_FILE {settings.inp_file}')
    # Write the slurm scripts
    SubFiles, ScratchFolders, JobLogs, i = WriteAugustaScripts(GausJobs, folder, settings, JobLogs, molecules_dict,
                                                            compound_and_sites, wholecompound)

    WriteJobLogFile(JobLogs)

    print(str(len(SubFiles)) + ' slurm scripts generated')

    # Upload .com files and slurm files to directory
    print("Uploading files to " + hostname + "...")
    if wholecompound is False:
        filelist = open('filelist', 'w')
        filelist.write('\n'.join(GausJobs))
        filelist.write('\n')

        filelist.write('\n'.join(SubFiles))
        filelist.write('\n')
        filelist.close()
    else:
        remote_filelist = open(
            str(Path('/')).join(settings.inp_file.split(str(Path("/")))[:-1]) + str(Path('/remote_filelist')), 'w')
        for jo in GausJobs:
            remote_filelist.write('/'.join(['/'.join(jo.split(str(Path("/")))[-3:-1]),
                                            jo.split(str(Path("/")))[-1].split('.')[0]]) + '.out' + '\n')
        remote_filelist.close()
        filelist_name = str(Path('/')).join(settings.inp_file.split(str(Path("/")))[:-1]) + str(Path('/filelist'))
        filelist = open(str(Path('/')).join(settings.inp_file.split(str(Path("/")))[:-1]) + str(Path('/filelist')), 'w')
        filelist.write('remote_filelist\n')
        filelist.write('\n')
        for file in SubFiles:
            filelist.write('/'.join(file.split(str(Path("/")))[-1:]) + '\n')
        filelist.close()
        filelist2 = open(str(Path('/')).join(settings.inp_file.split(str(Path("/")))[:-1]) + str(Path('/filelist2')), 'w')
        csv_dirs = []
        for job in GausJobs:
            csv_dir = job.split(str(Path("/")))[-4]
            if csv_dir not in csv_dirs:
                csv_dirs.append(csv_dir)
            filelist2.write('/'.join(job.split(str(Path("/")))[-3:]) + '\n')
        filelist2.close()
    if wholecompound is False:
        print('rsync', '-a', '--files-from=filelist', '.', hostname + ':' + fullfolder2 + '/')
        outp = subprocess.Popen(['rsync', '-a', '--files-from=filelist', '.', hostname + ':' + fullfolder2 + '/'],
                                stderr=subprocess.STDOUT,
                                stdout=subprocess.PIPE)
        output, error = outp.communicate()
        if outp.returncode != 0:
            print('ERROR' + str((outp.returncode, output, error)))
    else:
        remote_location = f'{hostname}:{fullfolder}'
        molecule_location = str(Path("/")).join(GausJobs[0].split(str(Path("/")))[1:-3])
        starting_loc = str(Path("/")).join(GausJobs[0].split(str(Path("/")))[6:-3])
        print(f'STARTING LOCATION {starting_loc}')
        print(str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1]))
        # outp = subprocess.Popen(
        #     ['cd', str(Path('/')).join(GausJobs[0].split("\\")[:-2]), ';rsync', '-a', '--files-from=filelist', '.',
        #      f'{hostname}:{fullfolder}/'], stderr=subprocess.STDOUT,
        #     stdout=subprocess.PIPE).communicate()[0]

        # subprocess.run(
        #     ['rsync', '-a', '--files-from=' + molecule_location + '\\filelist', molecule_location + '\\',
        #      remote_location], stderr=subprocess.STDOUT,
        #     stdout=subprocess.PIPE)
        print(Path.cwd())

        if not sulis:
            print('HERE', 'rsync', '-a', f'--files-from={str(Path("/"))}{molecule_location}{str(Path("/"))}filelist',
                 f'./{starting_loc}', remote_location)
            p = subprocess.Popen(
                ['rsync', '-a', f'--files-from={str(Path("/"))}{molecule_location}{str(Path("/"))}filelist',
                 f'./{starting_loc}', remote_location],
                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        else:

            print('rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a', f'--files-from=.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1])}{str(Path("/filelist"))}',
                  f'./{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1])}', remote_location)
            p = subprocess.Popen(
                ['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a', f'--files-from=.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1])}{str(Path("/filelist"))}',
                 f'./{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1])}', remote_location],
                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode != 0:
            print('ERROR' + str((p.returncode, output, error)))
    for di in csv_dirs:
        if sulis:

            print('rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
                  f'--files-from=.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1])}{str(Path("/filelist2"))}',
                  f'.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-2])}/{di}', remote_location)
            p = subprocess.Popen(
                ['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
                 f'--files-from=.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-1])}{str(Path("/filelist2"))}',
                 f'.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-5:-2])}/{di}', remote_location],
                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        else:
            print(Path.cwd())
            print(di)
            print('rsync', '-a', f'--files-from=.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-4:-1])}{str(Path("/filelist2"))}',
                  f'.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-4:-2])}/{di}',
                  remote_location)
            p = subprocess.Popen(
                ['rsync', '-a', f'--files-from=.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-4:-1])}{str(Path("/filelist2"))}',
                 f'.{str(Path("/"))}{str(Path("/")).join(settings.inp_file.split(str(Path("/")))[-4:-2])}/{di}',
                 remote_location],
                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode != 0:
            print('ERROR' + str((p.returncode, output, error)))
        # status = outp.decode().split('\n')[-2]
        # print(status + ' IS GOING WRONG')
    print(str(len(GausJobs)) + ' .com and ' + str(len(SubFiles)) + \
          ' slurm files uploaded to ' + hostname)

    """
    if settings.GuessChk:
        for f in GausJobs:
            outp = subprocess.Popen(['scp', f[:-4] + '.chk',
                                     'darwin:/home/' + settings.user + '/' + folder],
                                    stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print(str(len(GausJobs)) + ' .chk files uploaded to darwin')
    """
    """
    if (settings.KeepChk == True) and (settings.ResubsDone > 0):
        print('Copying .chk files to the new folder')
        for f in GausJobs:
            for c in settings.ChkFiles:
                if f[:-4] == c[0]:
                    print('ssh darwin cp /home/' + settings.user + '/' + settings.folders[-2] \
                          + '/' + c[1] + '.chk ' + fullfolder + '/')
                    outp = subprocess.Popen(['ssh', 'darwin', 'cp ' + '/home/' + settings.user + '/' \
                                             + settings.folders[-2] + '/' + c[1] + '.chk ' + fullfolder + '/', ], \
                                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    print('ssh darwin rm /home/' + settings.user + '/' + settings.folders[-2] \
                          + '/' + c[1] + '.chk')
                    outp = subprocess.Popen(['ssh', 'darwin', 'rm ' + '/home/' + settings.user + '/' \
                                             + settings.folders[-2] + '/' + c[1] + '.chk', ], \
                                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print(str(len(GausJobs)) + ' .chk files copied from previous folder on darwin, original files deleted.')
    """

    JobIDs = []
    batchfiles = []
    for f in SubFiles:
        if 'nwchem' in f.split(str(Path("/")))[-1]:
            batchfiles.append(f)
    # Launch the calculations
    for b in batchfiles:
        if wholecompound is False:
            print(' '.join(['ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + ';sbatch', f]))
            outp = subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + ';sbatch', b], \
                                    stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            status = outp.decode().split('\n')[-2]
            JobIDs.append(status.split('job ')[1])
        else:
            if not sulis:
                print('ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';dos2unix *')
                pre = subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';dos2unix *'],
                                       stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                print('ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch', b.split(str(Path("/")))[-1])
                outp = \
                subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch', b.split(str(Path("/")))[-1]],
                                 stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            else:
                print('ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ' + fullfolder + '/' + ';dos2unix *')
                pre = subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ' + fullfolder + '/' + ';dos2unix *'],
                                       stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                print('ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch',
                      b.split(str(Path("/")))[-1])
                outp = \
                    subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch',
                                      b.split(str(Path("/")))[-1]],
                                     stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

            status = outp.decode().split('\n')[-2]
            JobIDs.append(status.split('job ')[1])

    print(str(len(SubFiles)) + ' jobs submitted to the queue on ' + hostname)
    if settings.SubmitOnly:
        print('Submission only requested, exiting after successful submission of jobs to ' + hostname)
        quit()

    job_slurm_dict = {j: '' for j in GausJobs}
    sleep(60)
    OldQRes = CheckSlurmQueue(JobIDs)
    JobIDs = OldQRes[3]
    while OldQRes[0] < 0:
        OldQRes = CheckSlurmQueue(JobIDs)
        JobIDs = OldQRes[3]
        sleep(60)

    print('Pending: ' + str(OldQRes[0]) + ', Running: ' + str(OldQRes[1]) + ', Not in queue: ' + str(OldQRes[2]))

    TotalJobs = len(list(GausJobs))
    n2complete = TotalJobs
    files_time_dict = {}
    JobIDs_to_delete1 = []
    # Check and report on the progress of calculations
    while n2complete > 0:
        print('rsync', '-a', f'{remote_location}/slurm*',
              str(Path('/')).join(['.', str(Path('/')).join(GausJobs[0].split(str(Path('/')))[6:-3])]))
        if not sulis:
            outp = subprocess.Popen(['rsync', '-a', f'{remote_location}/slurm*', str(Path('/')).join(
                ['.', str(Path('/')).join(GausJobs[0].split(str(Path('/')))[6:-3])])],
                                    stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        else:
            outp = subprocess.Popen(['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a', f'{remote_location}/slurm*', str(Path('/')).join(
                ['.', str(Path('/')).join(GausJobs[0].split(str(Path('/')))[6:-3])])],
                                    stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        QRes = CheckSlurmQueue(JobIDs)
        print(QRes)
        JobIDs = QRes[3]
        job_slurm_dict = slurmjobpath(GausJobs[0], job_slurm_dict)
        files_time_dict, resub_jobs_dict1, xyz_del_jobs = write_int_xyz(GausJobs, fullfolder, files_time_dict, settings.nwchem, job_slurm_dict, settings,
                                        wholecompound)
        JobIDs_to_delete1.extend(xyz_del_jobs)
        print(JobIDs_to_delete1)
        for jo in JobIDs_to_delete1:
            print('Going through jobs to delete')
            if jo in QRes[3]:
                QRes[3].remove(jo)
                QRes[4].pop(jo)
        resub_jobs_dict2, JobIDs_to_delete, GausJobs = resubmission_check(GausJobs, jobs_dict=QRes[4], jobs_slurm_dict=job_slurm_dict, settings=settings,
                           compound_and_sites=compound_and_sites)

        resub_jobs_dict = resub_jobs_dict1.copy()
        print(resub_jobs_dict)
        for key, value in resub_jobs_dict2.items():
            TotalJobs += 1
            resub_jobs_dict[key] = value
        print(resub_jobs_dict)
        print(QRes)
        # Add new jobs to QRes array
        for key, value in resub_jobs_dict.items():
            print(key, value)
            QRes[3].append(key)
            print(QRes[3])
            QRes[4][key] = value
            print(QRes[4])
            if value == 'PD':
                QRes[0] += 1
            elif value == 'R':
                QRes[1] += 1
            elif value == '':
                QRes[2] += 1
        print(QRes)
        if QRes != OldQRes:
            if QRes[0] < 0:
                QRes = OldQRes
            else:
                OldQRes = QRes

                print('Pending: ' + str(OldQRes[0]) + ', Running: ' + str(OldQRes[1]) + ', Not in queue: ' + str(
                    OldQRes[2]))

        if (QRes[2] == len(JobIDs)) and (QRes[0] >= 0):
            # check each gaussian file to ascertain the status of individual gaus jobs
            print('No jobs left in ' + hostname + ' queue')
            break
        for jo in JobIDs_to_delete:
            print('Going through jobs to delete')
            if jo in QRes[3]:
                QRes[3].remove(jo)
                QRes[4].pop(jo)
        JobsFinished = IsAugustaGComplete(fullfolder, settings, compound_and_sites, wholecompound)
        JobsRemaining = TotalJobs - JobsFinished
        # If there are any resubmissions then add this number to the jumber of jobs remaining so the while loop
        # continues
        if JobsRemaining <= 0:
            if QRes[0] == 0:
                if QRes[1] == 0:
                    JobsRemaining = JobsRemaining
                elif QRes[1] != 0:
                    JobsRemaining = JobsRemaining + QRes[1]
            elif QRes[0] != 0:
                if QRes[1] == 0:
                    JobsRemaining = JobsRemaining + QRes[0]
                elif QRes[1] != 0:
                    JobsRemaining = JobsRemaining + QRes[0] + QRes[1]
        if n2complete != JobsRemaining:
            n2complete = JobsRemaining
        if QRes[0] == 0:
            if QRes[1] == 0:
                n2complete = 0
            print(str(n2complete) + " Gaussian jobs remaining.")
        print(f'JOBS LEFT {QRes[4]}')
        sleep(30)

    # When done, copy the results back
    print("\nCopying the output files back to localhost...")
    print('scp ' + hostname + ':' + fullfolder2 + '/*.out ' + str(start) + '/')

    # outp = subprocess.Popen(['scp', '-i', ssh_key, '-p', '22', hostname + ':' + fullfolder + '/*.out',
    #                          start], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    if not sulis:
        p = subprocess.Popen(
            ['rsync', '-a', f'--files-from={remote_location}/remote_filelist', f'{remote_location}', f'./{starting_loc}'],
            stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    else:
        p = subprocess.Popen(
            ['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a', f'--files-from={remote_location}/remote_filelist', f'{remote_location}',
             f'./{starting_loc}'],
            stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    if settings.KeepChk is True:
        print("\nCopying the checkpoint files back to localhost...")
        print('scp ' + hostname + ':' + fullfolder2 + '/*.chk ' + str(start) + '/')

        outp = subprocess.Popen(['scp', f'{ssh_key}', f'{hostname}' + ':' + fullfolder2 + '/*.chk',
                                 start], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    else:
        print("\nDeleting checkpoint files...")
        print('ssh ' + hostname + ' rm ' + fullfolder2 + '/*.chk')
        outp = subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'rm', fullfolder2 + '/*.chk'], \
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    if settings.Rega != '' and settings.nwchem is True:
        print("\nDeleting other NWChem files...")
        outp = \
        subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'rm', fullfolder2 + '/*/*.cphf* ' + fullfolder2
                          + '/*/*.db ' + fullfolder2 + '/*/*.fd* ' + fullfolder2 + '/*/*.hess ' + fullfolder2 +
                          '/*/*.movecs ' + fullfolder2 + '/*/*.nmode ' + fullfolder2 + '/*/*.stpr*'],
                         stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    if settings.WFN is True:
        print("\nCopying the wfn files back to localhost...")
        print('scp ' + hostname + ':' + fullfolder2 + '/*.wfn ' + start)

        outp = subprocess.Popen(['scp', 'darwin:' + fullfolder2 + '/*.wfn',
                                 start], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]


def WriteAugustaScripts(GausJobs, folder, settings, JobLogs, molecules_dict, compounds_and_sites, wholemolecule=False,
                        resubmission=False, distance_tweak=False):
    """Function that writes directs the writing of Slurm scripts to the correct function.
    @param GausJobs: List of file paths for each calculation to be done.
    @param folder: The folder to write the slurm scripts to.
    @param settings: Settings class
    @param JobLogs: List of jobs
    @param molecules_dict: Dictionary of hashed smiles string/molecule directory as key and the number of cores for each calculation as the value.
    @param compounds_and_sites: Dictionary of hashed smiles string/molecule directory as key and the list of sites within that molecule as the value.
    @param wholemolecule: Whether the wholemolecule is to be calculated or just a single site.
    @return: The list of names of the Slurm script files generated, The scratch folder and the JobLogs."""
    SubFiles = []
    ScratchFolders = []
    print('here!')
    print(molecules_dict)
    if settings.array is True:
        dire = settings.inp_file
        print(f'INP FILE IS {settings.inp_file}')
        cores8 = [k for k, v in molecules_dict.items() if v == 8]
        cores12 = [k for k, v in molecules_dict.items() if v == 12]
        cores16 = [k for k, v in molecules_dict.items() if v == 16]
        cores20 = [k for k, v in molecules_dict.items() if v == 20]
        cores24 = [k for k, v in molecules_dict.items() if v == 24]
        cores28 = [k for k, v in molecules_dict.items() if v == 28]
        cores32 = [k for k, v in molecules_dict.items() if v == 32]
        cores_dict = {8: cores8, 12: cores12, 16: cores16, 20: cores20, 24: cores24, 28: cores28, 32: cores32}
        for key in cores_dict.keys():
            if cores_dict[key]:
                if resubmission is False:
                    i = 0
                    jobs_list, bash_list = WriteSlurmArray(GausJobs, key, cores_dict, compounds_and_sites, folder, settings)
                    for j in jobs_list:
                        SubFiles.append(j)
                    for b in bash_list:
                        SubFiles.append(b)
                else:
                    i = 0
                    print('made it')
                    while os.path.exists(str(Path(
                            (str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{key}_jobs%s.tmp')) % i):
                        i += 1
                    WriteSlurmArray(GausJobs, key, cores_dict, compounds_and_sites, folder, settings, resubmission=True,
                                    i=i, distance_tweak=distance_tweak)
                    SubFiles.append(
                        str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{key}_jobs%s.tmp')) % i)
                    SubFiles.append(
                        str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{key}_nwchem%s.sh')) % i)
            else:
                continue
        if resubmission is False:
            arrayjob = open(str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + '/array_job.sh')), 'w')
            arrayjob.write('#!/bin/bash\n')
            if sulis is True:
                arrayjob.write('#SBATCH --account=su006-034\n')
            arrayjob.write('RUNLINE=$(cat $ARRAY_TASKFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)\n')
            arrayjob.write('eval $RUNLINE')
            arrayjob.close()
            SubFiles.append(str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + '/array_job.sh')))
    else:
        for GausJob in GausJobs:
            if not (os.path.exists(GausJob)):
                print("The input file " + GausJob + " does not exist. Exiting...")
                quit()

            SubFiles.append(WriteSlurmScript(GausJob.split(".")[0], folder, settings, molecules_dict, wholemolecule))

        ScratchFolders.append(GausJobs[0])
    return SubFiles, ScratchFolders, JobLogs, i


# Function to write slurm script
def WriteSlurmScript(GausJob, JobFolder, settings, molecules_dict, wholemolecule=False):
    """Function that writes the slurm bash script for single calculations
    @param GausJob: Full path of the file to be calculated.
    @param JobFolder: Folder the slurm script is written to.
    @param settings: Settings class.
    @param molecules_dict: Dictionary giving the number of cores/memory for this calculation.
    @param wholemolecule: Whether the wholemolecule is to be calculated or just a single site.
    @return: The name of the slurm script just written"""
    # Create the submission script
    SubScr = open(GausJob + "slurm", 'w')

    SubScr.write('#!/bin/bash\n\n')
    # Choose the queue/partition
    SubScr.write('#SBATCH --partition=defq\n')

    # Specify number of nodes (always 1), and number of cpus for this task
    if settings.Rega != '':
        SubScr.write('#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=' + str(
            molecules_dict[GausJob.split(str(Path("/")))[-3]]) + '\n')
    else:
        SubScr.write('#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=' + str(settings.nProc) + '\n')
    # Requested memory amount
    if settings.Rega != '':
        if molecules_dict[GausJob.split(str(Path("/")))[-3]] > 30:
            SubScr.write('#SBATCH --mem=' + str(molecules_dict[GausJob.split(str(Path("/")))[-3]] * MemPerCore) + 'g\n')
        else:
            SubScr.write('#SBATCH --mem=' + str(molecules_dict[GausJob.split(str(Path("/")))[-3]]) + 'g\n')
    else:
        SubScr.write('#SBATCH --mem=' + str(settings.nProc * MemPerCore) + 'g\n')
    # Specify time limit
    SubScr.write('#SBATCH --time=' + format(settings.TimeLimit, "02") + \
                 ':00:00\n\n')
    if settings.nwchem is False:
        # Purge and then load modules
        SubScr.write('echo "Running on `hostname`"\nmodule purge\n')
        SubScr.write('module load gaussian-uon/avx2/g16\n')
        # Create and export scratch folder
        if wholemolecule is False:
            SubScr.write('mkdir /gpfs01/home/' + user + '/scratch/' + settings.now + GausJob + '\n')
            SubScr.write('export GAUSS_SCRDIR=/gpfs01/home/' + user + '/scratch/' + settings.now + GausJob + '\n\n')
        else:

            SubScr.write('mkdir /gpfs01/home/' + user + '/scratch/' + settings.now + '\n')
            SubScr.write('mkdir /gpfs01/home/' + user + '/scratch/' + settings.now + '/' +
                         GausJob.split(str(Path('/')))[-2] + '\n')
            # SubScr.write('mkdir /gpfs01/home/' + user + '/scratch/' + settings.now + GausJob.split('/')[0] + '/' + GausJob.split('/')[-1] + '\n')
            SubScr.write('export GAUSS_SCRDIR=/gpfs01/home/' + user + '/scratch/' + settings.now + '/' +
                         GausJob.split(str(Path('/')))[-2] + '\n\n')

        # change into the job directory and start Gaussian
        SubScr.write('cd ${SLURM_SUBMIT_DIR}\n')
        if wholemolecule is False:
            SubScr.write('g16 < ' + GausJob + '.com > ' + GausJob + '.out\n\n')
        else:
            SubScr.write('g16 < ' + GausJob.split(str(Path('/')))[-1] + '.com > ' + GausJob.split(str(Path('/')))[
                -1] + '.out\n\n')
        SubScr.write('echo "Finished job now"\n')
    else:
        SubScr.write('module load nwchem-uoneasy/7.2.2-foss-2023.09\n\n')
        SubScr.write('cd $SLURM_SUBMIT_DIR\n\n')
        SubScr.write("export EXEC=`which nwchem`\n")
        SubScr.write('export INPUTFILE=' + GausJob.split(str(Path('/')))[-1] + '.in\n\n')
        SubScr.write('mpirun $EXEC $INPUTFILE > ' + GausJob.split(str(Path('/')))[-1] + '.out\n')
    # # Cleanup - delete the scratch folder
    # if wholemolecule is False:
    #     SubScr.write('rm -r /gpfs01/home/'+ user + '/scratch/' + settings.now + '/' + GausJob + '\n')
    # else:
    #     SubScr.write('rm -r /gpfs01/home/' + user + '/scratch/' + settings.now + '/' + GausJob.split(str(Path('/')))[-2] + '/' + GausJob.split(str(Path('/')))[-1] + '\n')
    #     SubScr.write('rm -r /gpfs01/home/' + user + '/scratch/' + settings.now + '\n')
    SubScr.close()

    return GausJob + "slurm"


def chunks(list_a, chunk_size):

    for i in range(0, len(list_a), chunk_size):
        yield list_a[i:i + chunk_size]


def WriteSlurmArray(GausJobs, core_count, cores_dict, compounds_and_sites, remote_folder, settings, resubmission=False, i=None, distance_tweak=False):
    """Function that writes the array script for each array for a given number of cores.
    @param GausJobs: List of calculations to be performed.
    @param core_count: The number of cores to give each of these calculations in this array.
    @param cores_dict: Dictionary with the core count as the key and the list of molecules requiring that core count as
    the value.
    @param compounds_and_sites: Dictionary with the hashed molecule SMILES string as the key and the list of sites as
    the value.
    @param remote_folder: The path of the destination on the remote server where calculations are going to be performed.
    """
    print('writing slurm array')
    chunked = []
    sh_list = []
    joby_list = []
    dire = settings.inp_file
    remote_folder = home + remote_folder
    if resubmission is False:
        job = str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_jobs.tmp'))
        joby_list.append(job)
    else:

        job = str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_jobs%s.tmp')) % i
        print('jobs written')
    jobs_list = []
    for molecule in cores_dict[core_count]:
        for site in compounds_and_sites[molecule]:
            for jo in GausJobs:
                if molecule == jo.split(str(Path('/')))[-3]:
                    if site == jo.split(str(Path('/')))[-2]:
                        job_of_interest = jo
            if not resubmission:
                jobs_list.append(f'cd {remote_folder}/{molecule}/{site} ; echo $PWD/hfoutp.in ; mpirun nwchem hfoutp.in > hfoutp.out')
            elif resubmission:
                if job_of_interest.split(str(Path('/')))[-1] == 'hfoutp.in':
                    jobs_list.append(f'cd {remote_folder}/{molecule}/{site} ; echo $PWD/hfoutp.in ; mpirun nwchem hfoutp.in > hfoutp.out')
                elif job_of_interest.split(str(Path('/')))[-1] == 'hf2outp.in':
                    jobs_list.append(f'cd {remote_folder}/{molecule}/{site} ; echo $PWD/hf2outp.in ; mpirun nwchem hf2outp.in > hf2outp.out')
    if len(jobs_list) <= 250:
        with open(job, mode='wt', encoding='utf-8') as myfile:
            myfile.write('\n'.join(jobs_list))
            myfile.write('\n')
    else:
        chunked = list(chunks(jobs_list, 250))

        i = 0
        for chunk in chunked:
            while os.path.exists(str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_jobs%s.tmp')) % i):
                i += 1
            joby = str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_jobs%s.tmp')) % i
            joby_list.append(joby)
            bash_path = str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_nwchem%s.sh')) % i
            sh_list.append(bash_path)
            with open(joby, 'w') as fils:
                fils.write('\n'.join(chunk))
                fils.write('\n')

    if core_count > 30:
        memory = core_count * 2
        time = 1080
    elif core_count > 22:
        memory = core_count
        time = 600
    else:
        memory = int(core_count / 2)
        time = 336
    print(f'MEMORY IS {memory}')
    if resubmission is False:
        if not chunked:
            shell = open(str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_nwchem.sh')),
                         'w')
            sh_list.append(str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_nwchem.sh')))
            shell.write('#!/bin/bash\n')
            if sulis is True:
                shell.write('#SBATCH --account=su006-034\n')
            shell.write('export ARRAY_JOBFILE=array_job.sh\n')
            shell.write(f'export ARRAY_TASKFILE={core_count}_jobs.tmp\n')
            shell.write('export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)\n')
            if not sulis:
                shell.write('module load nwchem-uoneasy/7.2.2-foss-2023.09\n')
            else:
                shell.write('module load iccifort/2019.5.281 OpenMPI/3.1.4\n')
                shell.write('module load NWChem/7.0.2-Python-3.7.4\n')
            shell.write('export EXEC=`which nwchem`\n')
            if sulis is True:
                shell.write(
                    f'CAL000=$(sbatch -J $EXEC --nodes=1 --ntasks-per-node={core_count} --mem={memory}g --account=su006-034 -t {time} --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)\n')
            else:
                shell.write(
                    f'CAL000=$(sbatch -J $EXEC --nodes=1  --ntasks-per-node={core_count} --mem={memory}g -p "defq" -t {time} --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)\n')
            shell.write('echo $CAL000')
            shell.close()
        else:
            for bas in sh_list:
                shell = open(bas, 'w')
                shell.write('#!/bin/bash\n')
                if sulis is True:
                    shell.write('#SBATCH --account=su006-034\n')
                shell.write('export ARRAY_JOBFILE=array_job.sh\n')
                shell.write(f'export ARRAY_TASKFILE={core_count}_jobs{sh_list.index(bas)}.tmp\n')
                shell.write('export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)\n')
                if not sulis:
                    shell.write('module load nwchem-uoneasy/7.2.2-foss-2023.09\n')
                else:
                    shell.write('module load iccifort/2019.5.281 OpenMPI/3.1.4\n')
                    shell.write('module load NWChem/7.0.2-Python-3.7.4\n')
                shell.write('export EXEC=`which nwchem`\n')
                if sulis is True:
                    shell.write(
                        f'CAL000=$(sbatch -J $EXEC --nodes=1 --ntasks-per-node={core_count} --mem={memory}g --account=su006-034 -t {time} --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)\n')
                else:
                    shell.write(
                        f'CAL000=$(sbatch -J $EXEC --nodes=1  --ntasks-per-node={core_count} --mem={memory}g -p "defq" -t {time} --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)\n')
                shell.write('echo $CAL000')
                shell.close()
    else:
        shell = open(
            str(Path((str(Path('/')).join(dire.split(str(Path('/')))[:-1])) + f'/{core_count}_nwchem%s.sh')) % i,
            'w')
        shell.write('#!/bin/bash\n')
        if sulis is True:
            shell.write('#SBATCH --account=su006-034\n')
        shell.write('export ARRAY_JOBFILE=array_job.sh\n')
        if resubmission is False:
            shell.write(f'export ARRAY_TASKFILE={core_count}_jobs.tmp\n')
        else:
            shell.write(f'export ARRAY_TASKFILE={core_count}_jobs{i}.tmp\n')
        shell.write('export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)\n')
        if not sulis:
            shell.write('module load nwchem-uoneasy/7.2.2-foss-2023.09\n')
        else:
            shell.write('module load iccifort/2019.5.281 OpenMPI/3.1.4\n')
            shell.write('module load NWChem/7.0.2-Python-3.7.4\n')
        shell.write('export EXEC=`which nwchem`\n')
        if sulis is True:
            shell.write(
                f'CAL000=$(sbatch -J $EXEC --nodes=1 --ntasks-per-node={core_count} --mem={memory}g --account=su006-034 -t {time} --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)\n')
        else:
            shell.write(
                f'CAL000=$(sbatch -J $EXEC --nodes=1  --ntasks-per-node={core_count} --mem={memory}g -p "defq" -t {time} --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)\n')
        shell.write('echo $CAL000')
        shell.close()
    print(f'JOBY LIST {joby_list}')
    print(f'SH LIST {sh_list}')
    return joby_list, sh_list

def CheckSlurmQueue(JobIDs):
    """Function that checks to see whether the calculqtions are still in the slurm queue. Also detects the new
    calculations that are started as part of the slurm array.
    @param JobIDs: The list of Slurm JobIDs that are assigned to each calculation
    @return: Array containing the number of calculations pending (being run), currently running, not in the slurm queue
    (and therefore complete) and the list of JobIDs back for when this function is called again"""
    if not sulis:
        outp = subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'squeue', '-u ' + user],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    else:
        outp = subprocess.Popen(
            ['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'squeue', '-u ' + user],
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    outp = outp.decode().split('\n')
    print(outp)
    QStart = -1000
    for i, line in enumerate(outp):
        if 'JOBID' in line:
            QStart = i + 1
            break

    if QStart < 0:
        return [-100, -100, -100, JobIDs]

    QueueReport = outp[QStart:-1]

    jobs_dict = {el: '' for el in JobIDs}
    for i, line in enumerate(QueueReport):
        Qjob = [_f for _f in line.split(' ') if _f][0]
        if Qjob in JobIDs:
            jobs_dict[Qjob] = [_f for _f in line.split(' ') if _f][4]
        elif '_' in Qjob:
            if int(Qjob.split('_')[0]) > int(JobIDs[0]):
                jobs_dict[Qjob] = [_f for _f in line.split(' ') if _f][4]
            else:
                continue

    PD = 'PD'
    R = 'R'
    NQ = ''
    JobIDs = list(jobs_dict.keys())
    Pending = sum(x == PD for x in jobs_dict.values())
    Running = sum(x == R for x in jobs_dict.values())
    NotInQueue = sum(x == NQ for x in jobs_dict.values())
    print(Pending, Running, NotInQueue, JobIDs, jobs_dict)
    return [Pending, Running, NotInQueue, JobIDs, jobs_dict]


def slurmjobpath(file, jobs_slurm_dict):
    slurm_file_location = str(Path((str(Path('/')).join(file.split(str(Path('/')))[:-3]))))
    slurm_files = glob(slurm_file_location + '/slurm*')
    for file in slurm_files:
        print(file)
        if os.path.getsize(file) != 0:
            with open(file, 'r+') as f:
                slurm_job_name = file.split(str(Path('/')))[-1]
                job_id = slurm_job_name.split('.')[0].split('-')[-1]
                remote_location = f.readlines()[0]
                print(job_id)
                print(remote_location)
                if '/' in remote_location:
                    molecule_and_site = '/'.join(remote_location.split('/')[-3:]).strip()
                    print(f'MOLECULE_AND_SITE {molecule_and_site}')
                    for key in jobs_slurm_dict.keys():
                        compound_site = '/'.join(key.split(str(Path('/')))[-3:]).strip()
                        print(f'COMPOUND SITE {compound_site}')
                        if compound_site == molecule_and_site:
                            print('MATCH')
                            jobs_slurm_dict[key] = job_id
                            print(jobs_slurm_dict)
                            break
        else:
            continue
    return jobs_slurm_dict


def resubmission_check(file_list, jobs_dict, jobs_slurm_dict, settings, compound_and_sites):
    resubmission_list = []
    resub_compounds_and_sites = {}
    JobIDs_to_delete = []
    print(f'JOBS DICT {jobs_dict}')
    print(f'JOBS SLURM DICT {jobs_slurm_dict}')
    for key in jobs_dict.keys():
        print(f'HERE: {key}')
        if '_' in key:
            if '[' not in key:
                if jobs_dict[key] == '':
                    print(f'FOUND SOMETHING')
                    try:
                        finished_compound_and_site = [k for k, v in jobs_slurm_dict.items() if v == key][0]
                        print(finished_compound_and_site)
                    except:
                        continue
                    if os.path.isfile(str(Path('/')).join(finished_compound_and_site.split(str(Path('/')))[:-1]) +
                                      str(Path('/hf2outp.out'))):
                        outfile = str(Path('/')).join(finished_compound_and_site.split(str(Path('/')))[:-1]) + \
                                  str(Path('/hf2outp.out'))
                    else:
                        outfile = str(Path('/')).join(finished_compound_and_site.split(str(Path('/')))[:-1]) + \
                                  str(Path('/hfoutp.out'))
                    print(outfile)
                    if settings.nwchem is False:
                        with open(outfile, 'r', encoding='utf-8') as f:
                            lines = f.readlines()
                            frequencies = []
                            for line in lines:
                                if 'Frequencies' in line:
                                    tmp = ' '.join(line.split())
                                    new = tmp.split(' ')[2:]
                                    frequencies.append(new)
                                else:
                                    continue
                            if frequencies:
                                frequencies = [item for sublist in frequencies for item in sublist]
                            else:
                                frequencies.append('-9999.999')
                            print(frequencies)
                            print('First Frequency is: ' + str(frequencies[0]))
                    else:
                        with open(outfile, 'r', encoding='utf-8') as f:
                            lines = f.readlines()
                            frequencies = []
                            for line in lines:
                                if 'P.Frequency' in line:
                                    tmp = ' '.join(line.split())
                                    new = tmp.split(' ')[1:]
                                    frequencies.append(new)
                                else:
                                    continue
                            if frequencies:
                                frequencies = [item for sublist in frequencies for item in sublist]
                            else:
                                frequencies.append('-9999.999')
                            print(frequencies)
                            print('First Frequency is: ' + str(frequencies[0]))
                    print('IM HERE')
                    if finished_compound_and_site.split(str(Path('/')))[-2] == 'reagent':
                        JobIDs_to_delete.append(key)
                        continue
                    elif finished_compound_and_site.split(str(Path('/')))[-2] == 'fukui':
                        JobIDs_to_delete.append(key)
                        continue

                    else:
                        print('LOOKING AT TS FREQUENCIES')
                        print(finished_compound_and_site.split(str(Path("/")))[-1])
                        if - 800 <= float(frequencies[0]) <= - 300:
                            if 0.00 <= float(frequencies[1]):
                                print('Transition State')
                                redo_geometry_generator(finished_compound_and_site, settings.nwchem, transition_state=True)
                                JobIDs_to_delete.append(key)
                                continue
                            elif - 100 <= float(frequencies[1]) <= - 0.01:
                                if 0.00 <= float(frequencies[2]):
                                    if finished_compound_and_site.split(str(Path("/")))[-1] == 'hf2outp.in':
                                        print('Close to TS')
                                        redo_geometry_generator(finished_compound_and_site, settings.nwchem,
                                                                transition_state=True)
                                        JobIDs_to_delete.append(key)
                                        continue
                                    else:
                                        redo_geometry_generator(finished_compound_and_site, settings.nwchem,
                                                                transition_state=False)
                                        resubmission_list.append(
                                            f'{str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1])}{str(Path("/hf2.sdf"))}')
                                        print(file_list)
                                        print(finished_compound_and_site)
                                        if finished_compound_and_site in file_list:
                                            print('removing old file from file_list')
                                            file_list.remove(finished_compound_and_site)
                                        file_list.append(
                                            f'{str(Path("/")).join([str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1]), str(Path("hf2outp.in"))])}')
                                        jobs_slurm_dict[str(Path("/")).join([str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1]), str(Path("hf2outp.in"))])] = ''
                                        resub_compounds_and_sites[finished_compound_and_site] = key
                                        JobIDs_to_delete.append(key)
                                        continue
                                elif - 100 <= float(frequencies[2]) <= - 0.01:
                                    if finished_compound_and_site.split(str(Path("/")))[-1] == 'hf2outp.in':
                                        print('Fail')
                                        JobIDs_to_delete.append(key)
                                        continue
                                    else:
                                        redo_geometry_generator(finished_compound_and_site, settings.nwchem,
                                                                transition_state=False)
                                        resubmission_list.append(
                                            f'{str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1])}{str(Path("/hf2.sdf"))}')
                                        print(finished_compound_and_site)
                                        if finished_compound_and_site in file_list:
                                            print('removing old file from file_list')
                                            file_list.remove(finished_compound_and_site)
                                        file_list.append(
                                            f'{str(Path("/")).join([str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1]), str(Path("hf2outp.in"))])}')
                                        jobs_slurm_dict[str(Path("/")).join(
                                            [str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1]),
                                             str(Path("hf2outp.in"))])] = ''
                                        resub_compounds_and_sites[finished_compound_and_site] = key
                                        JobIDs_to_delete.append(key)
                                        continue
                            elif - 600 <= float(frequencies[1]) <= - 101.1:
                                if finished_compound_and_site.split(str(Path("/")))[-1] == 'hf2outp.in':
                                    print('Close to TS')
                                    redo_geometry_generator(finished_compound_and_site, settings.nwchem,
                                                            transition_state=True)
                                    JobIDs_to_delete.append(key)
                                    continue
                                else:
                                    redo_geometry_generator(finished_compound_and_site, settings.nwchem,
                                                            transition_state=False)
                                    resubmission_list.append(
                                        f'{str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1])}{str(Path("/hf2.sdf"))}')
                                    if finished_compound_and_site in file_list:
                                        print('removing old file from file_list')
                                        file_list.remove(finished_compound_and_site)
                                    file_list.append(
                                        f'{str(Path("/")).join([str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1]), str(Path("hf2outp.in"))])}')
                                    jobs_slurm_dict[str(Path("/")).join(
                                        [str(Path("/")).join(finished_compound_and_site.split(str(Path("/")))[:-1]),
                                         str(Path("hf2outp.in"))])] = ''
                                    resub_compounds_and_sites[finished_compound_and_site] = key
                                    JobIDs_to_delete.append(key)
                                    continue
                        else:
                            print('Fail')
                            JobIDs_to_delete.append(key)
                            continue
            else:
                continue
        else:
            continue
    if resub_compounds_and_sites:
        for r in resubmission_list:
            if r.split(str(Path('/')))[-3] in resub_compounds_and_sites:
                if r.split(str(Path('/')))[-2] in resub_compounds_and_sites[r.split(str(Path('/')))[-3]]:
                    continue
                else:
                    resub_compounds_and_sites[r.split(str(Path('/')))[-3]] += [r.split(str(Path('/')))[-2]]
            else:
                resub_compounds_and_sites[r.split(str(Path('/')))[-3]] = [r.split(str(Path('/')))[-2]]
        resub_GausInpFiles, resub_molecules_dict = Gaussian.SetupGaussian(resubmission_list, settings, wholemolecule=True
                                                                          , resubmission=True)
    else:
        return {}, JobIDs_to_delete, file_list
    QRun = False
    Files2Run = Gaussian.GetFiles2Run(resub_GausInpFiles, settings)
    print(f'FILES TO RUN:{Files2Run}')
    print(file_list)
    if len(Files2Run) == 0:
        QRun = True

    if settings.GenOnlyQ:
        print("Gaussian input files generated, quitting as instructed.")
        quit()

    if QRun:
        print('DFT has already been run for these inputs. Skipping...')
    else:
        print('running retry')
        QRes = RunCalcs(Files2Run, settings, resub_molecules_dict, resubmission=True,
                 resub_compound_and_sites=resub_compounds_and_sites, distance_tweak=False)
        return QRes[4], JobIDs_to_delete, file_list


def redo_geometry_generator(out_file, nwchem, transition_state=False):
    """Function that creates either the finished transition state geometry or the geoemtry for the next round of HF TS
    search from the HF output file.
    @param out_file: The full file path of the HF output file.
    @param nwchem: Boolean on whether the NWChem software package is being used.
    @param transition_state: Boolean on whether the geometry should be labelled as the true HF transition state or a
    starting point for the second round of HF TS searching with tightened convergence criteria."""
    if nwchem is False:
        atom_dict = {'1': 'H', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl',
                     '35': 'Br',
                     '53': 'I'}
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hfoutp.out"))}', 'r',
                  encoding='utf-8') as rf:
            print('opening output file')
            rev_lines = reversed(rf.readlines())
            print(rev_lines)
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}', 'w',
                  encoding='utf-8') as wf:
            print('opening reversed_file')
            for line in rev_lines:
                wf.write(line)
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}',
                  'r', encoding='utf-8') as file:
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
                          'w+', encoding='utf-8') as xyz:
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
                print('written redo files')
            else:
                with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf_ts.xyz"))}',
                          'w+', encoding='utf-8') as xyz:
                    xyz.write(str(system_size) + '\n\n')
                    for i in coordinates:
                        xyz.write(str(i) + '\n')
                os.remove(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}')
    else:
        try:
            with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2outp.out"))}',
                      'r', encoding='utf-8') as rf:
                rev_lines = reversed(rf.readlines())
        except FileNotFoundError:
            with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hfoutp.out"))}',
                      'r', encoding='utf-8') as rf:
                rev_lines = reversed(rf.readlines())
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}', 'w', encoding='utf-8') as wf:
            print('opening reversed_file')
            for line in rev_lines:
                wf.write(line)
            print('writing reversed file')
        print('written reversed file')
        with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}',
                  'r', encoding='utf-8') as file:
            print('reading reversed file')
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
            with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf2.xyz"))}',
                      'w+', encoding='utf-8') as xyz:
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
            print('written redo files')
        else:
            with open(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/hf_ts.xyz"))}',
                      'w+', encoding='utf-8') as xyz:
                xyz.write(str(system_size) + '\n\n')
                for i in coordinates:
                    xyz.write(str(i) + '\n')
            os.remove(f'{str(Path("/")).join(out_file.split(str(Path("/")))[:-1])}{str(Path("/out_rev.txt"))}')


def ResubmitOnAugusta(folder, GausJobs, settings, molecules_dict, resub_compound_and_sites, distance_tweak = False, wholecompound=False):
    """
    Function that executes calculations on the Augusta cluster.
    @param folder: The name of the remote folder where files will be uploaded/downloaded from and where calcualtions are
     executed.
    @param GausJobs: List of calculations to be run on Augusta
    @param settings: Settings Class
    @param molecules_dict: Dictionary with the hashed SMILES string of each compound as the key and the list of sites as
     the value.
    @param resub_compound_and_sites: Dictionary with the hashed molecule SMILES string as the key and the list of sites as
    the value.
    @param wholecompound: Whether the full compound is to be calculated or just a single calculation.
    @return: Error if folder already found on the cluster.
    """
    print('writing slurm scripts for ')
    print(GausJobs)
    fullfolder = home + folder

    JobLogs = InitJobLog('/home/' + settings.user + '/' + folder, GausJobs, hostname, settings.nwchem)

    # Write the slurm scripts
    SubFiles, ScratchFolders, JobLogs, i = WriteAugustaScripts(GausJobs, folder, settings, JobLogs, molecules_dict,
                                                            resub_compound_and_sites, wholecompound, resubmission=True, distance_tweak=distance_tweak)

    WriteJobLogFile(JobLogs)

    print(str(len(SubFiles)) + ' slurm scripts generated')

    # Upload .com files and slurm files to directory
    print("Uploading files to " + hostname + "...")

    remote_filelist = str(Path('/')).join(GausJobs[0].split(str(Path("/")))[:-3]) + str(Path('/remote_filelist%s')) % i
    rem_file = 'remote_filelist%s' % i
    with open(remote_filelist, 'w', encoding='utf-8') as rem:
        for jo in GausJobs:
            rem.write('/'.join(['/'.join(jo.split(str(Path("/")))[-3:-1]),
                                            jo.split(str(Path("/")))[-1].split('.')[0]]) + '.out' + '\n')
    print('written remote_filelist')
    filelist = str(Path('/')).join(GausJobs[0].split(str(Path("/")))[:-3]) + str(Path('/filelist%s')) % i
    with open(filelist, 'w', encoding='utf-8') as fil:
        fil.write(f'{rem_file}\n')
        for job in GausJobs:
            print('/'.join(job.split(str(Path("/")))[-3:]) + '\n')
            fil.write('/'.join(job.split(str(Path("/")))[-3:]) + '\n')
        fil.write('\n')
        for file in SubFiles:
            print('/'.join(file.split(str(Path("/")))[-1:]) + '\n')
            fil.write('/'.join(file.split(str(Path("/")))[-1:]) + '\n')
        fil.write('\n')
    print('written filelist')
    remote_location = f'{hostname}:{fullfolder}'
    print(remote_location)
    molecule_location = str(Path("/")).join(GausJobs[0].split(str(Path("/")))[1:-3])
    print(molecule_location)
    starting_loc = str(Path("/")).join(GausJobs[0].split(str(Path("/")))[6:-3])
    print(starting_loc)
    # outp = subprocess.Popen(
    #     ['cd', str(Path('/')).join(GausJobs[0].split("\\")[:-2]), ';rsync', '-a', '--files-from=filelist', '.',
    #      f'{hostname}:{fullfolder}/'], stderr=subprocess.STDOUT,
    #     stdout=subprocess.PIPE).communicate()[0]

    # subprocess.run(
    #     ['rsync', '-a', '--files-from=' + molecule_location + '\\filelist', molecule_location + '\\',
    #      remote_location], stderr=subprocess.STDOUT,
    #     stdout=subprocess.PIPE)
    print('rsync', '-a', f'--files-from={str(Path("/"))}{molecule_location}{str(Path("/"))}filelist%s' % i,
          f'./{starting_loc} {remote_location}')
    if not sulis:
        p = subprocess.Popen(
            ['rsync', '-a', f'--files-from={str(Path("/"))}{molecule_location}{str(Path("/"))}filelist%s' % i,
             f'./{starting_loc}', remote_location],
            stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    else:
        p = subprocess.Popen(
            ['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a',
             f'--files-from={str(Path("/"))}{molecule_location}{str(Path("/"))}filelist%s' % i,
             f'./{starting_loc}', remote_location],
            stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0:
        print('ERROR' + str((p.returncode, output, error)))

    # status = outp.decode().split('\n')[-2]
    # print(status + ' IS GOING WRONG')
    print(str(len(GausJobs)) + ' .com and ' + str(len(SubFiles)) + \
          ' slurm files uploaded to ' + hostname)

    JobIDs = []
    batchfiles = []
    for f in SubFiles:
        if 'nwchem' in f.split(str(Path("/")))[-1]:
            batchfiles.append(f)
    # Launch the calculations
    for b in batchfiles:
        if not sulis:
            print('ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';dos2unix *')
            pre = subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';dos2unix *'],
                                   stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            print('ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch', b.split(str(Path("/")))[-1])
            outp = \
                subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch',
                                  b.split(str(Path("/")))[-1]],
                                 stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        else:
            print('ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}',
                  'cd ' + fullfolder + '/' + ';dos2unix *')
            pre = subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}',
                                    'cd ' + fullfolder + '/' + ';dos2unix *'],
                                   stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            print('ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ' + fullfolder + '/' + ';sbatch',
                  b.split(str(Path("/")))[-1])
            outp = \
                subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}',
                                  'cd ' + fullfolder + '/' + ';sbatch',
                                  b.split(str(Path("/")))[-1]],
                                 stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

        status = outp.decode().split('\n')[-2]
        JobIDs.append(status.split('job ')[1])

    print(str(len(SubFiles)) + ' jobs submitted to the queue on ' + hostname)
    if settings.SubmitOnly:
        print('Submission only requested, exiting after successful submission of jobs to ' + hostname)
        quit()

    job_slurm_dict = {j: '' for j in GausJobs}

    OldQRes = CheckSlurmQueue(JobIDs)
    JobIDs = OldQRes[3]
    return OldQRes


def write_int_xyz(files_list, folder, files_time_dict, nwchem, job_slurm_dict, settings, wholecompound=False):
    """Function that writes the .xyz files for intermediate steps in the transition state searches to separate files.

    @param files_list: List of files to be scrubbed and grab the geometries from.
    @param folder: Folder to write the files to.
    @param files_time_dict: Dictionary with the sites as the key and the latest time as the value, prevents writing the
     same geometry twice if optimisation is slow.
    @param nwchem: Whether the NWChem software package is being used or not.
    @param job_slurm_dict: Dictionary with the local input file path as the key and the JobID as the value.
    @param settings: Settings class
    @param wholecompound: Whether the wholemolecule is to be calculated or just a single site.
    @return: the files_time dict for use when this function called again."""
    jobs_to_restart = []
    jobs_to_delete = []
    resub_compounds_and_sites = {}
    if not files_time_dict:
        for file in files_list:
            site = file.split('/')[0]
            files_time_dict[site] = []
    for file in files_list:
        site = file.split('/')[0]
        b = []
        opt_no = 1
        print(file)
        print('/'.join(file.split(str(Path('/')))[-3:-1]))
        print(file.split(str(Path('/')))[-1].split('.')[0])
        if wholecompound is True:
            direc = file.split(str(Path('/')))[0] + '/'
        else:
            direc = ''
        atom_dict = {'1': 'H', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl',
                     '35': 'Br',
                     '53': 'I'}
        remote_site = folder + '/' + '/'.join(file.split(str(Path('/')))[-3:-1])
        hfoutp = file.split(str(Path('/')))[-1].split('.')[0]
        if not sulis:
            print('rsync', '-a', '-v', hostname + ':' + folder + '/' + '/'.join(
                ['/'.join(file.split(str(Path('/')))[-3:-1]), file.split(str(Path('/')))[-1].split('.')[0]]) + '.out',
                                     str(Path('/')).join(['.', str(Path('/')).join(file.split(str(Path('/')))[6:-1]),
                                                          file.split(str(Path('/')))[-1].split('.')[0]]) + '.out')
            outp = subprocess.Popen(['rsync', '-a', '-v', hostname + ':' + folder + '/' + '/'.join(
                ['/'.join(file.split(str(Path('/')))[-3:-1]), file.split(str(Path('/')))[-1].split('.')[0]]) + '.out',
                                     str(Path('/')).join(['.', str(Path('/')).join(file.split(str(Path('/')))[6:-1]),
                                                          file.split(str(Path('/')))[-1].split('.')[0]]) + '.out'],
                                    stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            status = outp.decode().split('\n')[-2]
            print(status + ' IS GOING WRONG')
        else:
            print('rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a', hostname + ':' + folder + '/' + '/'.join(
                ['/'.join(file.split(str(Path('/')))[-3:-1]), file.split(str(Path('/')))[-1].split('.')[0]]) + '.out',
                                     str(Path('/')).join(['.', str(Path('/')).join(file.split(str(Path('/')))[6:-1]),
                                                          file.split(str(Path('/')))[-1].split('.')[0]]) + '.out')
            print('heere')
            outp = subprocess.Popen(['rsync', '-e', r'ssh -i C:\Users\pcypw1\.ssh\sulis_key', '-a', hostname + ':' + folder + '/' + '/'.join(
                ['/'.join(file.split(str(Path('/')))[-3:-1]), file.split(str(Path('/')))[-1].split('.')[0]]) + '.out',
                                     str(Path('/')).join(['.', str(Path('/')).join(file.split(str(Path('/')))[6:-1]),
                                                          file.split(str(Path('/')))[-1].split('.')[0]]) + '.out'],
                                    stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            print('command attempted')
            status = outp.decode()
            print(status + ' IS GOING WRONG')
            sleep(2)
        try:
            if nwchem is False:
                with open(
                        f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.out', encoding='utf-8') as rf:
                    rev_lines = reversed(rf.readlines())
                with open(
                        Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt'), 'w', encoding='utf-8') as wf:
                    print('opening reversed_file')
                    for line in rev_lines:
                        wf.write(line)
                site_dir = str(Path("/")).join(file.split(str(Path("/")))[:-1])
                if os.path.getsize(Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt')) == 0:
                    print('empty')
                    continue
                else:
                    print('not empty')
                with open(Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt'), encoding='utf-8') as fil:
                    for line in fil:
                        lines = line.strip('\n')
                        if 'Leave Link  202' in lines:
                            time = lines[30:40]
                            if time in files_time_dict[site]:
                                continue
                            else:
                                files_time_dict[site].append(time)
                                # collect block-related lines
                                while True:
                                    try:
                                        lines = next(fil)
                                    except StopIteration:
                                        # there is no lines left
                                        break
                                    if 'Number' in lines:
                                        # we've reached the end of block
                                        break
                                    b.append(lines)
                                # stop iterating over file
                                break
                if not b:
                    continue
                else:
                    b.reverse()
                    trim = b[1:-2]
                    coordinates = []

                    # This leaves ATOM LABEL, X, Y Z
                    for atom in trim:
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
                    no_of_atoms = len(coordinates)
                    if 'C' in coordinates[-4]:
                        if 'H' in coordinates[-2]:
                            radical = 'cf2h'
                        else:
                            radical = 'cf3'
                    else:
                        radical = 'ipr'
                    while os.path.isfile(Path(f'{site_dir}/geo-opt_{opt_no}.xyz')) is True:
                        opt_no += 1
                    with open(Path(f'{site_dir}/geo-opt_{opt_no}.xyz'), 'w+', encoding='utf-8') as xyz:
                        xyz.write(f'{no_of_atoms}\n\n')
                        for i in coordinates:
                            xyz.write(str(i) + '\n')
                distance_check = hf_distance_checks(Path(f'{site_dir}/geo-opt_{opt_no}.xyz'), radical,
                                                    site.split(str(Path('/')))[-2])
                if distance_check == 3:
                    print('distance okay')
                    continue
                else:
                    print('distance not fine')
                    print(job_slurm_dict)
                    jobs_to_delete.append(job_slurm_dict[str(Path(f'{site_dir}/hfoutp.in'))])
                    print(jobs_to_delete)
                    print('ssh', f'{ssh_key}', f'{hostname}', 'scancel',
                          job_slurm_dict[str(Path(f'{site_dir}/hfoutp.in'))])
                    outp2 = subprocess.Popen(
                        ['ssh', f'{ssh_key}', f'{hostname}', 'scancel',
                         job_slurm_dict[str(Path(f'{site_dir}/hfoutp.in'))]],
                        stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    print('ssh', f'{ssh_key}', f'{hostname}', 'cd ', str(remote_site), ';rm ',
                          str(hfoutp) + '*')
                    outp = subprocess.Popen(
                        ['ssh', f'{ssh_key}', f'{hostname}', 'cd ', str(remote_site), ';rm ',
                         str(hfoutp) + '*'],
                        stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    fileList = glob(f'{site_dir}/geo-opt_*.xyz')
                    for filePath in fileList:
                        try:
                            os.remove(filePath)
                        except:
                            print("Error while deleting file : ", filePath)
                    try:
                        os.remove(Path(f'{site_dir}/hf2outp.out'))
                        os.remove(Path(f'{site_dir}/hf2.sdf'))
                        os.remove(Path(f'{site_dir}/hf2.xyz'))
                        os.remove(Path(f'{site_dir}/hf2outp.in'))
                        with open(Path(f'{site_dir}/hf2.dat'), 'a', encoding='utf-8') as e:
                            lin = e.readlines()
                            if radical == 'cf3':
                                li = lin[-4].split()
                            elif radical == 'cf2h':
                                li = lin[-4].split()
                            elif radical == 'ipr':
                                li = lin[-10].split()
                            li[1] = '2.080000\t'
                            t = " "
                            t = t.join(li)
                            if radical == 'cf3' or radical == 'cf2h':
                                lin[-4] = t
                            else:
                                lin[-10] = t
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf2.dat'), '-oxyz',
                                    '-O', Path(f'{site_dir}/hf2.xyz')], stdout=DEVNULL, stderr=STDOUT)
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf2.dat'), '-osdf',
                                    '-O', Path(f'{site_dir}/hf2.sdf')], stdout=DEVNULL, stderr=STDOUT)
                        jobs_to_restart.append(str(Path(f'{site_dir}/hf2.sdf')))

                    except FileNotFoundError:
                        os.remove(Path(f'{site_dir}/hfoutp.out'))
                        os.remove(Path(f'{site_dir}/hf.sdf'))
                        os.remove(Path(f'{site_dir}/hf.xyz'))
                        os.remove(Path(f'{site_dir}/hfoutp.in'))
                        with open(Path(f'{site_dir}/hf.dat'), 'r+', encoding='utf-8') as e:
                            lin = e.readlines()
                            if radical == 'cf3':
                                li = lin[-4].split()
                            elif radical == 'cf2h':
                                li = lin[-4].split()
                            elif radical == 'ipr':
                                li = lin[-10].split()
                            if distance_check == 1:
                                li[1] = '2.080000\t'
                            elif distance_check == 2:
                                li[1] = '2.160000\t'
                            t = " "
                            t = t.join(li)
                            if radical == 'cf3' or radical == 'cf2h':
                                lin[-4] = t
                            else:
                                lin[-10] = t
                        os.remove(Path(f'{site_dir}/hf.dat'))
                        with open(Path(f'{site_dir}/hf.dat'), 'w+', encoding='utf-8') as f:
                            for line in lin:
                                f.write(str(line))
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf.dat'), '-oxyz',
                                    '-O', Path(f'{site_dir}/hf.xyz')], stdout=DEVNULL, stderr=STDOUT)
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf.dat'), '-osdf',
                                    '-O', Path(f'{site_dir}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
                        jobs_to_restart.append(str(Path(f'{site_dir}/hf.sdf')))
                        continue
            else:
                if os.path.isfile(
                        Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt')):
                    os.remove(
                        Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt'))
                print(f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.out')
                with open(
                        f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.out', encoding='utf-8') as rf:
                    print('opening output file')
                    rev_lines = reversed(rf.readlines())
                    print(rev_lines)
                with open(
                        Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt'), 'w', encoding='utf-8') as wf:
                    print('opening reversed_file')
                    for line in rev_lines:
                        wf.write(line)
                    print('file written')
                    site_dir = str(Path("/")).join(file.split(str(Path("/")))[:-1])
                    print(site_dir)
                if os.path.getsize(Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt')) == 0:
                    print('empty')
                    continue
                else:
                    print('not empty')
                with open(Path(f'{str(Path("/")).join(file.split(str(Path("/")))[:-1])}/out_rev.txt'), encoding='utf-8') as fil:
                    for line in fil:
                        lines = line.strip('\n')
                        if 'Atomic Mass' in lines:
                            print('collecting coordinates')
                            # collect block-related lines
                            while True:
                                try:
                                    lines = next(fil)
                                except StopIteration:
                                    # there is no lines left
                                    break
                                if 'Output' in lines:
                                    # we've reached the end of block
                                    break
                                b.append(lines)
                            # stop iterating over file
                            break
                    print(f'HERE {b}')
                    if not b:
                        continue
                    else:
                        b.reverse()
                        trim = b[3:-1]
                        coordinates = []
                        for atom in trim:
                            strip = atom.strip('\n')
                            line_elements = re.split(' +', str(strip))
                            element = line_elements[2]
                            coord = line_elements[4:]
                            co = [element]
                            for coor in coord:
                                co.append(coor)
                            coordinate = '\t'.join(map(str, co))
                            coordinates.append(coordinate)
                        no_of_atoms = len(coordinates)
                        if 'C' in coordinates[-4]:
                            if 'H' in coordinates[-2]:
                                radical = 'cf2h'
                            else:
                                radical = 'cf3'
                        else:
                            radical = 'ipr'
                        while os.path.isfile(Path(f'{site_dir}/geo-opt_{opt_no}.xyz')) is True:
                            opt_no += 1
                        with open(Path(f'{site_dir}/geo-opt_{opt_no}.xyz'), 'w+', encoding='utf-8') as xyz:
                            xyz.write(f'{no_of_atoms}\n\n')
                            for i in coordinates:
                                xyz.write(str(i) + '\n')
                distance_check = hf_distance_checks(Path(f'{site_dir}/geo-opt_{opt_no}.xyz'), radical, site.split(str(Path('/')))[-2])
                if distance_check == 3:
                    print('distance okay')
                    continue
                else:
                    print('distance not okay')
                    jobs_to_delete.append(job_slurm_dict[f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.in'])
                    if not sulis:
                        print('ssh', f'{ssh_key}', f'{hostname}', 'scancel',
                                 job_slurm_dict[f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.in'])
                        outp2 = subprocess.Popen(
                                ['ssh', f'{ssh_key}', f'{hostname}', 'scancel',
                                 job_slurm_dict[f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.in']],
                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                        print('ssh', f'{ssh_key}', f'{hostname}', 'cd ', str(remote_site), ';rm ',
                             str(hfoutp) + '*')
                        outp = subprocess.Popen(
                            ['ssh', f'{ssh_key}', f'{hostname}', 'cd ', str(remote_site), ';rm ',
                             str(hfoutp) + '*'],
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    else:
                        print('ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}',
                              'scancel',
                                 job_slurm_dict[f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.in'])
                        pre = subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}',
                              'scancel',
                                 job_slurm_dict[f'{str(Path("/")).join([str(Path("/")).join(file.split(str(Path("/")))[:-1]), file.split(str(Path("/")))[-1].split(".")[0]])}.in']],
                                               stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                        print('ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ', str(remote_site), ';rm ',
                              str(hfoutp) + '*')
                        outp = subprocess.Popen(
                            ['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", f'{hostname}', 'cd ', str(remote_site), ';rm ',
                             str(hfoutp) + '*'],
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    fileList = glob(f'{site_dir}/geo-opt_*.xyz')
                    for filePath in fileList:
                        try:
                            os.remove(filePath)
                        except:
                            print("Error while deleting file : ", filePath)
                    try:
                        with open(Path(f'{site_dir}/hf2.dat'), 'r+', encoding='utf-8') as e:
                            lin = e.readlines()
                            if radical == 'cf3':
                                li = lin[-4].split()
                            elif radical == 'cf2h':
                                li = lin[-4].split()
                            elif radical == 'ipr':
                                li = lin[-10].split()
                            if distance_check == 1:
                                print('DISTANCE IS ')
                                print(float(li[1].strip('\t')))
                                if float(li[1].strip('\t')) >= 2.120000:
                                    li[1] = '2.080000\t'
                                elif float(li[1].strip('\t')) <= 2.119999:
                                    print('SITE STILL FLYING APART')
                                    continue
                            elif distance_check == 2:
                                print('DISTANCE IS ')
                                print(float(li[1].strip('\t')))
                                if float(li[1].strip('\t')) <= 2.120000:
                                    li[1] = '2.160000\t'
                                elif float(li[1].strip('\t')) >= 2.160000:
                                    print('SITE STILL SNAPPING TOGETHER')
                                    continue
                            if li[-1] != '\n':
                                li.append('\n')
                            t = " "
                            t = t.join(li)
                            if radical == 'cf3' or radical == 'cf2h':
                                lin[-4] = t
                            else:
                                lin[-10] = t
                        os.remove(Path(f'{site_dir}/hf2outp.out'))
                        os.remove(Path(f'{site_dir}/hf2.sdf'))
                        os.remove(Path(f'{site_dir}/hf2.xyz'))
                        os.remove(Path(f'{site_dir}/hf2outp.in'))
                        os.remove(Path(f'{site_dir}/hf2.dat'))
                        with open(Path(f'{site_dir}/hf2.dat'), 'w+', encoding='utf-8') as f:
                            for line in lin:
                                f.write(str(line))
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf2.dat'), '-oxyz',
                                    '-O', Path(f'{site_dir}/hf2.xyz')], stdout=DEVNULL, stderr=STDOUT)
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf2.dat'), '-osdf',
                                    '-O', Path(f'{site_dir}/hf2.sdf')], stdout=DEVNULL, stderr=STDOUT)
                        jobs_to_restart.append(str(Path(f'{site_dir}/hf2.sdf')))

                    except FileNotFoundError:
                        with open(Path(f'{site_dir}/hf.dat'), 'r+', encoding='utf-8') as e:
                            lin = e.readlines()
                            if radical == 'cf3':
                                li = lin[-4].split()
                            elif radical == 'cf2h':
                                li = lin[-4].split()
                            elif radical == 'ipr':
                                li = lin[-10].split()
                            if distance_check == 1:
                                print('DISTANCE IS ')
                                print(float(li[1].strip('\t')))
                                if float(li[1].strip('\t')) >= 2.120000:
                                    li[1] = '2.080000\t'
                                elif float(li[1].strip('\t')) <= 2.119999:
                                    print('SITE STILL FLYING APART')
                                    continue
                            elif distance_check == 2:
                                if float(li[1].strip('\t')) <= 2.120000:
                                    li[1] = '2.160000\t'
                                elif float(li[1].strip('\t')) >= 2.160000:
                                    print('SITE STILL SNAPPING TOGETHER')
                                    continue
                            if li[-1] != '\n':
                                li.append('\n')
                            t = " "
                            t = t.join(li)
                            if radical == 'cf3' or radical == 'cf2h':
                                lin[-4] = t
                            else:
                                lin[-10] = t
                        os.remove(Path(f'{site_dir}/hfoutp.out'))
                        os.remove(Path(f'{site_dir}/hf.sdf'))
                        os.remove(Path(f'{site_dir}/hf.xyz'))
                        os.remove(Path(f'{site_dir}/hfoutp.in'))
                        os.remove(Path(f'{site_dir}/hf.dat'))
                        with open(Path(f'{site_dir}/hf.dat'), 'w+', encoding='utf-8') as f:
                            for line in lin:
                                f.write(str(line))
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf.dat'), '-oxyz',
                                    '-O', Path(f'{site_dir}/hf.xyz')], stdout=DEVNULL, stderr=STDOUT)
                        check_call(['obabel', '-imopin', Path(f'{site_dir}/hf.dat'), '-osdf',
                                    '-O', Path(f'{site_dir}/hf.sdf')], stdout=DEVNULL, stderr=STDOUT)
                        jobs_to_restart.append(str(Path(f'{site_dir}/hf.sdf')))
                        continue

        except FileNotFoundError:
            continue
    if jobs_to_restart:
        for r in jobs_to_restart:
            if r.split(str(Path('/')))[-3] in resub_compounds_and_sites:
                if r.split(str(Path('/')))[-2] in resub_compounds_and_sites[r.split(str(Path('/')))[-3]]:
                    continue
                else:
                    resub_compounds_and_sites[r.split(str(Path('/')))[-3]] += [r.split(str(Path('/')))[-2]]
            else:
                resub_compounds_and_sites[r.split(str(Path('/')))[-3]] = [r.split(str(Path('/')))[-2]]
    else:
        return files_time_dict, {}, jobs_to_delete
    resub_GausInpFiles, resub_molecules_dict = Gaussian.SetupGaussian(jobs_to_restart, settings,
                                                                      wholemolecule=True
                                                                      , resubmission=False)
    QRun = False
    Files2Run = Gaussian.GetFiles2Run(resub_GausInpFiles, settings)
    print(f'FILES TO RUN:{Files2Run}')
    if len(Files2Run) == 0:
        QRun = True

    if settings.GenOnlyQ:
        print("Gaussian input files generated, quitting as instructed.")
        quit()

    if QRun:
        print('DFT has already been run for these inputs. Skipping...')
    else:
        print('running retry')
        QRes = RunCalcs(Files2Run, settings, resub_molecules_dict, resubmission=True,
                        resub_compound_and_sites=resub_compounds_and_sites, distance_tweak=True)
        return files_time_dict, QRes[4], jobs_to_delete


def IsAugustaGComplete(folder, settings, compounds_and_sites, wholemolecule):
    """
    Function that greps the number of times the Calculation complete keyphrase in each Job is seen.
    @param folder: remote folder on the server where calcualations are performed.
    @param settings: Settings class
    @param compounds_and_sites: Dictionary with compound hashed SMILES as key and list of sites as the value.
    @param wholemolecule: Whether the whole compound is being calculated or just a single calculation.
    @return: The number of jobs that are complete.
    """
    JobsComplete = 0
    for compound in compounds_and_sites.keys():
        path = folder + '/' + compound + '/'
        try:
            # outp = subprocess.Popen(['ssh ', hostname + ' grep Normal ' + path + '*out | wc -l'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            if wholemolecule is False:
                outp = \
                subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}', ' grep Normal ' + path + '*out | wc -l'], stderr=subprocess.STDOUT,
                                 stdout=subprocess.PIPE).communicate()[0]
            else:
                if settings.nwchem is False:
                    outp = subprocess.Popen(['ssh', f'{ssh_key}', f'{hostname}',
                                             ' cd ' + path + ';for f in $(find . -name "hfoutp.out"); do grep Normal $f | tail -1; done | wc -l'],
                                            stderr=subprocess.STDOUT,
                                            stdout=subprocess.PIPE).communicate()[0]
                else:
                    if not sulis:
                        outp = subprocess.Popen(['ssh ', hostname, ' grep "Total times" ' + path + '*/*out | wc -l'],
                                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    else:
                        outp = subprocess.Popen(['ssh', '-i', r"C:\Users\pcypw1\.ssh\sulis_key", hostname, ' grep "Total times" ' + path + '*/*out | wc -l'],
                                                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            JobsComplete += int(outp.decode().split(str(Path('\n')))[-2])

        except ValueError:
            JobsComplete = 0
    print('JOBS COMPLETE ' + str(JobsComplete))
    return JobsComplete


def hf_distance_checks(geometry, carbon_radical, site):
    """
    Function that checks whether the transition state geometry has a C-C bond length within a typical range.
    @param geometry: Cartesian coordinates file of transition state.
    @param carbon_radical: Which Functional group is being investigated.
    @param site: Which site within the compound is being investigated.
    @return: Boolean that decides whether the distance is within the typical range.
    """
    print('here')
    print(type(site))
    if site == 'reagent':
        return 3
    elif site == 'fukui':
        return 3
    else:
        site = int(site)
        with open(f'{geometry}', 'r+', encoding='utf-8')as f:
            print('opened file')
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
            print(site)
            if carbon_radical == 'cf3':
                carbon = cart[-4]
            elif carbon_radical == 'cf2h':
                carbon = cart[-4]
            elif carbon_radical == 'ipr':
                carbon = cart[-10]
            site_of_interest = cart[site-1]
            print('comparing atoms')
            [a, x1, y1, z1] = carbon
            [b, x2, y2, z2] = site_of_interest
            dist = float((((float(x2) - float(x1)) ** 2) + ((float(y2) - float(y1)) ** 2) + (
                        (float(z2) - float(z1)) ** 2)) ** (1 / 2))
            print(dist)

            if dist >= 2.80:
                print('DISTANCE TOO FAR')
                return 1
            elif dist <= 1.8:
                print('DISTANCE TOO CLOSE')
                return 2
            else:
                print('DISTANCE GOOD')
                return 3


