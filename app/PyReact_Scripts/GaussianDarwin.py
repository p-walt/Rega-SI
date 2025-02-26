#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Jul 25 2021

@author: K. Ermanis

Contains all of the code specific to running Gaussian on Darwin/CSD3 in Cambridge.
Calls a lot of code from Gaussian.py  for input generation and calculation
execution. Called by PyReact.py.
"""

import Gaussian

import subprocess
import os
import time
import shutil
import math
import datetime
from openbabel import *

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

def RunCalcs(Files2Run, settings):

    print('\nRunning Gaussian on Darwin...')
    # Run Gaussian jobs on Darwin cluster in folder named after date
    # and title and wait until the last file is completed
    now = datetime.datetime.now()
    MaxCon = settings.MaxConcurrentJobsDarwin

    if len(Files2Run) < MaxCon:
        settings.folders.append(now.strftime('%d%b%H%M') + settings.Title)
        RunOnDarwin(now.strftime('%d%b%H%M') + settings.Title,
                             Files2Run, settings)
    else:
        print("The DFT calculations will be done in " + \
              str(math.ceil(len(Files2Run) / MaxCon)) + " batches")
        i = 0
        while (i + 1) * MaxCon < len(Files2Run):
            print("Starting batch nr " + str(i + 1))
            RunOnDarwin(now.strftime('%d%b%H%M') + str(i + 1) +
                                 settings.Title, Files2Run[(i * MaxCon):((i + 1) * MaxCon)],
                                 settings)
            i += 1
        print("Starting batch nr " + str(i + 1))
        RunOnDarwin(now.strftime('%d%b%H%M') + str(i + 1) +
                             settings.Title, Files2Run[(i * MaxCon):], settings)


def RunOnDarwin(folder, GausJobs, settings):
    
    print("Darwin GAUSSIAN job submission script\n")

    #Check that folder does not exist, create job folder on ziggy
    outp = subprocess.Popen(['ssh', 'darwin', 'ls'], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print(folder)
    
    if folder in outp.decode():
        print("Folder exists on Darwin, choose another folder name.")
        return

    outp = subprocess.Popen(['ssh', 'darwin', 'mkdir', folder], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    JobLogs = InitJobLog('/home/' + settings.user + '/' + folder, GausJobs, 'darwin')

    #Write the slurm scripts
    SubFiles, ScratchFolders, JobLogs = WriteDarwinScripts(GausJobs, settings, JobLogs)

    WriteJobLogFile(JobLogs)

    print(str(len(SubFiles)) + ' slurm scripts generated')

    for scrf in ScratchFolders:
        fullscrf = '/home/' + settings.user + settings.t2scrf + settings.user + '/' + scrf + '\n'
        outp = subprocess.Popen(['ssh', 'darwin', 'mkdir', fullscrf], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print(str(len(ScratchFolders)) + ' scratch folders created')

    #Upload .com files and slurm files to directory
    print("Uploading files to darwin...")
    for f in GausJobs:
        outp = subprocess.Popen(['scp', f,
        'darwin:/home/' + settings.user + '/' + folder],
        stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

        
    for f in SubFiles:
        
        outp = subprocess.Popen(['scp', f,
            'darwin:/home/' + settings.user + '/' + folder], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    print(str(len(GausJobs)) + ' .com and ' + str(len(SubFiles)) +\
        ' slurm files uploaded to darwin')

    if settings.GuessChk:
        for f in GausJobs:
            outp = subprocess.Popen(['scp', f[:-4] + '.chk',
            'darwin:/home/' + settings.user + '/' + folder],
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print(str(len(GausJobs)) + ' .chk files uploaded to darwin')

    fullfolder = '/home/' + settings.user + '/' + folder

    if (settings.KeepChk == True) and (settings.ResubsDone > 0):
        print('Copying .chk files to the new folder')
        for f in GausJobs:
            for c in settings.ChkFiles:
                if f[:-4] == c[0]:
                    print('ssh darwin cp /home/' + settings.user + '/' + settings.folders[-2] \
                                             + '/' + c[1] + '.chk ' + fullfolder + '/')
                    outp = subprocess.Popen(['ssh', 'darwin', 'cp ' + '/home/' + settings.user + '/' \
                                             + settings.folders[-2] + '/' + c[1] + '.chk ' + fullfolder + '/',], \
                                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
                    print('ssh darwin rm /home/' + settings.user + '/' + settings.folders[-2] \
                          + '/' + c[1] + '.chk')
                    outp = subprocess.Popen(['ssh', 'darwin', 'rm ' + '/home/' + settings.user + '/' \
                                             + settings.folders[-2] + '/' + c[1] + '.chk', ], \
                                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print(str(len(GausJobs)) + ' .chk files copied from previous folder on darwin, original files deleted.')

    JobIDs = []
    #Launch the calculations
    for f in SubFiles:
        outp = subprocess.Popen(['ssh', 'darwin', 'cd ' + fullfolder + ';sbatch', f], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        status = outp.decode().split('\n')[-2]
        print(status)
        JobIDs.append(status.split('job ')[1])

    print(str(len(SubFiles)) + ' jobs submitted to the queue on darwin ' + \
        'containing ' + str(len(GausJobs)) + ' Gaussian jobs')

    time.sleep(60)

    OldQRes = CheckDarwinQueue(JobIDs, settings)

    while OldQRes[0] < 0:
        OldQRes = CheckDarwinQueue(JobIDs, settings)
        time.sleep(60)

    print('Pending: ' + str(OldQRes[0]) + ', Running: ' + str(OldQRes[1]) + ', Not in queue: ' + str(OldQRes[2]))

    Jobs2Complete = list(GausJobs)
    n2complete = len(Jobs2Complete)
    
    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JobFinished = IsDarwinGComplete(Jobs2Complete, folder, settings)
        
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not JobFinished[job[:-3] + 'out']]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print(str(n2complete) + " remaining.")

        QRes = CheckDarwinQueue(JobIDs, settings)
        if QRes != OldQRes:
            if QRes[0] < 0:
                QRes = OldQRes
            else:
                OldQRes = QRes
                print('Pending: ' + str(OldQRes[0]) + ', Running: ' + str(OldQRes[1]) + ', Not in queue: ' + str(OldQRes[2]))


        if (QRes[2] == len(JobIDs)) and (QRes[0] >= 0):
            #check each gaussian file to ascertain the status of individual gaus jobs
            print('No jobs left in Darwin queue')
            break

        time.sleep(180)

    #When done, copy the results back
    print("\nCopying the output files back to localhost...")
    print('scp darwin:' + fullfolder + '/*.out ' + os.getcwd() + '/')

    outp = subprocess.Popen(['scp', 'darwin:' + fullfolder + '/*.out',
            os.getcwd() + '/'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    if settings.KeepChk == True:
        print("\nCopying the checkpoint files back to localhost...")
        print('scp darwin:' + fullfolder + '/*.chk ' + os.getcwd() + '/')

        outp = subprocess.Popen(['scp', 'darwin:' + fullfolder + '/*.chk',
                                 os.getcwd() + '/'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    else:
        print("\nDeleting checkpoint files...")
        print('ssh darwin rm ' + fullfolder + '/*.chk')
        outp = subprocess.Popen(['ssh', 'darwin', 'rm', fullfolder + '/*.chk'], \
                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    if settings.WFN == True:
        print("\nCopying the wfn files back to localhost...")
        print('scp darwin:' + fullfolder + '/*.wfn ' + os.getcwd() + '/')

        outp = subprocess.Popen(['scp', 'darwin:' + fullfolder + '/*.wfn',
                                 os.getcwd() + '/'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]


def CheckDarwinQueue(JobIDs, settings):

    outp = subprocess.Popen(['ssh', 'darwin', 'squeue', '-u ' + settings.user], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    outp = outp.decode().split('\n')
    QStart = -1000
    for i, line in enumerate(outp):
        if 'JOBID' in line:
            QStart = i+1
            break

    if QStart < 0:
        return -100, -100, -100

    QueueReport = outp[QStart:-1]
    JobStats = []

    for job in JobIDs:
        status = ''
        for i, line in enumerate(QueueReport):
            if job in line:
                status = [_f for _f in line.split(' ') if _f][4]
        JobStats.append(status)

    Pending = JobStats.count('PD')
    Running = JobStats.count('R')
    NotInQueue = JobStats.count('')

    return Pending, Running, NotInQueue


def ReplaceLine(File, LineN, Line):
    gausf = open(File, 'r+')
    gauslines = gausf.readlines()
    gauslines[LineN] = Line
    gausf.truncate(0)
    gausf.seek(0)
    gausf.writelines(gauslines)
    gausf.close()


def WriteDarwinScripts(GausJobs, settings, JobLogs):
    
    SubFiles = []
    ScratchFolders = []
    AdjNodeSize = int(math.floor(settings.DarwinNodeSize/settings.nProc))
    
    if len(GausJobs) == AdjNodeSize:
        SubFiles.append(WriteSlurm(GausJobs, settings))
        ScratchFolders.append(GausJobs[0])
        for i, GausJob in enumerate(GausJobs):
            JobLogs[i].ScratchFolder = '/home/' + settings.user + settings.t2scrf + settings.user + '/'\
                    + GausJobs[0]
        #if settings.nProc == 32:
        #    for j, GausJob in enumerate(GausJobs):
        #        line = '%CPU=' + str(j*settings.nProc) + '-' + str((j+1)*settings.nProc-1) + '\n'
        #        ReplaceLine(GausJob, 0, line)

    elif len(GausJobs) < AdjNodeSize:
        NewNProc = int(math.floor(settings.DarwinNodeSize/len(GausJobs)))
        SubFiles.append(WriteSlurm(GausJobs, settings, nProc=NewNProc))
        ScratchFolders.append(GausJobs[0])
        for LogIndex, GausJob in enumerate(GausJobs):
            JobLogs[LogIndex].ScratchFolder = '/home/' + settings.user + settings.t2scrf + settings.user + '/'\
                    + GausJobs[0]
        print("Jobs don't fill the Darwin node, nproc increased to " + str(NewNProc))
        for j, GausJob in enumerate(GausJobs):
            line = '%nprocshared=' + str(NewNProc) + '\n'
            ReplaceLine(GausJob, 0, line)
        #if NewNProc == 32:
        #    for j, GausJob in enumerate(GausJobs):
        #        line = '%CPU=' + str(j*NewNProc) + '-' + str((j+1)*NewNProc-1) + '\n'
        #        ReplaceLine(GausJob, 0, line)

    else:
        print("The Gaussian calculations will be submitted as " +
                    str(math.ceil(len(GausJobs)/AdjNodeSize)) +
                    " jobs")
        i = 0
        while (i+1)*AdjNodeSize < len(GausJobs):
            PartGausJobs = list(GausJobs[(i*AdjNodeSize):((i+1)*AdjNodeSize)])
            print("Writing script nr " + str(i+1))
            SubFiles.append(WriteSlurm(PartGausJobs, settings, str(i+1)))
            ScratchFolders.append(PartGausJobs[0])
            for LogIndex in range((i*AdjNodeSize),((i+1)*AdjNodeSize)):
                JobLogs[LogIndex].ScratchFolder = '/home/' + settings.user + settings.t2scrf + settings.user + '/' \
                                           + PartGausJobs[0]
            #if settings.nProc == 32:
            #    for j, GausJob in enumerate(PartGausJobs):
            #        line = '%CPU=' + str(j * settings.nProc) + '-' + str((j + 1) * settings.nProc - 1) + '\n'
            #S        ReplaceLine(GausJob, 0, line)
            i += 1
        
        PartGausJobs = list(GausJobs[(i*AdjNodeSize):])
        if len(PartGausJobs) < AdjNodeSize:
            NewNProc = int(math.floor(settings.DarwinNodeSize / len(PartGausJobs)))
            print("Jobs don't fill the last Darwin node, nproc increased to " + str(NewNProc))
            print("Writing script nr " + str(i + 1))
            SubFiles.append(WriteSlurm(PartGausJobs, settings, nProc=NewNProc))
            ScratchFolders.append(PartGausJobs[0])
            for LogIndex in range((i*AdjNodeSize),len(GausJobs)):
                JobLogs[LogIndex].ScratchFolder = '/home/' + settings.user + settings.t2scrf + settings.user + '/' \
                                           + PartGausJobs[0]
            for j, GausJob in enumerate(PartGausJobs):
                line = '%nprocshared=' + str(NewNProc) + '\n'
                ReplaceLine(GausJob, 0, line)
            #if NewNProc == 32:
            #    for j, GausJob in enumerate(PartGausJobs):
            #        line = '%CPU=' + str(j * NewNProc) + '-' + str((j + 1) * NewNProc - 1) + '\n'
            #        ReplaceLine(GausJob, 0, line)
        else:
            print("Writing script nr " + str(i+1))
            SubFiles.append(WriteSlurm(PartGausJobs, settings, str(i+1)))
            ScratchFolders.append(PartGausJobs[0])
            for LogIndex in range((i*AdjNodeSize),len(GausJobs)):
                JobLogs[LogIndex].ScratchFolder = '/home/' + settings.user + settings.t2scrf + settings.user + '/' \
                                           + PartGausJobs[0]
            #if settings.nProc == 32:
            #    for j, GausJob in enumerate(PartGausJobs):
            #        line = '%CPU=' + str(j * settings.nProc) + '-' + str((j + 1) * settings.nProc - 1) + '\n'
            #        ReplaceLine(GausJob, 0, line)

    return SubFiles, ScratchFolders, JobLogs


def WriteSlurm(GausJobs, settings, index='', nProc=-1):

    if nProc == -1:
        nProc = settings.nProc

    cwd = os.getcwd()
    filename = settings.Title + 'slurm' + index
    
    shutil.copyfile(settings.ScriptDir + '/Defaultslurm',
                    cwd + '/' + filename)
    slurmf = open(filename, 'r+')
    slurm = slurmf.readlines()
    slurm[12] = '#SBATCH -J ' + settings.Title + '\n'
    slurm[14] = '#SBATCH -A ' + settings.project + '\n'
    slurm[19] = '#SBATCH --ntasks=' + str(len(GausJobs)*nProc) + '\n'
    slurm[21] = '#SBATCH --time=' + format(settings.TimeLimit,"02") +\
        ':00:00\n'
    #slurm[57] = 'module load gaussian/16' + '\n'
    #slurm[60] = 'application="g16"' + '\n'
    #if settings.project == settings.t2project:
    #    slurm[62] = 'export GAUSS_SCRDIR=/home/' + settings.user + settings.t2scrf + settings.user + '/'\
    #                + GausJobs[0] + '\n'
    slurm[62] = 'export GAUSS_SCRDIR=/home/' + settings.user + settings.t2scrf + settings.user + '/'\
                    + GausJobs[0] + '\n'

    for f in GausJobs:
        slurm.append('srun --exclusive -n1 -c' + str(nProc) + ' $application < ' + f[:-3] + \
            'com > ' + f[:-3] + 'out 2> error &\n')
    slurm.append('wait\n')

    slurmf.truncate(0)
    slurmf.seek(0)
    slurmf.writelines(slurm)
    
    return filename


def IsDarwinGComplete(GausJobs, folder, settings):

    path = '/home/' + settings.user + '/' + folder + '/'
    results = {}
    
    for f in GausJobs:
        outp = subprocess.Popen(['ssh', 'darwin', 'cat', path + f[:-3] + 'out'],
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].decode()
        #outp = subprocess.check_output('ssh darwin cat ' + path + f,
        #                                    shell=True)
        if ("Normal termination" in outp[-90:]) or ('termination' in outp[-300:] and 'l9999.exe' in outp[-300:]):
            results[f[:-3] + 'out'] = True
        else:
            results[f[:-3] + 'out'] = False

    return results


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
