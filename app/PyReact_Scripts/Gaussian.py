#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:56:54 2017

@author: K. Ermanis

Contains all of the Gaussian specific code for input generation and calculation
execution. Called by PyReact.py.
"""

import subprocess
import os
import time
import shutil
import datetime
from openbabel import *
from pathlib import Path

start = Path.cwd()


class JobLogEntry:
    def __init__(self):
        self.LaunchTime = 0
        self.Host = ''
        self.InputFile = ''
        self.OutputFile = ''
        self.CheckpointFile = ''
        self.ReadWriteFile = ''
        self.ResultFolder = ''
        self.ScratchFolder = ''
        self.Completed = False
        self.CompletionTime = 0


def SetupGaussian(inpfiles, settings, wholemolecule=False, resubmission=False):
    GausFiles2Run = []
    ChkFiles = []
    Geoms = []
    Atomlists = []
    # If we are actually dealing with Gaussian output file as a geometry source
    if settings.nwchem is False:
        nwchem = False
    else:
        nwchem = True

    molecules_dict = {}
    if settings.Rega != '':
        for File in inpfiles:
            if File.split(str(Path("/")))[-3] in molecules_dict:
                continue
            else:
                molecules_dict[File.split(str(Path("/")))[-3]] = 0
                continue
    for f in inpfiles:
        # Gaussian input file
        if f[-4:] == '.out':
            (atoms, geom, charge) = ReadGeometry(f, settings)
            molecules_dict = (molecules_dict, f, atoms)
            charge = max(charge, settings.Charge)
            settings.Charge = max(charge, settings.Charge)
            ChkFiles.append([f[:-4] + '.chk'])
            Geoms.append(geom)
            Atomlists.append(atoms)
        # SDF input file
        elif f[-4:] == '.sdf':
            (atoms, geom, charge) = ReadSDFGeometry(f)
            molecules_dict = core_allocation(molecules_dict, f, atoms)
            ChkFiles.append([f[:-4] + '.chk'])
            Geoms.append(geom)
            Atomlists.append(atoms)
        else:
            print('Unrecognized file type for file ' + f)
            print('Quitting...')

    if settings.Restart:
        # Read previous run information from the job log file first for restarting jobs
        OldJobLogs = ReadJobLogFile()
    print(f'MOLECULES_DICT {molecules_dict}')
    for num, (f, geom, atomlist) in enumerate(zip(inpfiles, Geoms, Atomlists)):
        if settings.Scan == '':
            filename = f[:-4] + settings.suffix
            if not settings.Restart:
                if filename.split(str(Path('/')))[-1] == 'hfoutp':
                    WriteGausFileOpt(filename, geom, atomlist, charge, settings, molecules_dict, wholemolecule, nwchem,
                                     resubmission=False)
                elif filename.split(str(Path('/')))[-1] == 'hf2outp':
                    WriteGausFileOpt(filename, geom, atomlist, charge, settings, molecules_dict, wholemolecule, nwchem,
                                     resubmission=True)
            else:
                # Read previous run information from the job log file first
                WriteGausFileRestart(filename, settings, OldJobLogs, molecules_dict, wholemolecule)
            if settings.GuessChk:
                cwd = start
                shutil.copyfile(cwd / str(ChkFiles[num]), cwd / Path(f'{filename}.chk'))
            if nwchem is False:
                GausFiles2Run.append(filename + '.com')
            else:
                GausFiles2Run.append(filename + '.in')
        else:
            temp = settings.Scan.split(';')
            scancoord, scanvals = temp[0], temp[1]
            for i, val in enumerate(scanvals.split(',')):
                filename = f[:-4] + settings.suffix + str(i + 1)
                WriteGausFileOpt(filename, geom, atomlist, charge, settings, molecules_dict, wholemolecule, nwchem,
                                 resubmission, coord=scancoord, val=val)
                if settings.GuessChk:
                    cwd = start
                    shutil.copyfile(cwd / str(ChkFiles[num]),
                                    cwd / Path(f'{filename}.chk'))
                if nwchem is False:
                    GausFiles2Run.append(filename + '.com')
                else:
                    GausFiles2Run.append(filename + '.in')

    print(str(len(GausFiles2Run)) + " .com files written")

    return GausFiles2Run, molecules_dict


def InitJobLog(folder, GausJobs, host, nwchem):
    JobLogs = []
    now = datetime.datetime.now()
    for GausJob in GausJobs:
        entry = JobLogEntry()
        entry.LaunchTime = now
        entry.Host = host
        if nwchem is False:
            entry.InputFile = GausJob[:-3] + 'com'
        else:
            entry.InputFile = GausJob[:-3] + 'in'
        entry.OutputFile = GausJob[:-3] + 'out'
        entry.CheckpointFile = GausJob[:-3] + 'chk'
        entry.ReadWriteFile = GausJob[:-3] + 'rwf'
        entry.ResultFolder = folder
        JobLogs.append(entry)

    return JobLogs


def WriteJobLogFile(JobLogs):
    print('got here')
    logfile = open('jobs.log', 'a')
    for JobLog in JobLogs:
        for attr, value in list(JobLog.__dict__.items()):
            logfile.write(str(attr) + '=' + str(value) + '\n')
        logfile.write('\n')
    logfile.close()


def ReadJobLogFile():
    JobLogs = []
    jobsf = open('jobs.log', 'r')
    jobslines = jobsf.readlines()
    jobsf.close()
    lineidx = 0
    for x in range(jobslines.count('\n')):
        JobLogs.append(JobLogEntry())
        print(('Job log ' + str(x) + ':'))
        while jobslines[lineidx] != '\n':
            tmp = jobslines[lineidx][:-1].split('=')
            if hasattr(JobLogs[-1], tmp[0]):
                setattr(JobLogs[-1], tmp[0], tmp[1])
                print(('  ' + tmp[0] + ': ' + tmp[1]))
            lineidx += 1
        else:
            lineidx += 1

    return JobLogs


def WriteGausFileRestart(Gausinp, settings, JobLogs, molecules_dict, wholemolecule):
    # write the initial DFT geometry optimisation input file first
        f1 = open(Gausinp + '.com', 'w')
        f1.write('%nprocshared=' + str(settings.nProc) + '\n')

        if (settings.DFT == 'd'):
            f1.write('%mem=' + str(3400 * settings.nProc) + 'MB\n')
        elif (settings.DFT == 'z'):
            if (settings.nProc > 1):
                f1.write('%mem=' + str(1000 * settings.nProc) + 'MB\n')
            else:
                f1.write('%mem=2000MB\n')
        elif (settings.DFT == 'a'):
            if (settings.nProc > 1):
                if settings.Functional == 'HF' or 'UHF':
                    f1.write('%mem=' + str(1000 * settings.nProc) + 'MB\n')
                else:
                    f1.write('%mem=' + str(3800 * settings.nProc) + 'MB\n')
            else:
                f1.write('%mem=2000MB\n')
        else:
            f1.write('%mem=2000MB\n')

        # Find the right joblog data
        RelJobLog = None
        for joblog in JobLogs[::-1]:
            if joblog.InputFile == Gausinp + '.com':
                RelJobLog = joblog
                break
        if RelJobLog == None:
            print('Relevant entry in jobs log not found, restarting not possible, quitting.')
        else:
            print('Relevant entry in jobs log:')
            for attr, value in list(RelJobLog.__dict__.items()):
                print((str(attr) + '=' + str(value)))
        f1.write('%RWF=' + RelJobLog.ScratchFolder + '/' + RelJobLog.ReadWriteFile + '\n%NoSave\n')
        if wholemolecule is False:
            f1.write('%chk=' + Gausinp + '.chk\n')
        else:
            f1.write('%chk=' + Gausinp.split('/')[-1] + '.chk\n')
        f1.write('#P Restart\n\n')
        f1.close()


def WriteGausFileOpt(Gausinp, conformer, atoms, charge, settings, molecules_dict, wholemolecule, nwchem,
                     resubmission=False, coord='', val=''):
    # write the initial DFT geometry optimisation input file first
    if nwchem is False:

        f1 = open(Gausinp + '.com', 'w')
        if settings.Rega != '':
            f1.write('%nprocshared=' + str(molecules_dict[Gausinp.split(str(Path("/")))[-3]]) + '\n')
        else:
            f1.write('%nprocshared=' + str(settings.nProc) + '\n')

        if (settings.DFT == 'd'):
            f1.write('%mem=' + str(3400 * settings.nProc) + 'MB\n')
        elif (settings.DFT == 'z'):
            if (settings.nProc > 1):
                f1.write('%mem=' + str(1000 * settings.nProc) + 'MB\n')
            else:
                f1.write('%mem=2000MB\n')
        elif (settings.DFT == 'a'):
            if (settings.nProc > 1):
                if settings.Functional == 'HF' or 'UHF':
                    f1.write('%mem=' + str(1000 * settings.nProc) + 'MB\n')
                else:
                    f1.write('%mem=' + str(3800 * settings.nProc) + 'MB\n')
            else:
                f1.write('%mem=2000MB\n')
        else:
            f1.write('%mem=2000MB\n')

        if (settings.ResubsDone == 0) or (settings.KeepChk == False):
            if wholemolecule is False:
                f1.write('%chk=' + Gausinp + '.chk\n')
            else:
                f1.write('%chk=' + Gausinp.split(str(Path('/')))[-1] + '.chk\n')
            settings.ChkFiles.append([Gausinp, Gausinp])
        else:
            chk = ''
            for c in settings.ChkFiles:
                if Gausinp == c[0]:
                    chk = c[1]

            if chk != '':
                f1.write('%chk=' + chk + '.chk\n')
            else:
                print(Gausinp + ' does not match any chk entry, exiting...')
                quit()

        BasisFolder = os.path.join(settings.ScriptDir, 'BasisSets')
        CustomBasisSets = os.listdir(BasisFolder)

        try:
            fold = Gausinp.split(str(Path('/')))[-2]
            if fold == 'reagent':
                if settings.BasisSet.lower() in CustomBasisSets:
                    CustomBasisSet = True
                    route = ComposeGausRoute(settings, coord, settings.BasisSet.lower(), reagent=True)
                else:
                    CustomBasisSet = False
                    route = ComposeGausRoute(settings, coord, reagent=True)
            elif fold == 'fukui':
                if settings.BasisSet.lower() in CustomBasisSets:
                    CustomBasisSet = True
                    route = ComposeGausRoute(settings, coord, settings.BasisSet.lower(), reagent=True)
                else:
                    CustomBasisSet = False
                    route = ComposeGausRoute(settings, coord, reagent=True)
            else:
                if settings.BasisSet.lower() in CustomBasisSets:
                    CustomBasisSet = True
                    route = ComposeGausRoute(settings, coord, settings.BasisSet.lower(), reagent=False)
                else:
                    CustomBasisSet = False
                    route = ComposeGausRoute(settings, coord, reagent=False)
        except ValueError:
            if settings.BasisSet.lower() in CustomBasisSets:
                CustomBasisSet = True
                route = ComposeGausRoute(settings, coord, settings.BasisSet.lower(), reagent=False)
            else:
                CustomBasisSet = False
                route = ComposeGausRoute(settings, coord, reagent=False)

        f1.write(route)

        f1.write('\n' + Gausinp + '\n\n')
        ChargeMult = ''

        if settings.ManCharge < -100:
            ChargeMult += str(charge) + ' '
        else:
            ChargeMult += str(settings.ManCharge) + ' '
        if settings.ManMultiplicity < -10:
            ChargeMult += '1\n'
        else:
            ChargeMult += str(settings.ManMultiplicity) + '\n'

        f1.write(ChargeMult)

        if settings.Isotopes != '':
            isomap = InterpretIsotopes(atoms, settings.Isotopes)
        else:
            isomap = ['' for a in atoms]

        natom = 0
        if not settings.ONIOM:
            for atom in conformer:
                if isomap[natom] == '':
                    f1.write(atoms[natom] + '  ' + atom[0] + '  ' + atom[1] + '  ' +
                             atom[2] + '\n')
                else:
                    f1.write(atoms[natom] + '(Iso=' + isomap[natom] + ')  ' + atom[0] + '  ' + atom[1] + '  ' +
                             atom[2] + '\n')
                natom = natom + 1
        else:
            atomsection = GetONIOMatoms(settings, atoms, conformer)

            f1.write(atomsection)

        f1.write('\n')

        # Part for setting up constraints and/or coordinate scan
        if settings.Constraints != '':
            splitconstr = settings.Constraints.split(';')
            for line in splitconstr:
                temp = line.split(',')
                if len(temp) == 3:
                    if temp[-1] != 'F':
                        f1.write('B ' + ' '.join(temp) + ' B\n')
                        f1.write('B ' + ' '.join(temp[:-1]) + ' F\n')
                    else:
                        f1.write('B ' + ' '.join(temp) + '\n')
                elif len(temp) == 4:
                    if temp[-1] != 'F':
                        f1.write('A ' + ' '.join(temp) + ' B\n')
                        f1.write('A ' + ' '.join(temp[:-1]) + ' F\n')
                    else:
                        f1.write('A ' + ' '.join(temp) + '\n')
                elif len(temp) == 5:
                    if temp[-1] != 'F':
                        f1.write('D ' + ' '.join(temp) + ' B\n')
                        f1.write('D ' + ' '.join(temp[:-1]) + ' F\n')
                    else:
                        f1.write('D ' + ' '.join(temp) + '\n')

        if coord != '':
            splitconstr = [coord + ',' + val]
            for line in splitconstr:
                temp = line.split(',')
                if len(temp) == 3:
                    if temp[-1] != 'F':
                        f1.write('B ' + ' '.join(temp) + ' B\n')
                        f1.write('B ' + ' '.join(temp[:-1]) + ' F\n')
                    else:
                        f1.write('B ' + ' '.join(temp) + '\n')
                elif len(temp) == 4:
                    if temp[-1] != 'F':
                        f1.write('A ' + ' '.join(temp) + ' B\n')
                        f1.write('A ' + ' '.join(temp[:-1]) + ' F\n')
                    else:
                        f1.write('A ' + ' '.join(temp) + '\n')
                elif len(temp) == 5:
                    if temp[-1] != 'F':
                        f1.write('D ' + ' '.join(temp) + ' B\n')
                        f1.write('D ' + ' '.join(temp[:-1]) + ' F\n')
                    else:
                        f1.write('D ' + ' '.join(temp) + '\n')

        if settings.FCcoords != '':
            splitFCcoords = settings.FCcoords.split(';')
            for line in splitFCcoords:
                temp = line.split(',')
                if len(temp) == 2:
                    f1.write('B ' + ' '.join(temp) + ' B\n')
                    f1.write('B ' + ' '.join(temp) + ' D\n')
                elif len(temp) == 3:
                    f1.write('A ' + ' '.join(temp) + ' B\n')
                    f1.write('A ' + ' '.join(temp) + ' D\n')
                elif len(temp) == 4:
                    f1.write('D ' + ' '.join(temp) + ' B\n')
                    f1.write('D ' + ' '.join(temp) + ' D\n')

        if (settings.Constraints != '') or (coord != '') or (settings.FCcoords != ''):
            f1.write('\n')

        if settings.SolventSurface != '':
            f1.write('Surface=' + settings.SolventSurface + ' AddSph\n')
            f1.write('\n')

        if settings.WFN:
            f1.write(Gausinp + '.wfn\n')
            f1.write('\n')

        if settings.ECP != '':
            f1.write(ComposeECP(settings, atoms))
            f1.write('\n\n')

        if CustomBasisSet == True:
            print('CUSTOM BASIS SET WILL BE USED!')
            basisfile = open(os.path.join(settings.ScriptDir, 'BasisSets', (settings.BasisSet).lower()), 'r')
            basisinp = basisfile.readlines()
            basisfile.close()
            nbasisatoms = basisinp.count('****\n')
            basis = [[] for x in range(nbasisatoms)]
            i = 0
            for line in basisinp:
                basis[i].append(line)
                if line == '****\n':
                    i += 1

            for atombasis in basis:
                if atombasis[0].split(' ')[0] in atoms:
                    f1.write(''.join(atombasis))
            f1.write('\n')
        f1.close()
    else:
        fold = Gausinp.split(str(Path('/')))[-2]
        f1 = open(Gausinp + '.in', 'w')
        f1.write('echo\n')
        f1.write('START\n')
        if molecules_dict[Gausinp.split(str(Path('/')))[-3]] > 30:
            f1.write('memory 2000 mb\n')
        elif molecules_dict[Gausinp.split(str(Path('/')))[-3]] > 22:
            f1.write('memory 1000 mb\n')
        else:
            f1.write('memory 500 mb\n')
        if charge != 0:
            f1.write(f'charge {settings.ManCharge}\n')
        f1.write('geometry noautosym noautoz\n')
        if settings.Isotopes != '':
            isomap = InterpretIsotopes(atoms, settings.Isotopes)
        else:
            isomap = ['' for a in atoms]
        natom = 0
        if not settings.ONIOM:
            for atom in conformer:
                if isomap[natom] == '':
                    f1.write(atoms[natom] + '  ' + atom[0] + '  ' + atom[1] + '  ' +
                             atom[2] + '\n')
                else:
                    f1.write(atoms[natom] + '(Iso=' + isomap[natom] + ')  ' + atom[0] + '  ' + atom[1] + '  ' +
                             atom[2] + '\n')
                natom = natom + 1
        else:
            atomsection = GetONIOMatoms(settings, atoms, conformer)
            core_allocation(atomsection)
            f1.write(atomsection)

        f1.write('end\n')
        f1.write('basis\n')
        if 'I' in atoms:
            f1.write(' * library 6-311G*\n')
        else:
            f1.write(' * library 6-31G*\n')
        f1.write('end\n')
        if settings.Functional == 'UHF':
            f1.write('scf\n')
            if fold != 'fukui':
                if settings.ManMultiplicity == 2:
                    f1.write(' doublet ; uhf\n')
            f1.write(' direct\n')
            f1.write(' maxiter 100\n')
            f1.write('end\n')
        f1.write('driver\n')
        f1.write(' xyz geo-opt\n')
        f1.write(' maxiter 1000\n')
        f1.write('end\n')
        if fold == 'reagent':
            f1.write('stepper\n')
            f1.write(' MIN\n')
            f1.write(' maxiter 1000\n')
            f1.write('end\n')
            f1.write('task scf optimize\n')
            f1.write('task scf freq')
        elif fold == 'fukui':
            f1.write('stepper\n')
            f1.write(' MIN\n')
            f1.write(' maxiter 1000\n')
            f1.write('end\n')
            f1.write('task scf optimize\n')
            f1.write('dft\n')
            f1.write(' xc b3lyp\n')
            f1.write(' fukui\n')
            f1.write(' print "Fukui information"\n')
            f1.write(' odft\n')
            f1.write('end\n')
            f1.write('task dft\n')
        else:
            if resubmission is False:
                f1.write('stepper\n')
                f1.write(' TS\n')
                f1.write(' maxiter 1000\n')
                f1.write('end\n')
                f1.write('task scf saddle\n')
                f1.write('task scf freq')
            else:
                f1.write('stepper\n')
                f1.write(' TS\n')
                f1.write(' maxiter 1000\n')
                f1.write(' convggm 1.0d-04\n')
                f1.write('end\n')
                f1.write('task scf saddle\n')
                f1.write('task scf freq')
        f1.close()


def InterpretIsotopes(atoms, isotopes):
    atomisotopes = ['' for x in atoms]
    atomsections = isotopes.split(';')
    atomsections = [asec.split(',') for asec in atomsections]
    for atomsection in atomsections:
        atomisotopes[int(atomsection[0]) - 1] = atomsection[1]

    return atomisotopes


def GetONIOMatoms(settings, atoms, conformer):
    if settings.OniomLinks != []:
        low = GetOniomLinkLow(settings, atoms, conformer)
        settings.OniomLow = low

    if (settings.OniomHigh == '') and (settings.OniomLow == ''):
        print("No ONIOM atom levels defined, quitting!")
        quit()

    if settings.OniomHigh != []:
        high = settings.OniomHigh
        low = [i + 1 for i in range(len(atoms)) if i + 1 not in high]
    else:
        low = settings.OniomLow
        high = [i + 1 for i in range(len(atoms)) if i + 1 not in low]

    linkatoms, linkstrings = GenLinks(atoms, conformer, high, low)
    atomstrings = ''
    natom = 0

    for idx, atom in enumerate(conformer):
        atomstring = atoms[natom] + '  0 ' + atom[0] + '  ' + atom[1] + '  ' + atom[2]
        if (idx + 1) in high:
            atomstring += ' H'
        else:
            atomstring += ' L'
            if (idx + 1) in linkatoms:
                linkindex = linkatoms.index(idx + 1)
                atomstring += linkstrings[linkindex]

        atomstring += '\n'

        atomstrings += atomstring
        natom = natom + 1
    return atomstrings


def GetOniomLinkLow(settings, atoms, conformer):
    obmol = BuildOBMol(atoms, conformer)

    molgraph = []
    low = []

    linkhighs = [link[0] for link in settings.OniomLinks]
    # generate the molgraph
    for atom in OBMolAtomIter(obmol):
        idx = atom.GetIdx()
        molgraph.append([])
        molgraph[idx - 1].append(idx)

        for NbrAtom in OBAtomAtomIter(atom):
            molgraph[idx - 1].append(NbrAtom.GetIdx())

    for link in settings.OniomLinks:
        # map the low level molecular structure starting at link[1]
        explored = [link[1]]
        seen = molgraph[link[1] - 1]

        for a in seen:
            if (a not in explored) and (a not in linkhighs):
                # nbs - neighbours of the atom a
                nbs = molgraph[a - 1][1:]
                seen.extend([x for x in nbs if x not in explored])
                explored.append(a)

        low.extend(explored)

    return low


def GenLinks(atoms, conformer, high, low):
    obmol = BuildOBMol(atoms, conformer)

    linkatoms = []
    linkstrings = []

    for bond in OBMolBondIter(obmol):

        beginatom = bond.GetBeginAtomIdx()
        endatom = bond.GetEndAtomIdx()

        if (beginatom in high) and (endatom in low):
            linkatoms.append(endatom)
            linkstrings.append(' H ' + str(beginatom))

        if (endatom in high) and (beginatom in low):
            linkatoms.append(beginatom)
            linkstrings.append(' H ' + str(endatom))

    return linkatoms, linkstrings


def BuildOBMol(atoms, coords):
    mol = OBMol()
    for anum, acoords in zip(atoms, coords):
        atom = OBAtom()
        atom.thisown = False
        atom.SetAtomicNum(GetAtomNum(anum))
        atom.SetVector(float(acoords[0]), float(acoords[1]), float(acoords[2]))
        mol.AddAtom(atom)

    # Restore the bonds
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()

    # mol.Kekulize()

    return mol


def core_allocation(molecules_dict, file, atoms):
    print(file)
    system_size = len(atoms)
    print(system_size)

    if file.split(str(Path('/')))[-3] in molecules_dict:
        if system_size <= 10:
            molecules_dict[file.split(str(Path('/')))[-3]] = 8
        elif system_size <= 15:
            molecules_dict[file.split(str(Path('/')))[-3]] = 12
        elif system_size <= 20:
            molecules_dict[file.split(str(Path('/')))[-3]] = 16
        elif system_size <= 25:
            molecules_dict[file.split(str(Path('/')))[-3]] = 20
        elif system_size <= 30:
            molecules_dict[file.split(str(Path('/')))[-3]] = 24
        elif system_size <= 35:
            molecules_dict[file.split(str(Path('/')))[-3]] = 28
        elif system_size >= 36:
            molecules_dict[file.split(str(Path('/')))[-3]] = 32
    return molecules_dict


def ComposeGausRoute(settings, coord, basisfile='', reagent=False):
    route = '#p '
    if not settings.ONIOM:
        route = route + settings.Functional
        if settings.ECP != '':
            route = route + '/genECP'
        elif basisfile != '':
            route = route + '/gen'
        else:
            route = route + '/' + settings.BasisSet
    else:
        route = route + 'oniom(' + settings.Functional + '/' + settings.BasisSet + ':uff)'

    route = route + ' Int=UltraFine'

    if settings.IRC:
        route = route + ' IRC=(' + settings.IRCopts + ')'
    elif settings.DFTOpt:
        if (settings.ResubsDone > 0) and (settings.KeepChk == True):
            route = route + ' Opt=(maxcycles=' + str(settings.MaxDFTOptCycles * (settings.ResubsDone + 1))
        else:
            route = route + ' Opt=(maxcycles=' + str(settings.MaxDFTOptCycles)

        if (settings.Constraints != '') or (coord != '') or (settings.FCcoords != ''):
            route = route + ',ModRedundant'
        if (settings.OptStepSize != 30):
            route = route + ',MaxStep=' + str(settings.OptStepSize)
        if (settings.KeepChk == True) and (settings.ResubsDone > 0):
            route = route + ',Restart'
        if settings.CalcFC == True:
            route = route + ',CalcFC'
        if settings.TS == True:
            if reagent == False:
                route = route + ',ts'
            else:
                route = route
        if settings.NoEigen == True:
            route = route + ',noeigentest'
        if settings.RFO == True:
            route = route + ',RFO'
        route = route + ')'

    if settings.SCF != '':
        route = route + ' SCF=' + settings.SCF

    if settings.GuessChk:
        route = route + ' Guess=Read'

    if settings.WFN:
        route = route + ' Output=wfn'

    if settings.CheckStability:
        route = route + ' Stable=Opt'

    if settings.EmpDisp:
        route = route + ' EmpiricalDispersion=GD3'

    if settings.Solvent != '':
        route = route + ' scrf=(solvent=' + settings.Solvent
        if settings.SCRFopts != '':
            route = route + ',' + settings.SCRFopts
        if settings.SolventSurface != '':
            route = route + ',' + 'read'
        route = route + ')'
    if settings.Freq == True:
        route = route + ' Freq'

    route = route + ' MaxDisk=' + str(settings.MaxDisk) + 'GB'
    route = route + ' scf(xqc) Freq'
    route = route + '\n'

    return route


def ComposeECP(settings, atoms):
    lightels = []
    heavyels = []

    ECPlines = []

    for a in atoms:
        if (GetAtomNum(a) <= 18) and (a not in lightels):
            lightels.append(a)
        if (GetAtomNum(a) > 18) and (a not in heavyels):
            heavyels.append(a)

    ECPlines.append(' '.join(lightels) + ' 0')
    ECPlines.append(settings.BasisSet)
    ECPlines.append('****')
    ECPlines.append(' '.join(heavyels) + ' 0')
    ECPlines.append(settings.ECP)
    ECPlines.append('****')
    ECPlines.append('')
    ECPlines.append(' '.join(heavyels) + ' 0')
    ECPlines.append(settings.ECP)

    return '\n'.join(ECPlines)


def GetFiles2Run(GinpFiles, settings):
    # GinpFiles contain the names of all relevant input files

    Files2Run = []
    # for every input file check that there is a completed output file,
    # delete the incomplete outputs and add the inputs to be done to Files2Run
    for filename in GinpFiles:
        print(filename)
        if not os.path.exists(filename[:-3] + 'out'):
            Files2Run.append(filename)
        else:
            if IsGausCompleted(filename[:-3] + 'out'):
                # print filename[:-3]+'out already exists'
                continue
            else:
                os.remove(filename[:-3] + 'out')
                Files2Run.append(filename)

    return Files2Run


def IsGausCompleted(f):
    Gfile = open(f, 'r')
    outp = Gfile.readlines()
    Gfile.close()
    if len(outp) < 10:
        return False
    if ("Normal termination" in outp[-1]) or ('termination' in outp[-3] and 'l9999.exe' in outp[-3]):
        return True
    else:
        return False


def RunCalcs(GausJobs, settings):
    print('\nRunning Gaussian locally...')

    NCompleted = 0
    GausPrefix = "g16 < "

    for f in GausJobs:
        time.sleep(3)
        print(GausPrefix + f + ' > ' + f[:-3] + 'out')
        try:
            outp = subprocess.check_output(GausPrefix + f + ' > ' + f[:-3] + 'out', shell=True)
        except Exception as e:
            outp = str(e.output)
        NCompleted += 1
        print("Gaussian job " + str(NCompleted) + " of " + str(len(GausJobs)) + \
              " completed.")


def CheckOutputFiles(GausJobs, functional):
    for f in GausJobs:
        GOutpFile = f[:-3] + '.out'
        gausfile = open(GOutpFile, 'r')

        GOutp = gausfile.readlines()
        gausfile.close()

        # Gather data from the output file
        Energies = []
        ConvRep = []
        OutOfSteps = ''
        StatPFound = ''
        StepN = ''
        Term = ''
        for i, line in enumerate(GOutp):
            if ('E(R' + functional.upper() + ')' in line) or \
                    ('E(U' + functional.upper() + ')' in line):
                temp = [_f for _f in line.split(" ") if _f]
                Energies.append(' '.join(temp[2:5]))
            if 'Item' in line:
                ConvRep.extend(GOutp[i:i + 6])
            if 'Number of steps exceeded' in line:
                OutOfSteps = line
            if 'Stationary' in line:
                StatPFound = line
            if 'Step number' in line:
                StepN = line
            if 'termination' in line:
                Term = line

        # Format and print the data
        title = "Report for " + GOutpFile
        print('\n' + '-' * len(title) + '\n' + title + '\n' + '-' * len(title) + '\n')
        if StepN != '':
            print(StepN)
        print(''.join(ConvRep[-6:]))
        print("Last 3 energies: \n" + '\n'.join(Energies[-3:]))
        # print "Last 3 energies: \n" + ''.join(Energies)
        if OutOfSteps != '':
            print(OutOfSteps)
        if StatPFound != '':
            print(StatPFound)
        if Term != '':
            print(Term)


def CheckConvergence(GausInpFiles):
    # FilesRun - list of files of the form input.com
    # we are looking for the corresponding optimization output files
    # in the form inputtemp.out
    GoutpFiles = []
    for filename in GausInpFiles:
        GoutpFiles.append(filename[:-3] + '.out')

    Nunconverged = 0
    unconverged = []
    for outfile in GoutpFiles:
        f = open(outfile, 'r')
        ginp = '\n'.join(f.readlines())
        f.close()
        if not 'Stationary point found' in ginp:
            print('STATIONARY POINT NOT FOUND')
            Nunconverged += 1
            unconverged.append(outfile)
    return unconverged


def ReadSDFGeometry(SDfile):
    from openbabel import OBMol, OBConversion, OBMolAtomIter

    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()

    print("Reading " + SDfile)

    obconversion.ReadFile(obmol, SDfile)

    obmol.ConnectTheDots()

    atoms = []
    coords = []
    charge = obmol.GetTotalCharge()

    for atom in OBMolAtomIter(obmol):
        x_str = format(atom.GetX(), '.6f')
        y_str = format(atom.GetY(), '.6f')
        z_str = format(atom.GetZ(), '.6f')
        coords.append([x_str, y_str, z_str])
        atoms.append(GetAtomSymbol(atom.GetAtomicNum()))
    return atoms, coords, charge


def ReadGeometry(GOutpFile, settings):
    gausfile = open(GOutpFile, 'r')
    print("Reading " + GOutpFile)
    GOutp = gausfile.readlines()

    index = 0
    atoms = []
    coords = []
    gindexes = []
    chindex = None

    # Find the geometry section and charge section
    for index in range(len(GOutp)):
        if 'Charge =' in GOutp[index]:
            chindex = index
        if ('Input orientation:' in GOutp[index]) or ("Standard orientation:" in GOutp[index]):
            gindexes.append(index + 5)

    # Read geometries
    print('Reading geometry ' + str(settings.GOutpGN) + ' for file ' + GOutpFile)
    for line in GOutp[gindexes[settings.GOutpGN]:]:
        if '--------------' in line:
            break
        else:
            data = [_f for _f in line[:-1].split(' ') if _f]
            atoms.append(GetAtomSymbol(int(data[1])))
            coords.append(data[3:])

    if chindex != None:
        line = GOutp[chindex].split('Charge = ')
        line = line[1].split(' Multiplicity = ')
        charge = int(line[0])
    else:
        charge = -1000

    gausfile.close()

    return atoms, coords, charge


def CopyGoutpFiles(Unconverged, settings):
    cwd = os.getcwd()
    outfiles = []

    for f in Unconverged:
        if settings.confs != []:
            number = f.split(settings.suffix)[1]
            number = number.split('.')[0]
            number = number.lstrip('0')
        else:
            number = ''
        if settings.ResubsDone > 0:
            base = f.split(settings.suffix)[0][:-1]
            base = base + number
        else:
            base = f.split(settings.suffix)[0]
            base = base + number
        newfilename = base + chr(98 + settings.ResubsDone) + '.out'
        print("Gaussian.CopyGoutpFiles: newfilename = " + newfilename)
        shutil.copyfile(cwd + '/' + f, cwd + '/' + newfilename)
        outfiles.append(newfilename)
        print("Gaussian.CopyGoutpFiles: outfiles = " + str(outfiles))
        for i in range(len(settings.ChkFiles)):
            if settings.ChkFiles[i][0] == f[:-4]:
                settings.ChkFiles[i][0] = newfilename[:-4] + settings.suffix + '001'

    return outfiles, settings.ResubsDone + 1, settings.ChkFiles


PTable = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', \
          'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', \
          'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
          'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', \
          'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', \
          'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
          'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']


def GetAtomSymbol(AtomNum):
    if AtomNum > 0 and AtomNum < len(PTable):
        return PTable[AtomNum - 1]
    else:
        print("No such element with atomic number " + str(AtomNum))
        return 0


def GetAtomNum(AtomSymbol):
    if AtomSymbol in PTable:
        return PTable.index(AtomSymbol) + 1
    else:
        print("No such element with symbol " + str(AtomSymbol))
        return 0


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
