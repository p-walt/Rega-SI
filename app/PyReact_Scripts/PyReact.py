#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyReact integrated workflow for the running DFT calculations
v0.5

Please cite as:
PyReact v0.5, K. Ermanis, Cambridge, 2021

Copyright (c) 2017-2018 Kristaps Ermanis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Created on Wed Oct 09 15:26:32 2017

@author: K. Ermanis

The main file, that should be called to start the PyReact workflow.
Interprets the arguments and takes care of the general workflow logic.
"""

try:
    import Gaussian
except ModuleNotFoundError:
    import PyReact_Scripts.Gaussian as Gaussian

import sys
import os
import argparse
import getpass
import datetime
import importlib

DFTpackages = [['g', 'd', 'a'],['Gaussian', 'GaussianDarwin', 'GaussianAugusta']]

#Assigning the config default values
class Settings:
    DFT = 'z'
    G16 = True
    Title = 'Gausmol'
    DarwinNodeSize = 32
    Solvent = ''
    SCRFopts = ''
    SolventSurface=''
    DFTOpt = True
    MaxDFTOptCycles = 50
    OptStepSize = 30
    TimeLimit = 24
    Resub = 0
    ResubsDone = 0
    SubmitOnly = False
    queue = ''
    project = 'GOODMAN-SL3-CPU'
    OBPath = '/home/ke291/Tools/openbabel-install/lib/python2.7/site-packages/'
    SCHRODINGER = '/usr/local/shared/schrodinger/current'
    folders = [] #Keep track of job folders, useful for reusing chk in auto-resubmissions
    ChkFiles = [] #Mapping between input files and chk files [inp, chk], useful for reusing chk in auto-resubmissions
    nProc = 1
    ScriptDir = ''
    user = 'ke291'
    t2scrf = '/rds/rds-t2-cs098/'
    t2project = 'T2-CS098-CPU'
    suffix = 'outp'
    GeomSource = ''
    MaxConcurrentJobs = 75
    MaxConcurrentJobsDarwin = 320
    GenOnly = False
    GenOnlyQ = False
    ManCharge = -1000
    ManMultiplicity = -1000
    Charge = -1000
    BasisSet = "6-31g(d,p)"
    ECP = ''
    Functional = "b3lyp"
    SCF = ''
    GuessChk = False
    CheckStability = False
    EmpDisp = ''
    GOutpInputs = False
    GOutpGN = -1
    SDFInputs = False
    GausInpFiles = []
    confs = []
    Constraints = ''
    Isotopes = ''
    Scan = ''
    IRC = False
    IRCopts = 'calcfc,maxpoints=10,recalc=12,maxcycle=40,cartesian'
    ONIOM = False
    OniomHigh = []
    OniomLow = []
    OniomLinks = []
    Restart=False
    Freq=False
    KeepChk = False
    MaxDisk = 100   # Max scratch disk available per job in gb
    WFN = False
    CalcFC = False
    FCcoords = []
    TS = False
    NoEigen = False
    RFO = False
    Rega = False
    inp_file = ''

settings = Settings()


def PyReact(filenames, settings):
    print('FILENAMES')
    inpfiles = filenames
    if settings.Rega != '':
        print(inpfiles)
        wholemolecule = True
    else:

        wholemolecule = False

    #Allowed inputs - .sdf files for calcs from scratch,
    # and Gaussian *.out files to use Gaussian geometries for other calcs
    for f in inpfiles:
        if (f[-4:] != '.sdf') and (f[-4:] != '.out') and (f[-4:] != '.txt'):
            print(f + ' is an invalid file type. Exiting...')
            quit()
        elif f[-4:] == '.txt':
            with open(f, 'r') as file:
                files_to_run = file.readlines()
                for line in files_to_run:
                    li = line.strip('\n')
                    inpfiles.append(li)
            inpfiles.pop(0)
            break
    print(f'INPFILES ARE {inpfiles}')
    # Import the appropriate DFT software interface code file -
    # one of Gaussian.py, GaussianDarwin.py, GaussianAugusta.py etc
    DFT = ImportDFT(settings.DFT)

    #Generate DFT input files
    print('\nRunning DFT setup...')
    settings.GausInpFiles, molecules_dict = DFT.SetupGaussian(inpfiles, settings, wholemolecule)
    QRun = False

    Files2Run = DFT.GetFiles2Run(settings.GausInpFiles, settings)
    print(f'FILES TO RUN:{Files2Run}')
    if len(Files2Run) == 0:
        QRun = True

    if settings.GenOnlyQ:
        print("Gaussian input files generated, quitting as instructed.")
        quit()

    if QRun:
        print('DFT has already been run for these inputs. Skipping...')
    else:
        DFT.RunCalcs(Files2Run, settings, molecules_dict)

    #Report on the outcome of the calculations
    DFT.CheckOutputFiles(settings.GausInpFiles, settings.Functional)

# Selects which DFT package to import, returns imported module
def ImportDFT(dft):
    if dft in DFTpackages[0]:
        DFTindex = DFTpackages[0].index(dft)
        DFT = importlib.import_module(DFTpackages[1][DFTindex])
    else:
        print("Invalid DFT package selected")
        quit()

    return DFT


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


if __name__ == '__main__':

    print("==========================")
    print("PyReact script\nv0.5")
    print("\nCopyright (c) 2017-2021 Kristaps Ermanis")
    print("Distributed under GNU GPLv3 license")
    print("Please cite as:")
    print("PyReact v0.5, K. Ermanis, Nottingham, 2021")
    print("==========================\n\n")

    parser = argparse.ArgumentParser(description='PyReact script to setup\
    and run Gaussian locally and on clusters')
    parser.add_argument('-d', '--dft', help="Select DFT program default is g (Gaussian locally)",
        choices=DFTpackages[0], default='g')
    parser.add_argument('StructureFiles', nargs='+', default=[], help=
    "One or more SDF file for the structures to be verified by DP4. At least one\
    is required, if automatic diastereomer and tautomer generation is used.\
    One for each candidate structure, if automatic generation is not used")
    parser.add_argument("-s", "--solvent", help="Specify solvent to use\
    for dft calculations")
    parser.add_argument('--SCRFopts', help="Supply extra options for implicit solvent models",
                        default=settings.SCRFopts)
    parser.add_argument("-q", "--queue", help="Specify queue for job submission\
    on ziggy", default=settings.queue)
    parser.add_argument("--project", help="Specify project for job submission\
    on darwin", default=settings.project)
    parser.add_argument("--suffix", help="Specify suffix added to gaussian input/output files",
                        default=settings.suffix)
    parser.add_argument("--TimeLimit", help="Specify job time limit for jobs\
    on ziggy or darwin", type=int)
    parser.add_argument("--nProc", help="Specify number of processor cores\
    to use for Gaussian calculations", type=int, default=1)
    parser.add_argument("--batch", help="Specify max number of jobs per batch",
    type=int, default=settings.MaxConcurrentJobs)
    parser.add_argument("-c", "--confs", help="Specify the conformers \
    from the conformation search to carry through to DFT, seperate numbers with commas,\
    different structure lists with semicolons. '1,2;1,2' will do DFT calculation on 2 \
    lowest energy conformations for 2 structures", default='')
    parser.add_argument("--constr", help="Specify distance/angle/dihedral constraints in g09 format, \
    with multiple inputs only make sense if the atom numbering is the same", default='')
    parser.add_argument("--ConstrF", help="Specify file containing constraints in g09 format, \
    with multiple inputs only make sense if the atom numbering is the same", default='')
    parser.add_argument("--scan", help="Specify 1 distance/angle/dihedral and a series of desired values, \
    example: --scan '1,2;2.0,2.2,2.4'", default='')
    parser.add_argument("--isotopes", help="Specify non-standard isotopes for atoms in the format " +
                                           "'AtomN1,AtomM1;AtomN2,AtomM2'", default='')
    parser.add_argument("--freq", help="Do DFT frequency calculation", action="store_true")
    parser.add_argument("--Restart", help="Restart previously failed single point job based on job log infor",
                        action="store_true")
    parser.add_argument("--SubmitOnly", help="Exit immediately after submitting a job to cluster " +
                        "instead of monitoring and downloading after finishing. Has no effect on local jobs.",
                        action="store_true")
    parser.add_argument("--KeepChk", help="Don't delete and copy back the .chk files", action="store_true")
    parser.add_argument("--WFN", help="Write WFN file", action="store_true")
    parser.add_argument("--GuessChk", help="Get SCF guess from previous chk file", action="store_true")
    parser.add_argument("--CheckStability", help="Check/reoptimize wavefunction", action="store_true")
    parser.add_argument("--SolventSurface", help="Define solvent surface option", default=settings.SolventSurface)
    parser.add_argument("--FC", help="Calculate force constants before optimization", action="store_true")
    parser.add_argument("--FCcoords", help="Specify coordinates to calculate numerical force constansts for," + \
        "example: --FCcoords '34,42;32,34'", default='')
    parser.add_argument("--IRC", help="Calculate IRC", action="store_true")
    parser.add_argument("--IRCopts", help="Specify additional options for IRC calculation", default=settings.IRCopts)
    parser.add_argument("--TS", help="Do transition state search", action="store_true")
    parser.add_argument("--NoEigen", help="Disable the check for the correct number of eigenvalues", action="store_true")
    parser.add_argument("--RFO", help="Use RFO geometry optimisation algorithm instead of GEDIIS", action="store_true")
    parser.add_argument("-l", "--ConfLimit", help="Specify maximum number of \
    conformers per structure. If above this, adaptive RMSD pruning will be \
    performed", type=int, default=100)
    parser.add_argument("-g", "--GenOnly", help="Only generate diastereomers\
    and tinker input files, but don't run any calculations", action="store_true")
    parser.add_argument("--GenOnlyQ", help="Only generate the DFT input files, then quit", action="store_true")
    parser.add_argument('--SinglePoint', help="Optimize geometries at DFT\
    level", action="store_true")
    parser.add_argument('--G16', help="Use Gaussian 16", action="store_true")
    parser.add_argument("--OptCycles", help="Specify max number of DFT geometry\
    optimization cycles", type=int, default=settings.MaxDFTOptCycles)
    parser.add_argument("--OptStep", help="Specify the max step size\
    Gaussian should take in optimization, default is 30", type=int, default=settings.OptStepSize)
    parser.add_argument("--GOutpGN", help="Specify which geometry should be used from the\
    Gaussian output file, default is -1 (the last one)", type=int, default=settings.GOutpGN)
    parser.add_argument("--MaxDisk", help="Specify the max available scratch disk space per job in gb",
                        type=int, default=settings.MaxDisk)
    parser.add_argument("--Resub", help="Specify the max number of times\
    the optimizations can be resubmitted, if unconverged", type=int, default=settings.Resub)
    parser.add_argument("--EmpDisp", help="Specify the use of empirical dispersion correction", action="store_true")
    parser.add_argument('--SCF', help="Supply directives for SCF procedure (in case of SCF convergence problems)",
                        default=settings.SCF)
    parser.add_argument('-B', '--BasisSet', help="Selects the basis set for\
    DFT calculations", default='6-31g(d,p)')
    parser.add_argument('-F', '--Functional', help="Selects the functional for\
    DFT calculations", default='b3lyp')
    parser.add_argument('--ECP', help="Selects the ECP basis set for heavy elements", default=settings.ECP)
    parser.add_argument("--ONIOMlow",
                        help="Comma-separated list of the atoms to put in the UFF layer in ONIOM calulation",
                        default='')
    parser.add_argument("--ONIOMhigh",
                        help="Comma-separated list of the atoms to put in the DFT layer in ONIOM calulation",
                        default='')
    parser.add_argument("--ONIOMlinks",
                        help="List of bonds linking high and low level regions in ONIOM calculation, format: '12,13;34,35'",
                        default='')
    parser.add_argument("--ONIOMfile",
                        help="file containing comma-separated lists for atoms to put in high and low levels" + \
                        " of ONIOM calculation", default='')
    parser.add_argument("--charge", help="Manually specify overall charge of the" +\
    " molecule for DFT calcs", type=int, default=settings.ManCharge)
    parser.add_argument("--mult", help="Manually specify overall multiplicity of the" + \
                                         " molecule for DFT calcs", type=int, default=settings.ManMultiplicity)
    parser.add_argument("--Mol", help="the name of the molecule being calculated", default='')
    parser.add_argument("--Rega", help="flag to start calculation as part of the Rega software package", action='store_true')
    parser.add_argument("--nwchem", help="specifies that the NWChem software package should be used for calculations "
                                         "rather than the default Gaussian", action="store_true")
    parser.add_argument("--array", help="flag to submit a slurm array rather than individual slurm scripts for each "
                                        "file", action='store_true')
    args = parser.parse_args()

    settings.suffix = args.suffix
    settings.compound = args.Mol
    settings.Rega = args.Rega
    settings.array = args.array
    settings.Title = settings.compound
    settings.DFT = args.dft
    settings.queue = args.queue
    settings.project = args.project
    settings.ScriptDir = getScriptPath()
    settings.BasisSet = args.BasisSet
    settings.ECP = args.ECP
    settings.Functional = args.Functional
    settings.SCF = args.SCF
    settings.SCRFopts = args.SCRFopts
    settings.SolventSurface = args.SolventSurface
    settings.IRCopts = args.IRCopts
    settings.nProc = args.nProc
    settings.MaxConcurrentJobs = args.batch
    settings.MaxDFTOptCycles = args.OptCycles
    settings.OptStepSize = args.OptStep
    settings.GOutpGN = args.GOutpGN
    settings.Resub = args.Resub
    settings.ManCharge = args.charge
    settings.ManMultiplicity = args.mult
    settings.MaxDisk = args.MaxDisk

    settings.Constraints = args.constr
    settings.Scan = args.scan
    settings.FCcoords = args.FCcoords

    settings.Isotopes = args.isotopes

    if args.ONIOMhigh != '':
        settings.ONIOM = True
        settings.OniomHigh = [int(x) for x in args.ONIOMhigh.split(',')]

    if args.ONIOMlow != '':
        settings.ONIOM = True
        settings.OniomLow = [int(x) for x in args.ONIOMlow.split(',')]

    if args.ONIOMlinks != '':
        settings.ONIOM = True
        temp = args.ONIOMlinks.split(';')
        for link in temp:
            settings.OniomLinks.append([int(x) for x in link.split(',')])

    if args.ONIOMfile != '':
        settings.ONIOM = True
        OniomFile = open(args.ONIOMfile, 'r')
        OniomData = OniomFile.readlines()
        OniomFile.close()
        if len(OniomData) < 2:
            print("Invalid ONIOM file, quitting...")
            quit()
        if len(OniomData[0][:-1]) != 0:
            settings.OniomHigh = [int(x) for x in OniomData[0][:-1].split(',')]

        if len(OniomData[1][:-1]) != 0:
            settings.OniomLow = [int(x) for x in OniomData[1][:-1].split(',')]

        if len(OniomData) > 2:
            temp = OniomData[2][:-1].split(';')
            for link in temp:
                settings.OniomLinks.append([int(x) for x in link.split(',')])

    if settings.DFT == 'd' and not args.TimeLimit:
        print("For calculations on Darwin explicit time limit in hours " + \
            "must be specified, exiting...")
        quit()
    if args.TimeLimit:
        settings.TimeLimit = args.TimeLimit
    if args.NoEigen:
        settings.NoEigen = True
    if args.TS:
        settings.TS = True
    if args.FC:
        settings.CalcFC = True
    if args.IRC:
        settings.IRC = True
    if args.RFO:
        settings.RFO = True
    if args.KeepChk:
        settings.KeepChk = True
    if args.WFN:
        settings.WFN = True
    if args.GuessChk:
        settings.GuessChk = True
    if args.CheckStability:
        settings.CheckStability = True
    if args.freq:
        settings.Freq = True
    if args.SubmitOnly:
        settings.SubmitOnly = True
    if args.Restart:
        settings.Restart = True
    if args.EmpDisp:
        settings.EmpDisp = True
    if args.SinglePoint:
        settings.DFTOpt = False
    if args.G16:
        settings.G16 = True
    if args.GenOnly:
        settings.GenOnly = True
    if args.GenOnlyQ:
        settings.GenOnlyQ = True
    if args.solvent:
        settings.Solvent = args.solvent
    if args.nwchem:
        settings.nwchem = True
    else:
        settings.nwchem = False

    now = datetime.datetime.now()
    settings.now = now.strftime('%d%b%H%M')

    settings.user = getpass.getuser()

    inpfiles = args.StructureFiles
    jobs_txt = inpfiles.copy()
    settings.inp_file = jobs_txt[0]
    print('INITIAL')
    print(inpfiles)
    print(f'SETTINGS.INP_FILE {settings.inp_file}')
    with open('cmd.log', 'a') as f:
        f.write(' '.join(sys.argv) + '\n')
    print(settings.Rega)
    PyReact(inpfiles, settings)

    unconverged = Gaussian.CheckConvergence(settings.GausInpFiles)  #Check convergence

    print(str(len(unconverged)) + " files have not converged.")

    while ((settings.Resub - settings.ResubsDone) > 0) and (len(unconverged) > 0):

        inpfiles, settings.ResubsDone, settings.ChkFiles = Gaussian.CopyGoutpFiles(unconverged, settings)
        print("PyReact: inpfiles = " + str(inpfiles))
        settings.confs = []
        settings.GOutpGN = -1
        now = datetime.datetime.now()
        settings.now = now.strftime('%d%b%H%M')
        PyReact(inpfiles, settings)
        unconverged = Gaussian.CheckConvergence(settings.GausInpFiles)
        print(str(len(unconverged)) + " files have not converged.")

    if (settings.Resub > 0) and (len(unconverged) > 0):
        print("Max number of resubmissions reached.")
