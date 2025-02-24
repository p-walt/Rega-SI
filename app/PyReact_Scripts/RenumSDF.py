#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyReact integrated workflow for the running of MM, DFT calculations
v0.1

Copyright (c) 2020 Kristaps Ermanis

Script that takes in a atom renumbering map and an sdf file(s),
outputs a renumbered copy of the sdf file(s)
"""



import argparse
import glob
import shutil
import os

from openbabel import *
from openbabel import OBMol, OBAtom, OBMolAtomIter, OBAtomAtomIter, OBAtomBondIter, OBMolBondIter, OBConversion

def PrintObMolAtoms(obmol):
    atoms = []
    for atom in OBMolAtomIter(obmol):
        atoms.append(GetAtomSymbol(atom.GetAtomicNum()) + str(atom.GetIdx()))
    print(','.join(atoms))


def RenumSDF(f, renum_map):

    #Read molecule from file
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    oldmol = OBMol()
    obconversion.ReadFile(oldmol, f)

    natoms = oldmol.NumAtoms()

    if len(renum_map) < natoms:
        print('Warning: Renumbering map has less atoms than the molecule in ' + f)
        print('The atoms not in the map will be added at the end.')

    temp = []
    anums = []

    added = []
    newmol = OBMol()
    for anum in renum_map:
        if anum > natoms:
            print('Error: ' + f + ' does not have atom with index ' + str(anum))
            quit()
        newmol.AddAtom(oldmol.GetAtom(anum))
        added.append(anum)

    if oldmol.NumAtoms() > newmol.NumAtoms():
        for obatom in OBMolAtomIter(oldmol):
            if obatom.GetIdx() not in added:
                newmol.AddAtom(obatom)

    #Restore the bonds
    newmol.ConnectTheDots()
    newmol.PerceiveBondOrders()

    print('\nOLD NUMBERING FOR ' + f + ':')
    PrintObMolAtoms(oldmol)
    print(oldmol.NumAtoms())
    print('\nNEW NUMBERING FOR ' + f + ':')
    PrintObMolAtoms(newmol)
    print(newmol.NumAtoms())

    #Write renumbered molecule to file
    obconversion.SetOutFormat("sdf")
    obconversion.WriteFile(newmol, f[:-4] + 'r.sdf')


PTable = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', \
          'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', \
          'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
          'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', \
          'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', \
          'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
          'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']


def GetAtomSymbol(AtomNum):

    if AtomNum > 0 and AtomNum < len(PTable):
        return PTable[AtomNum-1]
    else:
        print("No such element with atomic number " + str(AtomNum))
        return 0


def GetAtomNum(AtomSymbol):

    if AtomSymbol in PTable:
        return PTable.index(AtomSymbol)+1
    else:
        print("No such element with symbol " + str(AtomSymbol))
        return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for renaming files')
    parser.add_argument('RenumMap', help="Renumbering map in the form '3,2,1' " +
                                         "meaning that old atom 3 will be the new atom 1 etc.")

    parser.add_argument('StructFiles', nargs='+', default=[], help=
    "One or more SDF file for the structures to be renumbered. The atom numbering should be the same")

    args = parser.parse_args()
    print(args.RenumMap)
    print(args.StructFiles)

    for f in args.StructFiles:
        if f[-4:] != '.sdf':
            print('Error: ' + f + ' has invalid type')
            quit()
        if not os.path.exists(f):
            print('Error: ' + f + ' does not exist')
            quit()

    RenumMap = [_f for _f in args.RenumMap.split(',') if _f]
    RenumMap = [int(x) for x in RenumMap]
    for f in args.StructFiles:
        RenumSDF(f, RenumMap)
