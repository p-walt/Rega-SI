#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyReact integrated workflow for the running of MM, DFT calculations
v0.1

Copyright (c) 2017 Kristaps Ermanis

Script that reads MacroModel conf search files and groups the
conformations by the chiral configuration of any 4 given atoms
"""

#Assigning the config default values

# A mini settings class, only what's needed for ReadMacromodel to work
class Settings:
    Rot5Cycle = False
    SCHRODINGER = '/usr/local/shared/schrodinger/current'
    MaxCutoffEnergy = 1000.0        #We set it very high so that all conformations are read
    confs = ''

settings = Settings()

import MacroModel
import argparse
import os
from math import sqrt, acos

def main(filename, ChirAtoms, AtomNums, CoordVal):

    atoms, conformers, charge = MacroModel.ReadMacromodel(filename, settings)

    print("Molecule has " + str(len(atoms)) + " atoms")
    print(str(len(conformers)) + " conformers read")

    for an in ChirAtoms:
        print(atoms[an-1])

    Rconfs = []
    Sconfs = []
    NDconfs = []

    tempAtomNums = [x-1 for x in AtomNums]

    for i, conf in enumerate(conformers):

        res = AssignConfig(conf, ChirAtoms)
        if res > 0 and CoordBelowVal(conf, tempAtomNums, CoordVal):
            Rconfs.append(i+1)
        elif res < 0 and CoordBelowVal(conf, tempAtomNums, CoordVal):
            Sconfs.append(i+1)
        elif CoordBelowVal(conf, tempAtomNums, CoordVal):
            NDconfs.append(i+1)

    print("Rconfs: " + str(Rconfs) + "\n")
    print("Sconfs: " + str(Sconfs) + "\n")
    print("NDconfs: " + str(NDconfs) + "\n")

    print("There are " + str(len(Rconfs)) + " R conformations")
    print("There are " + str(len(Sconfs)) + " S conformations")
    print("There are " + str(len(NDconfs)) + " conformations with undetermined configuration")


def AssignConfig(conf, ChirAtoms):

    AtomCoords = []
    for an in ChirAtoms:
        temp = conf[an-1]
        temp = [float(x) for x in temp[1:]]
        AtomCoords.append(temp)
    A, B, C, D = AtomCoords
    Cent = [(A[0] + B[0] + C[0])/3, (A[1] + B[1] + C[1]) / 3, (A[2] + B[2] + C[2]) / 3]

    RefV = P2Vec(D, Cent)
    V1 = P2Vec(Cent, A)
    V2 = P2Vec(Cent, B)
    V3 = P2Vec(Cent, C)

    cp1 = crossproduct(V1, V2)
    cp2 = crossproduct(V2, V3)

    if (dotproduct(cp1, RefV) > 0) and (dotproduct(cp2, RefV)>0):
        return 1
    elif (dotproduct(cp1, RefV) < 0) and (dotproduct(cp2, RefV) < 0):
        return -1
    else:
        return 0


def CoordBelowVal(conf, coords, val):

    if coords != []:
        if len(coords) == 2:
            res = CalcDist(conf, coords)
        elif len(coords) == 3:
            res = CalcAngle(conf, coords)
        else:
            res = CalcDihed(conf, coords)

        if res < val:
            return True
        else:
            return False
    else:
        return True


def CalcDist(conf, atoms):

    geom = []
    for a in conf:
        geom.append([float(x) for x in a[1:]])

    a1 = geom[atoms[0]]
    a2 = geom[atoms[1]]

    vec = P2Vec(a1, a2)

    return length(vec)


def CalcAngle(geom, atoms):

    a1 = geom[atoms[0]]
    a2 = geom[atoms[1]]
    a3 = geom[atoms[2]]

    vec1 = P2Vec(a2, a1)
    vec2 = P2Vec(a2, a3)

    return 180*angle(vec1, vec2)/pi


def CalcDihed(geom, atoms):
    print("Not implemented yet")
    quit()
    return 0


def P2Vec(p1, p2):
    return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]

def crossproduct(v1, v2):
    product = [0, 0, 0]
    product[0] = v1[1]*v2[2]-v1[2]*v2[1]
    product[1] = v1[2]*v2[0]-v1[0]*v2[2]
    product[2] = v1[0]*v2[1]-v1[1]*v2[0]
    return product


def dotproduct(v1, v2):
        return sum((a*b) for a, b in zip(v1, v2))


def length(v):
    return sqrt(dotproduct(v, v))


def angle(v1, v2):
    return acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


#Given 3 atoms, finds a plane defined by a normal vector and d
def FindPlane(atom1, atom2, atom3):

    vector1 = [atom2.x() - atom1.x(), atom2.y() - atom1.y(),
               atom2.z() - atom1.z()]
    vector2 = [atom3.x() - atom1.x(), atom3.y() - atom1.y(),
               atom3.z() - atom1.z()]
    cross_product = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
                     -1 * vector1[0] * vector2[2] + vector1[2] * vector2[0],
                     vector1[0] * vector2[1] - vector1[1] * vector2[0]]

    d = cross_product[0] * atom1.x() - cross_product[1] * atom1.y() + \
        cross_product[2] * atom1.z()

    return cross_product, d


#Calculates distance from an atom to a plane
def PointPlaneDist(norm, d, atom):

    point = []

    point.append(atom.x())
    point.append(atom.y())
    point.append(atom.z())

    a = norm[0]*point[0] + norm[1]*point[1] + norm[2]*point[2] + d
    b = sqrt(norm[0]**2 + norm[1]**2 + norm[2]**2)

    return a/b


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Enantiomeric conformer grouping script')
    parser.add_argument('ConfSearchFile', help="Base name for Macromodel output files")
    parser.add_argument('AtomsList', help="Comma seperated list of 4 atom numbers, that should" +\
                        " be used for enantiomeric grouping. Numbering starts with 1, atom" + \
                        " priorities assumed to be the same as the order of the atoms in the list")
    parser.add_argument("--MaxC", help="Comma seperated list of 2 to 4 atom numbers, that define" +\
                        " the coordinate of interest, followed by a value. All conformers with the" +\
                        "coordinate value above this will be ignored", default="")
    args = parser.parse_args()
    print(args.ConfSearchFile)
    print(args.AtomsList)
    print(args.MaxC)

    if args.MaxC != "":
        MaxCList = args.MaxC.split(',')
        AtomNums = []
        for x in MaxCList[:-1]:
            try:
                AtomNums.append(int(x))
            except ValueError:
                print("Atom numbers must be integers.")
                quit()
        try:
            CoordVal = float(MaxCList[-1])
        except ValueError:
            print("Coordinate value must be float.")
            quit()

        if (len(AtomNums) > 4) or (len(AtomNums) < 2):
            print("Invalid number of atoms, must be between 2 and 4")
            quit()
    else:
        AtomNums = []
        CoordVal = -100000

    SchrodEnv = os.getenv('SCHRODINGER')
    if SchrodEnv != None:
        settings.SCHRODINGER = SchrodEnv

    ChirAtoms = [int(x) for x in args.AtomsList.split(',')]

    main(args.ConfSearchFile, ChirAtoms, AtomNums, CoordVal)
