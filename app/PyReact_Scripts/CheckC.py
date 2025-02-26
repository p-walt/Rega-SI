#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyReact integrated workflow for the running of MM, DFT calculations
v0.1

Copyright (c) 2017 Kristaps Ermanis

Script that reads Gaussian output files and outputs the value
of a specified internal coordinate for all geometries in the file
"""

import argparse
import os
from math import sqrt, acos, pi

def main(filename, coords):

    atoms, geometries = ReadGeometries(filename)

    print("Molecule has " + str(len(atoms)) + " atoms")
    print(str(len(geometries)) + " geometries read")

    results = []

    for g in geometries:
        if len(coords) == 2:
            results.append(CalcDist(g, [x-1 for x in coords]))
        elif len(coords) == 3:
            results.append(CalcAngle(g, [x-1 for x in coords]))
        else:
            results.append(CalcDihed(g, [x-1 for x in coords]))

    CoordStr = ','.join([atoms[x-1] + str(x) for x in coords])
    if len(coords) == 2:
        print("The distances between atoms " + CoordStr + ":")
        print(', '.join([format(x,"5.3f") for x in results]))
    elif len(coords) == 3:
        print("The angles between atoms " + CoordStr + ":")
        print(', '.join([format(x,"4.1f") for x in results]))
    else:
        print("The dihedrals between atoms " + CoordStr + ":")
        print(', '.join([format(x,"4.1f") for x in results]))


def ReadGeometries(GOutpFile):
    gausfile = open(GOutpFile, 'r')
    print("Reading " + GOutpFile)
    GOutp = gausfile.readlines()
    gausfile.close()

    index = 0
    atoms = []
    coords = []
    GeomIndexes = []

    # Find the geometry section and charge section
    for index in range(len(GOutp)):
        if ('Input orientation:' in GOutp[index]) or ("Standard orientation:" in GOutp[index]):
            GeomIndexes.append(index + 5)

    for i, gindex in enumerate(GeomIndexes):
        # Read coords
        coords.append([])
        for line in GOutp[gindex:]:
            if '--------------' in line:
                break
            else:
                data = [_f for _f in line[:-1].split(' ') if _f]
                if i == 0:
                    atoms.append(GetAtomSymbol(int(data[1])))
                coords[-1].append([float(x) for x in data[3:]])

    return atoms, coords


def CalcDist(geom, atoms):

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


def GetAtomSymbol(AtomNum):
    Lookup = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', \
              'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', \
              'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
              'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', \
              'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', \
              'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
              'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']

    if AtomNum > 0 and AtomNum < len(Lookup):
        return Lookup[AtomNum-1]
    else:
        print("No such element with atomic number " + str(AtomNum))
        return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for checking internal coordinates in a' +\
                                                 'molecule from a Gaussian output file')
    parser.add_argument('GausFile', help="Gaussian output file")
    parser.add_argument('IntCoords', help="Comma seperated list of 2 to 4 atom numbers, that define" +\
                        " the coordinate of interest. 2 atoms define distance, 3 - angle and" + \
                        " 4 - dihedral angle. Atom numbering starts with 1.")
    args = parser.parse_args()
    print(args.GausFile)
    print(args.IntCoords)

    StrIntCoords = args.IntCoords.split(',')
    IntCoords = []
    for x in StrIntCoords:
        try:
            IntCoords.append(int(x))
        except ValueError:
            print("Atom numbers must be integers.")
            quit()

    if (len(IntCoords) > 4) or (len(IntCoords) < 2):
        print("Invalid number of atoms, must be between 2 and 4")
        quit()

    if not os.path.isfile(args.GausFile):
        print(args.GausFile + " is not a valid file.")
        quit()

    main(args.GausFile, IntCoords)
