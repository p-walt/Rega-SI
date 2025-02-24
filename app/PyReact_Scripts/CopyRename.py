#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyReact integrated workflow for the running of MM, DFT calculations
v0.1

Copyright (c) 2020 Kristaps Ermanis

Script that given a filename fragment, finds matching files and replaces the fragment
with the second filename fragment
"""


import argparse
import glob
import shutil


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for renaming files')
    parser.add_argument('FileNameFragment1', help="Filename fragment to be replaced without wildcards")
    parser.add_argument('FileNameFragment2', help="Filename fragment to replace the first one with")
    args = parser.parse_args()
    print(args.FileNameFragment1)
    print(args.FileNameFragment2)

    filelist = glob.glob('*' + args.FileNameFragment1 + '*')

    if len(filelist) == 0:
        print('Error: No matching files found, exiting...')
        quit()

    for fname in filelist:
        shutil.copyfile(fname, fname.replace(args.FileNameFragment1, args.FileNameFragment2))
