#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno

Some GenBank files may have been submitted with wrong format. At the moment I don't have time to fix those so we will ignore them by adding to a blacklist.
"""

from Bio import SeqIO
import sys
import os

def verify_genbank_file(gb_file):
    try:
        with open(gb_file, "r") as handle:
            # trying to parse the file
            for _ in SeqIO.parse(handle, "genbank"):
                pass
        print(f"{gb_file} is a valid GenBank file.")
        return True
    except Exception as e:
        print(f"{gb_file} is not a valid GenBank file. Error: {e}")
        return False

workingDir = os.getcwd()

# Just verify them
for gb_file in os.listdir(workingDir):
    if gb_file.endswith(".gb"):
        verify_genbank_file(os.path.join(workingDir, gb_file)) 
