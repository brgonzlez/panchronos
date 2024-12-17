#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 13:58:15 2024

@author: bruno
"""

import sys
import pandas as pd
import os


# I need to make sure that every value in sys.argv[2] , sys.argv[3] and sys.argv[4] is a float variable.

data = sys.argv[1]
geneCompleteness = float(sys.argv[2])
lowerBound = float(sys.argv[3])
upperBound = float(sys.argv[4])


base_filename = os.path.basename(data)
output_file = f"{base_filename}_final.csv"

with open(data, 'r') as file:
        summary = pd.read_csv(data, sep='\t')


updated_series = []

for index, row in summary.iterrows():
    if row['completeness'] >= geneCompleteness and row['normalizedCoverage'] >= lowerBound and row['normalizedCoverage'] <= upperBound:
        updated_series.append(1)
    else:
        updated_series.append(0)


final_df = pd.DataFrame({'Gene': summary['Gene'], 'presenceAbsence': updated_series})
final_df.to_csv(output_file, index=False)
