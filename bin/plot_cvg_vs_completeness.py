#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno


TODO: 
    I NEED TO FIX BOX POSITION AND MOVE IT TO THE RIGHT, line 63
    I NEED TO INCORPORATE COMPLETENESS VALUES INTO THE GREY BOXES AND NOT HAVING 50 THERE
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle

data = sys.argv[1]
completeness_str = sys.argv[2]
coverage_lower_str = sys.argv[3]
coverage_upper_str = sys.argv[4]


print(f"Data file: {data}")
print(f"Completeness: {completeness_str}")
print(f"Coverage: {coverage_lower_str}")
print(f"Coverage: {coverage_upper_str}")

try:
    completeness = int(completeness_str)
    coverage_lower = float(coverage_lower_str)
    coverage_upper = float(coverage_upper_str)
except ValueError as e:
    print(f"Error converting arguments to float: {e}")
    sys.exit(1)


with open(data, 'r') as file:
        summary = pd.read_csv(data, sep='\t')
        print("geneNormalizedUpdated.tab file loaded successfully into python.")
        print(summary.head())

max_value = summary['normalizedGeneSimple'].max()

line_segments = [
    [(completeness, 0), (completeness, max_value)],  # Line 1: vertical line 
    [(0, coverage_lower), (100, coverage_lower)], # Line 2: horizontal line for lower boundary
    [(0, coverage_upper), (100, coverage_upper)]
]

plt.figure(figsize=[15,10])
ax = sns.scatterplot(data=summary, y='normalizedGeneSimple' , x='geneCompleteness' , hue='sampleID', s=12)

lines = LineCollection(line_segments, linewidths=2, colors='black', linestyle='dashed')
ax.add_collection(lines)


# Create the Rectangle patch
rec1 = Rectangle((0, 0), completeness, max_value, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3)
rec2 = Rectangle((completeness, 0), completeness, coverage_lower, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3)
rec3 = Rectangle((completeness, 0), completeness, coverage_upper, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3)

# Add the rectangle to the plot 
ax.add_patch(rec1)
ax.add_patch(rec2)
ax.add_patch(rec3)
plt.legend(markerscale=1, fontsize=10, loc='center right', bbox_to_anchor=(1.7,0.969))

sample_id = summary['sampleID'].iloc[0]
plt.savefig(f'plotCoverage_vs_Completeness_{sample_id}.png', bbox_inches='tight')
