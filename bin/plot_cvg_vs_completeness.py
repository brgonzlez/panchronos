#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno

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
normalised_boundary_str = sys.argv[5]

print(f"Data file: {data}")
print(f"Gene Completeness: {completeness_str}")
print(f"Normalised Coverage Lower Bound: {coverage_lower_str}")
print(f"Normalised Coverage Upper Bound: {coverage_upper_str}")
print(f"Normalised Coverage Plot Bound: {normalised_boundary_str}")


try:
    completeness = int(completeness_str)
    coverage_lower = float(coverage_lower_str)
    coverage_upper = float(coverage_upper_str)
    normalised_boundary = float(normalised_boundary_str)
except ValueError as e:
    print(f"Error converting arguments to float: {e}")
    sys.exit(1)


with open(data, 'r') as file:
        summary = pd.read_csv(data, sep='\t')
        print(summary.head())

max_value = summary['normalizedGeneSimple'].max()

line_segments = [
    [(completeness, 0), (completeness, coverage_upper)],  # Line 1: vertical line 
    [(0, coverage_lower), (100, coverage_lower)], # Line 2: horizontal line for lower boundary
    [(0, coverage_upper), (100, coverage_upper)] # Line 3: horizontal line for upper boundary
]

high_genes_count = (summary['normalizedGeneSimple'] > normalised_boundary).sum()
plot_summary = summary[summary['normalizedGeneSimple'] <= normalised_boundary]


plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern
plt.figure(figsize=[15,10])
ax = sns.scatterplot(data=plot_summary, y='normalizedGeneSimple' , x='geneCompleteness' , hue='sampleID', s=20)

lines = LineCollection(line_segments, linewidths=2, colors='black', linestyle='dashed')
ax.add_collection(lines)

# Adjust y-axis based on max_value and coverage thresholds
if max_value > coverage_upper:
    new_ymax = max_value
else:
    new_ymax = coverage_upper + 1


# Create the Rectangle patch
rec1 = Rectangle((0, 0), completeness, coverage_upper, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3)  # completeness box
rec2 = Rectangle((completeness, 0), 100 - completeness, coverage_lower, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3) # lower box
rec3 = Rectangle((0, coverage_upper), 100, new_ymax, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3) # upper box

# Add the rectangle to the plot 
ax.add_patch(rec1)
ax.add_patch(rec2)
ax.add_patch(rec3)

plt.legend(markerscale=1, fontsize=10, loc='center right', bbox_to_anchor=(1.2,0.969))

ax.set_ylim(0, new_ymax)

# Force the plot to start at 0
ax.set_xlim(left=0)

ax.set_xlabel("Gene Completeness", fontsize=14)
ax.set_ylabel("Normalized Gene Coverage", fontsize=14)
ax.set_title("Normalised Coverage per Gene vs Gene Completeness", fontsize=16, pad=20)

sample_id = summary['sampleID'].iloc[0]
plt.savefig(f'plotCoverage_vs_Completeness_{sample_id}.png', bbox_inches='tight')
