#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:25:42 2024

@author: ericsson
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Visualization Functions
def plot_thresholded_heatmap(p_values, threshold, output_path):
    """Generate heatmap only showing connections (p <= threshold)."""
    # Mask connections above the threshold
    significant_mask = p_values <= threshold
    masked_p_values = np.where(significant_mask, p_values, np.nan)

    # Plotting
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        masked_p_values,
        cmap="viridis",
        annot=True,
        vmin=0,
        vmax=threshold
    )

    num_rois = p_values.shape[0]
    xticks_interval = 1
    yticks_interval = 1

    plt.xticks(
        ticks=np.arange(0, num_rois, xticks_interval) + 0.5,
        labels=[str(i + 1) for i in range(0, num_rois, xticks_interval)],
        fontsize=8
    )
    plt.yticks(
        ticks=np.arange(0, num_rois, yticks_interval) + 0.5,
        labels=[str(i + 1) for i in range(0, num_rois, yticks_interval)],
        rotation=0,
        fontsize=8
    )
    
    plt.title(f"Significant P-values (p <= {threshold})")
    plt.xlabel("ROI")
    plt.ylabel("ROI")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Thresholded heatmap saved to {output_path}")

# Set the output directory
output_dir = "out_yeo"
threshold = 0.05
runs = ["run-01", "run-03"]

for run in runs:
    # Check for the group analysis results for this run
    t_values_file = os.path.join(output_dir, f"group_t_values_{run}_lsd_vs_plcb.csv")
    p_values_file = os.path.join(output_dir, f"group_p_values_{run}_lsd_vs_plcb.csv")

    if os.path.exists(t_values_file) and os.path.exists(p_values_file):
        # Load T-values and P-values
        t_values = np.loadtxt(t_values_file, delimiter=",")
        p_values = np.loadtxt(p_values_file, delimiter=",")

        # Generate a thresholded heatmap for the current run
        plot_thresholded_heatmap(
            p_values,
            threshold=threshold,
            output_path=os.path.join(output_dir, f"thresholded_p_values_heatmap_{run}.png")
        )

    else:
        print(f"T-values or P-values files not found for {run}. Ensure the main script has been executed correctly.")
