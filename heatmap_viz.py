#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 22:23:48 2024

@author: ericsson
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Visualization Functions
def plot_heatmap(matrix, title, output_path, vmin=None, vmax=None):
    """Plot a heatmap of a matrix."""
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        matrix,
        cmap="coolwarm",
        vmin=vmin,
        vmax=vmax,
        cbar=True,
        annot=True,  # Display the values
        fmt=".2f",  # Format the numbers to 2 decimal places
        annot_kws={"size": 8},  # Size of the annotations
        linewidths=0.5  # Linewidth between cells
    )
    
    # Adjust ticks to start from 1 instead of 0
    num_rois = matrix.shape[0]  # Number of ROIs (matrix size)
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
    
    plt.title(title)
    plt.xlabel("ROI")
    plt.ylabel("ROI")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Heatmap saved to {output_path}")


def plot_distribution(data, title, output_path):
    """Plot a distribution of Z-values or T-values."""
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=data)
    plt.title(title)
    plt.xlabel("ROI Pair")
    plt.ylabel("Value")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Distribution plot saved to {output_path}")

# Set the output directory
output_dir = "out_yeo"
runs = ["run-01", "run-03"]

for run in runs:
    # Check for the group analysis results for this run
    t_values_file = os.path.join(output_dir, f"group_t_values_{run}_lsd_vs_plcb.csv")
    p_values_file = os.path.join(output_dir, f"group_p_values_{run}_lsd_vs_plcb.csv")

    if os.path.exists(t_values_file) and os.path.exists(p_values_file):
        # Load T-values and P-values
        t_values = np.loadtxt(t_values_file, delimiter=",")
        p_values = np.loadtxt(p_values_file, delimiter=",")

        # Plot heatmaps for T-values and P-values
        plot_heatmap(
            t_values,
            f"Group-Level T-Values (Paired Test, {run})",
            os.path.join(output_dir, f"group_t_values_heatmap_{run}.png"),
            vmin=-5, vmax=5
        )

        plot_heatmap(
            p_values,
            f"Group-Level P-Values (Paired Test, {run})",
            os.path.join(output_dir, f"group_p_values_heatmap_{run}.png"),
            vmin=0, vmax=0.05
        )

    else:
        print(f"T-values or P-values files not found for {run}. Ensure the main script has been executed correctly.")
