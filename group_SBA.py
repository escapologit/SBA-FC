#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 21:29:28 2024

@author: ericsson
"""

import nibabel as nib
import numpy as np
import os
from scipy.stats import ttest_rel
import nibabel.freesurfer as fs

# Define directories and parameters
data_dir = "mri_surf2surf"  # Base directory for subject files
output_dir = "out_yeo"  # Directory to save extracted time series
sessions = ["ses-LSD", "ses-PLCB"]
runs = ["run-01", "run-03"]
hemispheres = ["lh", "rh"]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load a GIFTI file
def load_gii(file_path):
    """Load a GIFTI file."""
    gii_file = nib.load(file_path)
    # Collect data from all darrays
    data_list = [darray.data for darray in gii_file.darrays]
    # Stack the data along the second axis (timepoints)
    data_combined = np.column_stack(data_list)
    return data_combined

# Load an annotation file (Yeo 2011 Atlas)
def load_annot(atlas_file):
    """Load an annotation file."""
    # Load the annotation file using nibabel's FreeSurfer reader
    annot = fs.read_annot(atlas_file)
    labels = annot[0]  # ROI labels
    ctab = annot[1]    # Color table (optional, not needed for this analysis)
    return labels

# Extract mean time series for each ROI
def extract_mean_timeseries(bold_data, atlas_data, roi_labels):
    """Extract mean time series for each ROI."""
    n_timepoints = bold_data.shape[-1]
    n_rois = len(roi_labels)
    time_series = np.zeros((n_rois, n_timepoints))

    for i, roi in enumerate(roi_labels):
        roi_mask = atlas_data == roi
        roi_mask = roi_mask.flatten()  # Ensure mask is 1D
        roi_bold = bold_data[roi_mask, :]
        time_series[i, :] = np.mean(roi_bold, axis=0)

    return time_series

# Perform pairwise Pearson correlation
def compute_connectivity(timeseries):
    """Compute pairwise Pearson correlation for ROI time series."""
    corr_matrix = np.corrcoef(timeseries)
    return corr_matrix

# Perform Fisher R-to-Z transformation
def fisher_r_to_z(corr_matrix):
    """Perform Fisher R-to-Z transformation."""
    epsilon = 1e-6
    corr_matrix = np.clip(corr_matrix, -1 + epsilon, 1 - epsilon)
    z_matrix = np.arctanh(corr_matrix)
    np.fill_diagonal(z_matrix, 0.0)  # Reset diagonal to 0
    return z_matrix

# Group-level analysis functions
def load_z_matrices_for_run(subjects, run, session, output_dir):
    """Load Z-matrices for all subjects for a given run and session."""
    z_matrices = []
    for subject in subjects:
        subject_file = os.path.join(output_dir, f"{subject}_{session}_{run}_fisher_z_matrix.csv")
        if os.path.exists(subject_file):
            z_matrix = np.loadtxt(subject_file, delimiter=",")
            z_matrices.append(z_matrix)
        else:
            print(f"Z-matrix file not found for {subject}, {session}, {run}")
    return np.array(z_matrices)

def group_analysis(z_matrices_session1, z_matrices_session2):
    """Perform paired t-tests for group-level analysis."""
    t_values, p_values = ttest_rel(z_matrices_session1, z_matrices_session2, axis=0)
    return t_values, p_values

# Subjects
subjects = ["sub-001", "sub-002", "sub-004", "sub-006", "sub-009", "sub-010", "sub-011", "sub-013", "sub-017", "sub-018", "sub-019", "sub-020"]

# Subject-level processing
for subject in subjects:
    for session in sessions:
        for run in runs:
            bold_data_hemis = []
            atlas_data_hemis = []

            for hemi in hemispheres:
                bold_file = os.path.join(
                    data_dir, subject, session, f"{subject}_{session}_{run}_fsaverage_{hemi}.gii"
                )
                atlas_file = f"Yeo_JNeurophysiol11_FreeSurfer/fsaverage/label/{hemi}.Yeo2011_7Networks_N1000.annot"

                if not os.path.exists(bold_file) or not os.path.exists(atlas_file):
                    continue

                bold_data = load_gii(bold_file)
                atlas_data = load_annot(atlas_file)  # Use the new load_annot function

                bold_data_hemis.append(bold_data)
                atlas_data_hemis.append(atlas_data)

            if len(bold_data_hemis) != 2 or len(atlas_data_hemis) != 2:
                continue

            bold_data_combined = np.concatenate(bold_data_hemis, axis=0)
            atlas_data_combined = np.concatenate(atlas_data_hemis, axis=0)

            roi_labels = np.unique(atlas_data_combined)
            roi_labels = roi_labels[roi_labels > 0]

            roi_timeseries = extract_mean_timeseries(bold_data_combined, atlas_data_combined, roi_labels)

            connectivity_matrix = compute_connectivity(roi_timeseries)
            fisher_z_matrix = fisher_r_to_z(connectivity_matrix)

            # Save Fisher Z-matrix for the current session and run
            output_file = os.path.join(output_dir, f"{subject}_{session}_{run}_fisher_z_matrix.csv")
            np.savetxt(output_file, fisher_z_matrix, delimiter=",")
            print(f"Fisher Z-matrix saved to {output_file}")

# Group-level analysis
for run in runs:
    z_matrices_lsd = load_z_matrices_for_run(subjects, run, "ses-LSD", output_dir)
    z_matrices_plcb = load_z_matrices_for_run(subjects, run, "ses-PLCB", output_dir)

    if z_matrices_lsd.size > 0 and z_matrices_plcb.size > 0:
        t_values, p_values = group_analysis(z_matrices_lsd, z_matrices_plcb)

        # Save results for this run
        t_output_file = os.path.join(output_dir, f"group_t_values_{run}_lsd_vs_plcb.csv")
        p_output_file = os.path.join(output_dir, f"group_p_values_{run}_lsd_vs_plcb.csv")
        np.savetxt(t_output_file, t_values, delimiter=",")
        np.savetxt(p_output_file, p_values, delimiter=",")
        print(f"Group-level analysis complete for {run}. T and P values saved.")
