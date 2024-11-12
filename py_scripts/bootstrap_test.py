#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats

def parse_shape_file(path):
    print(f"Parsing shape file: {path}")
    data = []
    with open(path, 'r') as file:
        raw_lines = file.readlines()
    for line in raw_lines:
        values = line.strip().split()
        if len(values) > 401:
            values = values[:401]
        elif len(values) < 401:
            values += [np.nan] * (401 - len(values))
        data.append(values)
    data = pd.DataFrame(data).apply(pd.to_numeric, errors='coerce')
    print(f"Completed parsing {path}, shape: {data.shape}")
    return data

def bootstrap_confidence_interval(spec, prop, raw_shapes, hmm_shapes, n_bootstraps, n_samples, confidence_level=0.90):
    raw_shapes = stats.zscore(raw_shapes, axis = 1, nan_policy='omit')
    hmm_shapes = stats.zscore(hmm_shapes, axis = 1, nan_policy='omit')

    observed_diff = raw_shapes.mean(axis=0) - hmm_shapes.mean(axis=0)
    
    bootstrap_diffs = []
    for _ in range(n_bootstraps):
        if _ % 100 == 0:
            print(f"Bootstrap iteration {_}/{n_bootstraps}")
        raw_sample = raw_shapes.sample(n=n_samples, replace=True).mean(axis=0)
        hmm_sample = hmm_shapes.sample(n=n_samples, replace=True).mean(axis=0)
        
        bootstrap_diff = raw_sample - hmm_sample
        bootstrap_diffs.append(bootstrap_diff)
    
    
    bootstrap_diffs_df = pd.DataFrame(bootstrap_diffs)
    
    lower_bound = (1 - confidence_level) / 2
    upper_bound = 1 - lower_bound
    ci_lower = bootstrap_diffs_df.quantile(lower_bound, axis=0)
    ci_upper = bootstrap_diffs_df.quantile(upper_bound, axis=0)
    
    significant_positions = ~((ci_lower <= 0) & (ci_upper >= 0))

    
    return pd.DataFrame({
        "species": spec,
        "property": prop,
        "position": range(-200, 201),
        "raw_values": raw_shapes.mean(axis=0),
        "hmm_values": hmm_shapes.mean(axis=0),
        'observed_difference': observed_diff,
        "alpha": 1 - confidence_level,
        "lower_bound": ci_lower,
        "upper_bound": ci_upper,
        'significant': significant_positions
    })




species = ["athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae"]
properties = ["Buckle", "HelT", "MGW", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt"]

raw_seq_dir = "/local/home/quee4387/dna_shape/"
control_seq_dir = "/local/home/quee4387/dna_shape_control/"

n_bootstraps = 1000
n_samples = 100

result_df = pd.DataFrame(columns=["species", "property", "position", "raw_values", "hmm_values", 
                                  "observed_difference", "alpha", "lower_bound", "upper_bound", "significant"])

for spec in species:
    for prop in properties:
        print(f"\nProcessing species: {spec}, property: {prop}")
        path_to_raw = f"{raw_seq_dir}{spec}_{prop}_200.txt"
        path_to_hmm = f"{control_seq_dir}{spec}_{prop}_200_hmm.txt"
        confidence_level = 0.80
        
        try:
            raw_shapes = parse_shape_file(path_to_raw)
            hmm_shapes = parse_shape_file(path_to_hmm)
        except FileNotFoundError as e:
            print(f"File not found: {e}")
            continue
        
        temp_df = bootstrap_confidence_interval(spec, prop, raw_shapes, hmm_shapes, n_bootstraps, n_samples, confidence_level)
        print(temp_df)

        result_df = pd.concat([result_df, temp_df], ignore_index=True)
        print(f"Completed processing for {spec}, {prop}")

output_path = "/local/home/quee4387/raw_vs_hmm_bootstrap_test_cl_0.80.csv"
result_df.to_csv(output_path, index=False)
print(f"\nAll processing complete. Results saved to {output_path}")


