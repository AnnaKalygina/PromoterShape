#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats
import argparse 

# def parse_shape_file(path):
#     print(f"Parsing shape file: {path}")
#     data = []
#     with open(path, 'r') as file:
#         raw_lines = file.readlines()
#     for line in raw_lines:
#         values = line.strip().split()
#         if len(values) > 401:
#             values = values[:401]
#         elif len(values) < 401:
#             values += [np.nan] * (401 - len(values))
#         data.append(values)
#     data = pd.DataFrame(data).apply(pd.to_numeric, errors='coerce')
#     print(f"Completed parsing {path}, shape: {data.shape}")
#     return data

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

def main():
    parser = argparse.ArgumentParser(prog='Bootstrap test')
    parser.add_argument("--output_path", help="Output path", required=True)
    parser.add_argument("--confidence_level", type=float, default=0.90, help="Confidence level for bootstrap test (default: 0.90)")
    parser.add_argument("--sample_size", type=float, default=100, help="Sample size percentage relative to data (default: 100%)")
    args = parser.parse_args()

    species = ["athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae"]
    properties = ["Buckle", "HelT", "MGW", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt"]

    input_dir = "/local/home/quee4387/shapes_csv/"

    n_bootstraps = 1000
    result_df = pd.DataFrame(columns=["species", "property", "position", "raw_values", "hmm_values", 
                                      "observed_difference", "alpha", "lower_bound", "upper_bound", "significant"])

    for spec in species:
        for prop in properties:
            print(f"\nProcessing species: {spec}, property: {prop}")
            
            path_to_raw = f"{input_dir}{spec}_{prop}_raw.csv"
            path_to_hmm = f"{input_dir}{spec}_{prop}_hmm.csv"

            try:
                raw_shapes = pd.read_csv(path_to_raw, header=None)
                hmm_shapes = pd.read_csv(path_to_hmm, header=None)
            except FileNotFoundError as e:
                print(f"File not found: {e}")
                continue

            n_samples = int(raw_shapes.shape[0] * (args.sample_size / 100))
            print(f"Sample size for {spec}, {prop}: {n_samples}")
            confidence_level = args.confidence_level
            print(f"Sample size is {n_samples} and a confidence level is {confidence_level}. Output is {args.output_path}.")

            temp_df = bootstrap_confidence_interval(spec, prop, raw_shapes, hmm_shapes, n_bootstraps, n_samples, confidence_level)
            print(temp_df)

            result_df = pd.concat([result_df, temp_df], ignore_index=True)
            print(f"Completed processing for {spec}, {prop}")

    result_df.to_csv(args.output_path, index=False)
    print(f"\nAll processing complete. Results saved to {args.output_path}")


if __name__ == "__main__":
    main()
