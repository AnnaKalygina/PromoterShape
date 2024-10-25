#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.stats as stats

'''
This script combines predicted shapes in a long data format.
'''
species = ["athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae"]
properties = ["Buckle", "HelT", "MGW", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt"]
sources = ["raw", "shuffled", "hmm"]

raw_seq_dir = "/local/home/quee4387/dna_shape"
control_seq_dir = "/local/home/quee4387/dna_shape_control"

result_df = pd.DataFrame(columns=["species", "property", "source", "position", "value"])

for spec in species:
    for prop in properties:
        for scr in sources:
            if scr == "raw":
                path = raw_seq_dir + "/" + spec + "_" + prop + "_200.txt"
            else:
                path = control_seq_dir + "/" + spec + "_" + prop + "_200_" + scr + ".txt"
            try: 
                print(f'Parsing {prop} for {spec} from {scr}.')

                # Handles every row individually
                with open(path, 'r') as file:
                    raw_lines = file.readlines()
                
                data = []
                for line in raw_lines:
                    values = line.strip().split()
                    if len(values) > 401:
                        values = values[:401] 
                    elif len(values) < 401:
                        values += [np.nan] * (401 - len(values)) 
                    data.append(values)

                data = pd.DataFrame(data).apply(pd.to_numeric, errors='coerce')
                
                # Z-score the data and take the average.
                data_zscored = stats.zscore(data, axis = 1, nan_policy='omit')
                data_averaged = data_zscored.mean(axis = 0)

                positions = list(range(-200, 201))
                temp_df = pd.DataFrame({
                        "species": spec,
                        "property": prop,
                        "source": scr,
                        "position": positions,
                        "value": data_averaged
                    })

                result_df = pd.concat([result_df, temp_df], ignore_index=True)

            except Exception as e:
                print(f"Error processing file {path}: {e}")
        
result_df.to_csv("combined_dna_shape_table.csv", index=False)