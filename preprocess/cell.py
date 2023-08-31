import os

import numpy as np

import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import pearsonr
import torch
import candle
from tqdm import tqdm

def preprocess_cell(args):
    """
    caculate cell line similarity matrix
    :return:
    """
    cell_id_file = args.cell_id_file
    cell_sim_file = args.cell_sim_file
    cell_sim_top10_file = args.cell_sim_top10_file

    # List of cell line IDs that are common across all data files
    cell_id_set = pd.read_csv(cell_id_file, sep=',', header=None, index_col=[0])
    cell_num = cell_id_set.shape[0]

    # Resulting similarity matrix
    cell_sim = torch.zeros(size=(cell_num, cell_num))

    for data_file in args.cell_data_files:
        print(f'Processing {data_file}')

        #load microRNA expression data, DNA methylation data of cell line
        features_full = pd.read_csv(data_file, sep=',', header=0, index_col=[0])

        features = features_full.loc[list(cell_id_set.index)].values


        # Normalization
        min_max=MinMaxScaler()
        features = min_max.fit_transform(features)


        #calculate similarity matrix
        features_sim = torch.zeros(size=(cell_num, cell_num))


        for i in tqdm(range(cell_num)):
            fi = features[i, :]
            for j in range(cell_num):
                temp_sim = pearsonr(fi, features[j, :])
                features_sim[i][j] = np.abs(temp_sim[0])

        cell_sim = cell_sim + features_sim

    #calculate cell line similarity matrix
    cell_sim = cell_sim / len(args.cell_data_files)

    # Simplify to only top ten closest cell lines for each
    _,cell_sim_top10 = torch.topk(cell_sim, 10, dim=1)

    torch.save(cell_sim, cell_sim_file)
    torch.save(cell_sim_top10, cell_sim_top10_file)


def main():
    """Run preprocessing with default values"""

    default_args = {
        'cell_id_file': "../data/cell_line/cell_index.csv",
        'cell_data_files': [
            '../data/cell_line/miRNA_470cell_734dim_clean.csv',
            '../data/cell_line/CpG_407cell_69641dim_clean.csv' # Simplifying loading, no header needed
        ],
        'cell_sim_file': "../data/cell_line/cell_sim.pt",
        'cell_sim_top10_file': "../data/cell_line/cell_sim_top10.pt",
    }

    preprocess_cell(args = candle.ArgumentStruct(**default_args))


if __name__ == '__main__':
    main()
