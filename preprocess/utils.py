import os
import candle
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import pearsonr
import torch
import numpy as np
from tqdm import tqdm


def calculate_similarity(id_file, input_files, output_file, output_top10_file):
    """
    Calculate similarities between different vectors
    
    Store the similarity matrix as well as a matrix of 10 closest entries for each one
    """

    if os.path.exists(output_file) and os.path.exists(output_top10_file):
        print(f'Preprocessed files already exist: \n\t{output_file}\n\t{output_top10_file}')
        return

    # List of cell line IDs that are common across all data files
    id_set = pd.read_csv(id_file, sep=',', header=None, index_col=[0])
    num = id_set.shape[0]

    # Resulting similarity matrix
    sim = torch.zeros(size=(num, num))

    for data_file in input_files:
        print(f'Processing {data_file}')

        #load microRNA expression data, DNA methylation data of cell line
        features_full = pd.read_csv(data_file, sep=',', header=0, index_col=[0])

        features = features_full.loc[list(id_set.index)].values


        # Normalization
        min_max=MinMaxScaler()
        features = min_max.fit_transform(features)


        #calculate similarity matrix
        features_sim = torch.zeros(size=(num, num))


        for i in tqdm(range(num)):
            fi = features[i, :]
            for j in range(num):
                temp_sim = pearsonr(fi, features[j, :])
                features_sim[i][j] = np.abs(temp_sim[0])

        sim = sim + features_sim

    #calculate cell line similarity matrix
    sim = sim / len(input_files)

    # Simplify to only top ten closest cell lines for each
    _,sim_top10 = torch.topk(sim, 10, dim=1)

    torch.save(sim, output_file)
    torch.save(sim_top10, output_top10_file)
