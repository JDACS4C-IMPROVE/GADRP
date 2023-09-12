import csv
import os
import random

from torch import nn as nn
from model.autoencoder import Auto_Encoder
import torch
from sklearn.preprocessing import scale, MinMaxScaler
import pandas as pd
import torch.utils.data as Data
import datetime
from preprocess.candle_original_defaults import default_args

# cell
# os.environ["CUDA_VISIBLE_DEVICES"]="1"

def train_ae(model,trainLoader,test_feature,autoencoder_lr):
    start = datetime.datetime.now()
    optimizer = torch.optim.Adam(model.parameters(), lr=autoencoder_lr)
    loss_func = nn.MSELoss()
    best_model=model
    best_loss=100
    for epoch in range(1, 2500 + 1 ):
        for x in trainLoader:
            y=x
            encoded, decoded = model(x)
            train_loss = loss_func(decoded, y)
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()
        with torch.no_grad():
            y = test_feature
            encoded, decoded = model(test_feature)
            test_loss = loss_func(decoded, y)
        if (test_loss.item() < best_loss):
            best_loss = test_loss
            best_model = model
        if epoch%100==0:
            end = datetime.datetime.now()
            print('epoch:' ,epoch, 'train loss = ' ,train_loss.item(),"test loss:",test_loss.item(), "time:",(end - start).seconds)
    return best_model

def prepare_ae(args):
    random.seed(args.autoencoder_rng_seed)
    device = torch.device(args.device)
    cell_RNAseq_file = args.cell_rnaseq_file
    cell_copynumber_file = args.cell_copynumber_file
    cell_index_file = args.cell_id_file
    cell_RNAseq_ae = args.cell_rnaseq_ae
    cell_copynumber_ae = args.cell_copynumber_ae

    cell_index = pd.read_csv(cell_index_file, sep=',', header=None, index_col=[0])

    #normalization
    min_max = MinMaxScaler()

    # dimension reduction(gene expression data)
    if os.path.exists(cell_RNAseq_ae):
        print('RNAseq autoencoder already computed')
    else:
        print('Training RNAseq autoencoder')
        RNAseq_feature = pd.read_csv(cell_RNAseq_file, sep=',', header=None, index_col=[0], skiprows=2)
        RNAseq_feature = RNAseq_feature.loc[list(cell_index.index)].values
        RNAseq_feature = torch.tensor(min_max.fit_transform(RNAseq_feature)).float().to(device)
        RNAseq_indim = RNAseq_feature.shape[-1]
        print(RNAseq_indim)

        RNAseq_ae = Auto_Encoder(device,RNAseq_indim, 400)
        train_list = random.sample((RNAseq_feature).tolist(), int(0.9 * len(RNAseq_feature)))
        test_list = [item for item in (RNAseq_feature).tolist() if item not in train_list]
        train=torch.tensor(train_list).float().to(device)
        test = torch.tensor(test_list).float().to(device)
        data_iter = Data.DataLoader(train, args.autoencoder_batch_size, shuffle=True)
        best_model=train_ae(RNAseq_ae,data_iter,test,autoencoder_lr=args.autoencoder_lr)
        torch.save(best_model.output(RNAseq_feature),cell_RNAseq_ae)

    # dimension reduction(DNA copy number data)
    if os.path.exists(cell_copynumber_ae):
        print('Copy number autoencoder already trained')
    else:
        print('Training copy number autoencoder')
        copynumber_feature = pd.read_csv(cell_copynumber_file, sep=',', header=None, index_col=[0], skiprows=5)
        copynumber_feature = copynumber_feature.loc[list(cell_index.index)].values
        copynumber_feature = torch.tensor(min_max.fit_transform(copynumber_feature)).float().to(device)      
        copynumber_indim = copynumber_feature.shape[-1]
        print(copynumber_indim)

        copynumber_ae = Auto_Encoder(device,copynumber_indim, 400)
        train_list = random.sample((copynumber_feature).tolist(), int(0.9 * len(copynumber_feature)))
        test_list = [item for item in (copynumber_feature).tolist() if item not in train_list]
        train = torch.tensor(train_list).float().to(device)
        test = torch.tensor(test_list).float().to(device)
        data_iter = Data.DataLoader(train, args.autoencoder_batch_size, shuffle=True)
        best_model = train_ae(copynumber_ae, data_iter, test, autoencoder_lr=args.autoencoder_lr)
        torch.save(best_model.output(copynumber_feature), cell_copynumber_ae)


def main():
    prepare_ae(default_args)

if __name__ == '__main__':
    main()
