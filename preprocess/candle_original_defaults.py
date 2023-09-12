import candle

default_args = candle.ArgumentStruct(**{
    # Autoencoder
    'autoencoder_rng_seed': 4,
    'autoencoder_lr': 0.0001,
    'autoencoder_batch_size': 388,

    # Main training
    'lr': 0.0001,
    'epochs': 150,
    'batch_size': 1024,

    # Cell inputs
    'cell_id_file': '../data/cell_line/cell_index.csv',
    'cell_data_files': [
        '../data/cell_line/miRNA_470cell_734dim_clean.csv',
        '../data/cell_line/CpG_407cell_69641dim_clean.csv'
    ],
    'cell_rnaseq_file': '../data/cell_line/RNAseq_462cell_48392dim.csv',
    'cell_copynumber_file': '../data/cell_line/copynumber_461cell_23316dim.csv',

    # Cell outputs
    'cell_sim_file': '../data/cell_line/cell_sim.pt',
    'cell_sim_top10_file': '../data/cell_line/cell_sim_top10.pt',
    'cell_rnaseq_ae': '../data/cell_line/cell_RNAseq400_ae.pt',
    'cell_copynumber_ae': '../data/cell_line/cell_copynumber400_ae.pt',

    # Drug inputs
    'drug_id_file': '../dcd ata/drug/drug_index.csv',
    'drug_data_files': [
        '../data/drug/269_dim_physicochemical.csv',
    ],
    'drug_fingerprint_file': '../data/drug/881_dim_fingerprint.csv',

    # Drug outputs
    'drug_sim_file': '../data/drug/drug_sim.pt',
    'drug_sim_top10_file': '../data/drug/drug_sim_top10.pt',

    # Drug-cell inputs
    'cell_index_file': '../data/cell_line/cell_index.csv',
    'drug_cell_sim_matrix_txt': '../data/pair/drug_cell_sim_matrix.txt',
    'drug_cell_label_file': '../data/pair/drug_cell_response.csv',

    # Drug-cell outputs
    'edge_idx_file': '../data/pair/edge_idx_file.pt',
    'drug_cell_sim_file': '../data/pair/drug_cell_sim.pt',
    'drug_cell_sim_top100_index_file': '../data/pair/drug_cell_sim_top100_index.pt',
    'drug_cell_label_index_file': '../data/pair/drug_cell_label.pt',
})