import candle

default_args = candle.ArgumentStruct(**{
    # Cell
    'cell_id_file': '../data/cell_line/cell_index.csv',
    'cell_data_files': [
        '../data/cell_line/miRNA_470cell_734dim_clean.csv',
        '../data/cell_line/CpG_407cell_69641dim_clean.csv'
    ],
    'cell_sim_file': '../data/cell_line/cell_sim.pt',
    'cell_sim_top10_file': '../data/cell_line/cell_sim_top10.pt',

    # Drug
    'drug_id_file': '../data/drug/drug_index.csv',
    'drug_data_files': [
        '../data/drug/269_dim_physicochemical.csv',
    ],
    'drug_sim_file': '../data/drug/drug_sim.pt',
    'drug_sim_top10_file': '../data/drug/drug_sim_top10.pt',

    # Drug-cell
    'cell_index_file': '../data/cell_line/cell_index.csv',

    'edge_idx_file': '../data/pair/edge_idx_file.pt',
    'drug_cell_sim_file': '../data/pair/drug_cell_sim.pt',
    'drug_cell_sim_matrix_txt': '../data/pair/drug_cell_sim_matrix.txt',
    'drug_cell_sim_top100_index_file': '../data/pair/drug_cell_sim_top100_index.pt',

    'drug_cell_label_file': '../data/pair/drug_cell_response.csv',
    'drug_cell_label_index_file': '../data/pair/drug_cell_label.pt',
})