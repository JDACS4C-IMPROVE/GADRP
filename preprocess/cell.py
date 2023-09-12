import numpy as np

from preprocess.utils import calculate_similarity
from preprocess.candle_original_defaults import default_args

def preprocess_cell(args):
    """
    caculate cell line similarity matrix
    :return:
    """
    calculate_similarity(id_file=args.cell_id_file,
                         input_files=args.cell_data_files,
                         output_file=args.cell_sim_file,
                         output_top10_file=args.cell_sim_top10_file)


def main():
    """Run preprocessing with default values"""

    preprocess_cell(args = default_args)


if __name__ == '__main__':
    main()
