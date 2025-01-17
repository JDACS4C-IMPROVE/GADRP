from preprocess.utils import calculate_similarity
from preprocess.candle_original_defaults import default_args

drug_physicochemical_file = "../data/drug/269_dim_physicochemical.csv"
drug_sim_file= "../data/drug/drug_sim.pt"
drug_sim_top10_file="../data/drug/drug_sim_top10.pt"


def preprocess_drug(args):
    """
    caculate cell line similarity matrix
    :return:
    """
    calculate_similarity(id_file=args.drug_id_file,
                         input_files=args.drug_data_files,
                         output_file=args.drug_sim_file,
                         output_top10_file=args.drug_sim_top10_file)


def main():
    """Run preprocessing with default values"""

    preprocess_drug(args = default_args)


if __name__ == '__main__':
    main()
