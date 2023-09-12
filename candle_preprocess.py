from preprocess.cell import preprocess_cell
from preprocess.drug import preprocess_drug
from preprocess.drug_cell import preprocess_drug_cell

import torch
import candle
import os

required = [
    'device',

    'cell_id_file',
    'cell_data_files',
    'cell_sim_file',
    'cell_sim_top10_file',

    'drug_id_file',
    'drug_data_files',
    'drug_sim_file',
    'drug_sim_top10_file',
]

additional_definitions = [
# {   
#     'name':'pool',
#     'nargs':'+',
#     'type': int,
#     'help':'network structure of shared layer'
# },
]

class PreprocessGADRP(candle.Benchmark):
    """ Preprocess GADRP data """

    def set_locals(self):
        """ Set parameters for the benchmark.

        Args:
            required: set of required parameters for the benchmark.
            additional_definitions: list of dictionaries describing the additional parameters for the
            benchmark.
        """
        if required is not None:
            self.required = set(required)
        if additional_definitions is not None:
            self.additional_definitions = additional_definitions

def run(gParameters):
    print("In Run Function:\n")
    args = candle.ArgumentStruct(**gParameters)

    preprocess_cell(args)
    preprocess_drug(args)
    preprocess_drug_cell(args)

    return 0

def initialize_parameters():
    """ Initialize the parameters for the GADRP preprocessing. """
    print("Initializing parameters\n")
    preprocess = PreprocessGADRP(
        filepath=os.path.dirname(__file__),
        defmodel="gadrp_default_model.txt",
        framework="pytorch",
        prog="GADRP",
        desc="CANDLE compliant GADRP",
    )
    gParameters = candle.finalize_parameters(preprocess)
    return gParameters


def main():
    gParameters = initialize_parameters()
    print(gParameters)

    run(gParameters)
    print("Done.")

if __name__ == "__main__":
    main()