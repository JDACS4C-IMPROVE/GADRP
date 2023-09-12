import torch
import candle
import os
import cell_ae
import train

required = [
    'epochs'
]
additional_definitions = [
# {   
#     'name':'pool',
#     'nargs':'+',
#     'type': int,
#     'help':'network structure of shared layer'
# },
]

class BenchmarkGADRP(candle.Benchmark):
    """ Benchmark for GADRP. """

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
    train.train(args)
    return 0

def initialize_parameters():
    """ Initialize the parameters for the GADRP benchmark. """
    print("Initializing parameters\n")
    benchmark = BenchmarkGADRP(
        filepath=os.path.dirname(__file__),
        defmodel="gadrp_default_model.txt",
        framework="pytorch",
        prog="GADRP",
        desc="CANDLE compliant GADRP",
    )
    gParameters = candle.finalize_parameters(benchmark)
    return gParameters


def main():
    gParameters = initialize_parameters()
    print(gParameters)

    scores = run(gParameters)
    print("Done.")

if __name__ == "__main__":
    main()