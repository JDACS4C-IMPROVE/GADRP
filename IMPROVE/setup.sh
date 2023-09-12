set -e -u

source ~/miniconda3/etc/profile.d/conda.sh
conda create -n GADRP -y python==3.7.10 scikit-learn==0.24.2 numpy==1.19.2 pandas=1.3.4
conda activate GADRP
conda install -y pytorch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0 -c pytorch


TORCH=1.9.0
CUDA=cu102
pip install torch-scatter==2.0.9 -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
pip install torch-sparse==0.6.12 -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
pip install torch-cluster==1.5.9 -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
pip install torch-spline-conv==1.2.1 -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
pip install torch-geometric==2.3.1

( cd .. ; git lfs install ; git lfs fetch --all ; git lfs pull )

pip install git+https://github.com/ECP-CANDLE/candle_lib@0d215ec47d48191220a8cdd6fef3fbd00d4bd852

( cd .. ; git clone git@github.com:JDACS4C-IMPROVE/benchmark-dataset-generator.git ) # side by side with the project