# FDR-Estimation-for-Metabolomics
Methods for estimating FDR in metabolomics. The methods presented are based on a target-decoy search strategy.

The first two methods are based on generating decoys via a naive and spectrum-based method, as discussed [here](https://www.nature.com/articles/s41467-017-01318-5). The third method is based on generating "knockoffs" by first fitting a spec2ve model on our target spectra, then creating knockoffs from the embedding space via a Gaussian mixture model (as discussed [here](https://arxiv.org/abs/1807.06214)). 

Assessment of the FDR approaches is done by producing QQ-plots of estimated q-values against true q-values.


## Create environment with necessary packages
```
conda create --name FDR_Metabolomics python=3.7
conda activate FDR_Metabolomics
conda install --channel nlesc --channel bioconda --channel conda-forge matchms
conda install --channel nlesc --channel bioconda --channel conda-forge spec2vec
pip install sklearn
pip install jupyter
```
