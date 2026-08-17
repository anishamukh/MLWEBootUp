# MLWEBootUp
This repository contains the proof-of-concept implementation of the following approaches for MLWE-CKKS:
1. Rank-adaptive encryption framework
2. Baseline bootstrapping (port from RLWE-CKKS to MLWE-CKKS domain)
3. Cross-rank optimization for the first bootstrapping after client upload


## Description

## How to build/run in SageMath
Different subroutines have been added as separate files with .sage extension, for example, the key generation sub-routine is under the file name *KeyGen.sage*. The file *example.sage* runs all the sub-routines for two plaintext messages (real numbers).
1. In order to run the entire scheme in a Jupyter notebook using a SageMath kernel, simply run the command, *load("example.sage")* in the notebook.
2. In order to run selected subroutines of the scheme in a Jupyter notebook using SageMath,  comment out the excess file names from the first 6 lines of *example.sage* following this order: *constants.sage, auxiliary_functions.sage, KeyGen.sage, Encrypt-Decrypt.sage, HMult.sage, Relin.sage*. For e.g., testing a simple encryption-decryption requires the files *constants.sage, auxiliary_functions.sage, KeyGen.sage, Encrypt-Decrypt.sage*, so *HMult.sage, Relin.sage* can be commented out.
**Note:** An online SageMath compiler, such as, https://sagecell.sagemath.org/, would tolerate parameter sizes upto N = 2^8 (change it in the parameter *POLDEG* in *constants.sage*) but might not run for higher parameters.

## How to build/run in OpenFHE

The README for the OpenFHE implementation can be found [here](openfhe-development/README.md).

## Acknowledgements

We express our gratitude towards Matthias Bergmann for contributing and integrating ModHE subroutines in OpenFHE. 
