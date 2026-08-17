# MLWEBootUp
This repository contains the proof-of-concept implementation of the following approaches for MLWE-CKKS:
1. Rank-adaptive encryption framework
2. Baseline bootstrapping (port from RLWE-CKKS to MLWE-CKKS domain)
3. Cross-rank optimization for the first bootstrapping after client upload

# Requirements

1. The implementation requires SageMath.
2. The code has been developed and tested using SageMath version 10.3 in a Conda environment under WSL.

 
## How to build/run the files in SageMath
 The file names containing 'example' contain small examples for experimenting with the different constructions present in the POC. 
 
 For instance, to run a simple example for the baseline leveled ModHE scheme in a SageMath environment, run:

```text
 sage example.sage
```

Similarly, to test a simple bootstrapping along with the rank-adaptive encryption, run:

```text
 sage boost_modhe_heaanboot_example.sage
```
File names containing 'test' have been used to report experimental measurements in the paper.
