# The-Bethe-Approximation-of-the-Pattern-Maximum-Likelihood-Distribution
This repository contains some MATLAB code, aiming to find the Bethe Approximation of the Pattern Maximum Likelihood distribution of a source when the alphabet is large and message is short.

## Content
Ensure all files are in the same directory. Some of the files are (maybe) not in this repo, please find them in my repo [Degree-M-Bethe-and-Sinkhorn-Permanent](https://github.com/Rmcher/Degree-M-Bethe-and-Sinkhorn-Permanent).

There are four different types of *computeBethePermanent* functions, with different suffix. *block* aims to make factor graphs's convergence faster, and *log* aims to avoid underflow of calculation, where some values are low to 1e-250. If you want to run the program directly and see the results, use *computeBethePermanent_double_block_log_all.m*.

The calculation logic and function are clearly stated. If not, please refer to [Degree-M-Bethe-and-Sinkhorn-Permanent](https://github.com/Rmcher/Degree-M-Bethe-and-Sinkhorn-Permanent).

## Copyright
Free to use.

Author: WU Binghong
Date: 2024.Nov.17
