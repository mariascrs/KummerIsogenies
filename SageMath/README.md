This directory contains the following files:
- `fp_arith.py` contains code needed for operations in **F**<sub>p</sub>  and **F**<sub>p<sup>2</sup></sub>. This code was taken from the [ApresSQI codebase](https://github.com/TheSICQ/ApresSQI).
- `hash.py` contains the KummerHash class, which uses our (3,3)-isogeny formulae to construct the hash function **KuHash**.
- `isogeny33.py` contains functions needed to compute (3,3)-isogenies as given in Section 4 of the paper.
- `kummer_arithmetic.py` contains functions needed to perform arithmetic on fast Kummer surfaces.
- `main.py` is the main file used for benchmarking our hash function on various parameters.
- `param_load.py` contains functions that load the parameters needed to run the hash function **KuHash** at different security levels $\lambda = 128, 192, 256$.