This directory contains the following files:
- `hash.m` contains functions needed to run the hash function **KuHash**, and is the main file used to get the timings for Section 6. 
- `isogeny33.m` contains functions needed to compute the $(3,3)$-isogeny formulae described in Section 4 of the paper. 
- `kummer_arithmetic.m` contains functions needed to perform arithmetic on fast Kummer surfaces.
- `params128.m` contains the precomputed parameters for security level $\lambda = 128$.
- `params192.m` contains the precomputed parameters for security level $\lambda = 192$.
- `params256.m` contains the precomputed parameters for security level $\lambda = 256$.

This directory also contains a sub-directory `section4`, which contains files that verify claims made in Section 4 of the paper. Namely:
- `lemma-4_1.m` verifies that our $(3,3)$-isogenies are of the form as in the statement of the lemma. More specifically, it verifies that the invariance under action by $2$-torsion points results in the correct form for the $(3,3)$-isogeny.
- `linear-transform.m` verifies that the formulae for the intersection is correct, and verifies that applying the linear transformation to this intersection will give a (3,3)-isogeny between fast Kummer surfaces.