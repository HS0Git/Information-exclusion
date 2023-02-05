## Code to accompany: [*Information exclusion in quantum measurements and uncertainty relations*](https://arxiv.org/abs/2210.00958)
**Shan Huang, Hua-Lei Yin, Zeng-Bing Chen, and Shengjun Wu**

This is a repository for MATLAB code which was written for the e-Print article *Information exclusion in quantum measurements and uncertainty relations. Shan Huang, Hua-Lei Yin, Zeng-Bing Chen, and Shengjun Wu. [arXiv:2210.00958 [quant-ph]](https://arxiv.org/abs/2210.00958).*

- [Inf_Bound_Unitary](https://github.com/HS0Git/Information-exclusion/blob/main/Inf_Bound_Unitary.m) - computes the sum of information gain of different finite-dimensional quantum systems over arbitrary set of orthonormal bases (unitaries), and makes comparisons with the upper bound given in Eq. (11) of the main text. This code requires the free MATLAB toolbox [QETLAB](https://qetlab.com/) to generate random unitaries.
- [mub](https://github.com/HS0Git/Information-exclusion/blob/main/mub.m)-a function outputs MUBs of dimensions ranging from 2 to 5 and a collection of traceless orthogonal Hermitian operators (placed in a matrix) that can be used to expand all the traceless Hermitian operators of dimensions up to 5.
- [Inf_Bound_MUB](https://github.com/HS0Git/Information-exclusion/blob/main/Inf_Bound_MUB.m)-Computes the sum of information gain of low-dimensional quantum systems over MUBs or other orthonormal bases (unitaries), and makes comparisons with the upper bound given in Eq. (11) of the main text. The unitaries are not generated at random, but denpend on how we generate the expansion parameters of Hermitian operators. Hence, one can use this code to generate MUBs, approximately MUBs by applying a unitary close to identity to MUBs, and random unitaries.
- [entropybound](https://github.com/HS0Git/Information-exclusion/blob/main/entropybound.m)-compares Eqs. (19,20) with the respective numerical optimal state-independent entropic lower bounds.
- [SentropyIC](https://github.com/HS0Git/Information-exclusion/blob/main/SentropyIC.m)-a function outputs the optimal Shannon entropic lower bound that can be obtained based on the sum of indexes of coincidence (IC) only, namely, the r.h.s. of Eq. (20) of the main text.
