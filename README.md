## Code to accompany: [*Information exclusion in quantum measurements and uncertainty relations*](https://arxiv.org/abs/2210.00958)
**Shan Huang, Hua-Lei Yin, Zeng-Bing Chen, and Shengjun Wu**

This is a repository for MATLAB code which was written for the article *Information exclusion in quantum measurements and uncertainty relations. Shan Huang, Hua-Lei Yin, Zeng-Bing Chen, and Shengjun Wu. [arXiv:2210.00958 [quant-ph]](https://arxiv.org/abs/2210.00958).*

- [Inf_Bound_Unitary](https://github.com/HS0Git/Information-exclusion/blob/main/Inf_Bound_Unitary.m) - computes the sum of information gain of different finite-dimensional quantum systems over arbitrary set of orthonormal bases (unitaries), and makes comparisons with the upper bound given in Eq. (8) of the main text. This code requires the free MATLAB toolbox [QETLAB](https://qetlab.com/) to generate random unitaries.
- [mub](https://github.com/HS0Git/Information-exclusion/blob/main/mub.m)-a function outputs MUBs of dimensions ranging from 2 to 5 and a collection of traceless orthogonal Hermitian operators (placed in a matrix) that can be used to expand all the traceless Hermitian operators of dimensions up to 5.
- [Inf_Bound_MUB](https://github.com/HS0Git/Information-exclusion/blob/main/Inf_Bound_MUB.m)-Computes the sum of information gain of low-dimensional quantum systems over MUBs or other orthonormal bases (unitaries), and makes comparisons with the upper bound given in Eq. (8) of the main text. One can use this code to generate not only random unitaries and MUBs, but also approximate MUBs.
- [entropybound](https://github.com/HS0Git/Information-exclusion/blob/main/entropybound.m)-compares Eqs. (14,16) with the respective numerical optimal state-independent entropic lower bounds.
- [SentropyIC](https://github.com/HS0Git/Information-exclusion/blob/main/SentropyIC.m)-a function outputs the optimal Shannon entropic lower bound that can be obtained based on the sum of indexes of coincidence (IC) only.
