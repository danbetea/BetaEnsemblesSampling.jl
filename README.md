# BetaEnsemblesSampling

## Introduction

This package provides methods to sample iid (independent and identically distributed) eigenvalues from the classical beta ensembles of random matrix theory (RMT). So far the following multivariate eigenvalue distributions are supported:

- Gaussian Beta Ensemble with probability distribution

$$P(\lambda_1, \dots, \lambda_N)d \lambda_1 \dots d \lambda_N \propto \prod_{1 \leq i < j \leq N} |\lambda_i - \lambda_j|^\beta \prod_{1 \leq i \leq N} e^{-\frac{(\lambda_i - \mu)^2}{2 \sigma^2}}  d \lambda_i$$

- Laguerre Beta Ensemble with probability distribution

$$P(\lambda_1, \dots, \lambda_N)d \lambda_1 \dots d \lambda_N \propto \prod_{1 \leq i < j \leq N} |\lambda_i - \lambda_j|^\beta \prod_{1 \leq i \leq N} \lambda_i^{\alpha-1} e^{-\lambda_i} d \lambda_i$$

- Jacobi Beta Ensemble with probability distribution

$$P(\lambda_1, \dots, \lambda_N)d \lambda_1 \dots d \lambda_N \propto \prod_{1 \leq i < j \leq N} |\lambda_i - \lambda_j|^\beta \prod_{1 \leq i \leq N} \lambda_i^{a-1} (1-\lambda_i)^{b-1}  d \lambda_i$$

For generic $beta$ these distributions are abbreviated G$\beta$E, L$\beta$E, and J$\beta$E respectively. For $\beta = 1, 2, 4$ they are known as G/L/JOE, G/L/JUE, G/L/JSE, where O stantds for orthogonal, U for unitary, and S for symplectic.

## Method

It uses the tridiagonal method of Dumitriu--Edelman and Killip--Nenciu, following the paper by Gautier--Bardenet--Valko (namely the conventions and notation of Theorems 2.1, 2.2, and 2.3 of the latter). See the references for links.

## Dependencies

It uses ```LinearAlgebra```, ```Random```, and ```Distributions``` 

## References

- Dumitriu, Edelman, *Matrix models for beta ensembles*, [arXiv link](https://arxiv.org/pdf/math-ph/0206043.pdf)
- Gautier, Bardenet, Valko, *Fast sampling from beta ensembles*, [arXiv link](https://arxiv.org/pdf/2003.02344.pdf)
- Killip, Nenciu, *Matrix models for circular ensembles*, [arXiv link](https://arxiv.org/pdf/math/0410034.pdf)


