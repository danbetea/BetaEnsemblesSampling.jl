# BetaEnsemblesSampling

## Introduction

This Julia package provides methods to sample iid (independent and identically
distributed) eigenvalues from the classical $\beta$ (beta) ensembles of random
matrix theory (RMT). The parameter $\beta > 0$ plays the role of inverse
temperature. So far the following multivariate eigenvalue distributions on $N$
ordered tuples $(\lambda_1, \dots, \lambda_N)$ are supported:

- Gaussian Beta Ensemble with probability distribution on
$-\infty < \lambda_1 < \dots < \lambda_N < \infty$

$$P(\lambda_1, \dots, \lambda_N)d \lambda_1 \dots d \lambda_N \propto
\prod_{1 \leq i < j \leq N} |\lambda_i - \lambda_j|^\beta
\prod_{1 \leq i \leq N} e^{-\frac{(\lambda_i-\mu)^2}{2 \sigma^2}} d \lambda_i$$

with $\mu$ real and $\sigma, \beta > 0.$

- Laguerre Beta Ensemble with probability distribution on
$0 < \lambda_1 < \dots < \lambda_N < \infty$

$$P(\lambda_1, \dots, \lambda_N)d \lambda_1 \dots d \lambda_N \propto
\prod_{1 \leq i < j \leq N} |\lambda_i - \lambda_j|^\beta
\prod_{1 \leq i \leq N} \lambda_i^{\alpha-1} e^{-\lambda_i/\theta} d \lambda_i$$

with $\alpha, \theta, \beta > 0.$

- Jacobi Beta Ensemble with probability distribution on
$0 < \lambda_1 < \dots < \lambda_N < 1$

$$P(\lambda_1, \dots, \lambda_N)d \lambda_1 \dots d \lambda_N \propto
\prod_{1 \leq i < j \leq N} |\lambda_i - \lambda_j|^\beta
\prod_{1 \leq i \leq N} \lambda_i^{a-1} (1-\lambda_i)^{b-1}  d \lambda_i$$

  with $a, b, \beta > 0.$

For generic $\beta$ these distributions are abbreviated $G \beta E$,
$L \beta E$, and $J \beta E$ respectively. For $\beta = 1, 2, 4$ they are known
 as $G/L/JOE$, $G/L/JUE$, $G/L/JSE$, where O stands for orthogonal,
 U for unitary, and S for symplectic.

## Usage

To sample ```num_samples = 10000``` iid samples from the Laguerre Beta Ensemble
above for matrix size ```N = 1000```, and take only the smallest ordered
```num_taken = 100``` eigenvalues, for ```beta = 4```, ```alpha = 2.5```,
```theta = 1``` (from the LSE ensemble), you would do:

`Lambda = LbetaE_eval_samples(N, num_taken, alpha, theta, beta, num_samples)`

or

`Lambda = LbetaE_eval_samples(1000, 100, 2.5, 1.0, 4.0, 10000)`

Each of the ```num_samples``` rows of ```Lambda``` is an iid sample containing
the first 100 eigenvalues.

## Idea

While other packages providing similar functionality exist (see below), the
 idea is to provide a package that is:

- focused solely on sampling for researchers who want to experiment with
 classical random matrix models;
- fast;
- easy to use;
- easy to read;
- easy to modify for users unfamiliar with Julia;
- useful for batch sampling for situations in which taking a large number of iid
 samples (e.g. to compute expectations) might be needed.

## Method

The sampling is done using the tridiagonal method of Dumitriu--Edelman and
 Killip--Nenciu, following the paper by Gautier--Bardenet--Valko
 (namely the conventions and notation of Thm 2.1, 2.2, and 2.3 of the latter).
 See the references for links.

## Dependencies

It uses `LinearAlgebra`, `Random`, and `Distributions`.

## Other similar packages

The Python package [DPPy](https://github.com/guilgautier/DPPy) was originally
implemented by Gautier--Bardenet--Valko as a companion to their paper
referenced below.

The [RandomMatrices.jl](https://github.com/JuliaMath/RandomMatrices.jl)
implements some of the same functionality as well.

## References

### Relevant references for the algorithms used

- Dumitriu, Edelman, *Matrix models for beta ensembles*,
[arXiv link](https://arxiv.org/pdf/math-ph/0206043.pdf)
- Gautier, Bardenet, Valko, *Fast sampling from beta ensembles*,
[arXiv link](https://arxiv.org/pdf/2003.02344.pdf)
- Killip, Nenciu, *Matrix models for circular ensembles*,
[arXiv link](https://arxiv.org/pdf/math/0410034.pdf)

### Relevant references on random matrix theory

- Anderson, Guionnet, Zeitouni, *An introduction to random matrices*,
[link to PDF](https://www.wisdom.weizmann.ac.il/~zeitouni/cupbook.pdf)
- Edelman, *Eigenvalues and condition numbers of random matrices*, PhD thesis,
[link to PDF](https://math.mit.edu/~edelman/publications/eigenvalues_and_condition_numbers.pdf)
- Livan, Novaes, Vivo, *An introduction to random matrices - theory and practice*,
[arXiv link](https://arxiv.org/abs/1712.07903)
