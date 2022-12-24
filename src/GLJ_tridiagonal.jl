"""

tridiagonal sampling routines for Hermite, Laguerre, and Jacobi beta ensembles

uses ideas from the Edelman--Dumitriu (and Killip--Nenciu for Jacobi ensembles) 
tridiagonal method described in (e.g. for Laguerre) Theorem 2.2 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko 
with arXiv reference: arXiv:2003.02344 [stat.CO]

each function returns a matrix, each line in that matrix is an iid sample of 
num_taken eigenvalues:
    e_1, ..., e_{num_taken}

"""

using Random, Distributions, LinearAlgebra

""" 
samples num_samples instances of num_taken (out of total of N) Gaussian Beta Ensemble points 
(inverse temperature = power of Vandermonde interaction = beta) 
from the distribution with univariate potential 
    
        exp( -(x - mu)^2 / (2 sigma_sq) ) 

(mean mu, variance sigma_sq normal distribution)

uses the Edelman--Dumitriu tridiagonal method described in Theorem 2.1 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko arXiv:2003.02344 [stat.CO]

input: N = size of GbetaE matrix (it is N x N)
       num_taken = the number of smallest (sorted) eigenvalues taken out of the N total ones
       mu = mean of Gaussian
       sigma_sq = variance of Gaussian
       beta = inverse temperature
       num_samples = number of independent samples

output: Lambda = num_samples x num_taken (usually N = all), each of its columns an independent GbetaE sample
"""
function GbetaE_eval_samples(N::Int64, num_taken::Int64, mu::Float64, sigma_sq::Float64, beta::Float64, num_samples::Int64)::Array{Float64,2}
    A = zeros(N, num_samples)
    B = zeros(N-1, num_samples)

    # TODO: optimize code below
    for n in 1:N 
        A[n, :] = rand(Normal(mu, sigma_sq), num_samples)
    end

    for n in 1:(N-1)
        B[n, :] = rand(Gamma((N-n)*beta/2, sigma_sq), num_samples)
    end

    # each row of Lambda will be one vector of num_taken GUE points ordered ascendingly 
    Lambda = zeros(num_samples, num_taken) 

    for k in 1:num_samples
        J = SymTridiagonal(A[:, k], sqrt.(B[:, k])) # construct the J(acobi) tridiagonal matrix
        Lambda[k, :] = eigvals(J, range(1, length=num_taken)) # TODO: bound checks on num_taken
    end

    return Lambda
end


"""
samples num_samples instances of num_taken (out of total of N) Gaussian Unitary Ensemble (GUE) points 
(inverse temperature = power of Vandermonde interaction = beta) 
from the distribution with univariate potential 
    
        exp( -x^2 / 2 ) 

(mean mu, variance sigma_sq normal distribution)

uses the Edelman--Dumitriu tridiagonal method described in Theorem 2.1 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko arXiv:2003.02344 [stat.CO]

input: N = size of GUE matrix (it is N x N)
       num_taken = the number of smallest (sorted) eigenvalues taken out of the N total ones
       num_samples = number of independent samples

output: Lambda = num_samples x num_taken (usually N = all), each of its columns an independent GUE sample
"""
function GUE_eval_samples(N::Int64, num_taken::Int64, num_samples::Int64)::Array{Float64,2}
    A = zeros(N, num_samples)
    B = zeros(N-1, num_samples)

    # TODO: optimize code below
    for n in 1:N 
        A[n, :] = rand(Normal(0.0, 1.0), num_samples)
    end

    for n in 1:(N-1)
        B[n, :] = rand(Gamma((N-n)*beta/2, sigma_sq), num_samples)
    end

    # each row of Lambda will be one vector of num_taken GUE points ordered ascendingly 
    Lambda = zeros(num_samples, num_taken) 

    for k in 1:num_samples
        J = SymTridiagonal(A[:, k], sqrt.(B[:, k])) # construct the J(acobi) tridiagonal matrix
        Lambda[k, :] = eigvals(J, range(1, length=num_taken)) # TODO: bound checks on num_taken
    end

    return Lambda
end


"""
samples num_samples instances of num_taken (out of total of N) Laguerre Beta Ensemble points 
(inverse temperature = power of Vandermonde interaction = beta) 
from the distribution with univariate potential 
    
        x^( alpha - 1 ) * exp( -x / theta ) --- NOTE the -1

uses the Edelman--Dumitriu tridiagonal method described in Theorem 2.2 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko arXiv:2003.02344 [stat.CO]

input: N = size of LbetaE matrix (it is N x N)
        num_taken = the number of smallest (sorted) eigenvalues taken out of the N total ones
        alpha = Laguerre parameter (technically, see above, it's alpha-1), also shape parameter for the gamma distribution
        theta = scale parameter
        beta = inverse temperature = power of the Vandermonde
        num_samples = number of independent samples

output: num_samples x num_taken (usually N = all), each of its columns an independent LbetaE sample
"""
function LbetaE_eval_samples(N::Int64, num_taken::Int64, alpha::Float64, theta::Float64, beta::Float64, num_samples::Int64)::Array{Float64,2}
    Xi = zeros(2*N-1, num_samples) 
    A = zeros(N, num_samples)
    B = zeros(N-1, num_samples)

    # each column of Xi corresponds to idependent vectors of independent entries xi[n] as in the ref. paper 
    for n in 1:(N-1) 
        Xi[2*n-1, :] = rand(Gamma((N-n)*beta/2 + alpha, theta), num_samples)
        Xi[2*n, :] = rand(Gamma((N-n)*beta/2, theta), num_samples)
    end
    Xi[2*N-1, :] = rand(Gamma(alpha, theta), num_samples)

    # TODO: optimize code below
    A[1, :] = Xi[1, :]
    for n in 2:N 
        A[n, :] = Xi[2*n-2, :] + Xi[2*n-1, :]
    end

    for n in 1:(N-1)
        B[n, :] = Xi[2*n-1, :] .* Xi[2*n, :]
    end

    # each row of Lambda will be one vector of num_taken GUE points ordered ascendingly 
    Lambda = zeros(num_samples, num_taken) 

    for k in 1:num_samples
        J = SymTridiagonal(A[:, k], sqrt.(B[:, k])) # construct the J(acobi) tridiagonal matrix
        Lambda[k, :] = eigvals(J, range(1, length=num_taken)) # TODO: bound checks on num_taken
    end

    return Lambda
end


"""
samples num_samples instances of num_taken (out of max N) points from the N x N LUE 
(Laguerre Unitary Ensemble, beta = 2) having distribution with univariate potential 

    x^(alpha-1) * exp(-x) --- NOTE the -1

uses the Edelman--Dumitriu tridiagonal method described in Theorem 2.2 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko arXiv:2003.02344 [stat.CO]

input: N = size of LUE matrix (it is N x N)
       num_taken = the number of smallest (sorted) eigenvalues taken out of the N total ones
       alpha = Laguerre parameter (technically, see above, it's alpha-1)
       num_samples = number of independent samples

output: Lambda = num_samples x num_taken (usually N = all), each of its columns an independent LUE sample
"""
function LUE_eval_samples(N::Int64, num_taken::Int64, alpha::Float64, num_samples::Int64)::Array{Float64,2}
    Xi = zeros(2*N-1, num_samples) 
    A = zeros(N, num_samples)
    B = zeros(N-1, num_samples)

    # each column of Xi corresponds to idependent vectors of independent entries xi[n] as in the ref. paper 
    for n in 1:(N-1) 
        Xi[2*n-1, :] = rand(Gamma(N-n+alpha, 1.0), num_samples)
        Xi[2*n, :] = rand(Gamma(N-n, 1.0), num_samples)
    end
    Xi[2*N-1, :] = rand(Gamma(alpha, 1.0), num_samples)

    # TODO: optimize code below
    A[1, :] = Xi[1, :]
    for n in 2:N 
        A[n, :] = Xi[2*n-2, :] + Xi[2*n-1, :]
    end

    for n in 1:(N-1)
        B[n, :] = Xi[2*n-1, :] .* Xi[2*n, :]
    end

    # each row of Lambda will be one vector of num_taken GUE points ordered ascendingly 
    Lambda = zeros(num_samples, num_taken) 

    for k in 1:num_samples
        J = SymTridiagonal(A[:, k], sqrt.(B[:, k])) # construct the J(acobi) tridiagonal matrix
        Lambda[k, :] = eigvals(J, range(1, length=num_taken)) # TODO: bound checks on num_taken
    end

    return Lambda
end


"""
samples num_samples instances of num_taken (out of total of N) Jacobi Beta Ensemble points 
(inverse temperature = power of Vandermonde interaction = beta) 
from the distribution with univariate potential 
    
        x^( a - 1 ) * ( 1 - x )^( b - 1 ) --- NOTE the -1 in the exponents

uses the Killip--Nenciu tridiagonal method described in Theorem 2.3 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko arXiv:2003.02344 [stat.CO]

input: N = size of JbetaE matrix (it is N x N)
       num_taken = the number of smallest (sorted) eigenvalues taken out of the N total ones
       a = first Jacobi parameter (technically, see above, it's a-1)
       b = second Jacobi parameter (technically, see above, it's b-1)
       beta = inverse temperature = power of the Vandermonde
       num_samples = number of independent samples

output: Lambda = num_samples x num_taken (usually N = all), each of its columns an independent JbetaE sample
"""
function JbetaE_eval_samples(N::Int64, num_taken::Int64, a::Float64, b::Float64, beta::Float64, num_samples::Int64)::Array{Float64,2}
    C = zeros(2*N-1, num_samples) 
    A = zeros(N, num_samples)
    B = zeros(N-1, num_samples)
    unos = ones(num_samples)

    # each column of C corresponds to idependent vectors of independent entries c[n] as in the ref. paper 
    for n in 1:(N-1) 
        C[2*n-1, :] = rand(Beta((N-n)*beta/2 + a, (N-n)*beta/2 + b), num_samples)
        C[2*n, :] = rand(Beta((N-n)*beta/2, (N-n-1)*beta/2 + a + b), num_samples)
    end
    C[2*N-1, :] = rand(Beta(a, b), num_samples)

    # TODO: optimize code below
    A[1, :] = C[1, :]
    for n in 2:N 
        A[n, :] = (unos - C[2*n-3, :]) .* C[2*n-2, :] + (unos - C[2*n-2, :]) .* C[2*n-1, :]
    end

    B[1, :] = C[1, :] .* (unos - C[1, :]) .* C[2, :]
    for n in 1:(N-1)
        B[n, :] = (unos - C[2*n-2, :]) .* C[2*n-1, :] .* (unos - C[2*n-1, :]) .* C[2*n, :]
    end

    # each row of Lambda will be one vector of num_taken GUE points ordered ascendingly 
    Lambda = zeros(num_samples, num_taken) 

    for k in 1:num_samples
        J = SymTridiagonal(A[:, k], sqrt.(B[:, k])) # construct the J(acobi) tridiagonal matrix
        Lambda[k, :] = eigvals(J, range(1, length=num_taken)) # TODO: bound checks on num_taken
    end

    return Lambda
end


"""
samples num_samples instances of num_taken (out of total of N) Jacobi Unitary Ensemble points 
(inverse temperature = power of Vandermonde interaction = 2) 
from the distribution with univariate potential 
    
        x^( a - 1 ) * ( 1 - x )^( b - 1 ) --- NOTE the -1 in the exponents

uses the Killip--Nenciu tridiagonal method described in Theorem 2.3 of the paper 
'Fast sampling from beta-ensembles' by Gautier--Bardenet--Valko arXiv:2003.02344 [stat.CO]

input: N = size of JUE matrix (it is N x N)
       num_taken = the number of smallest (sorted) eigenvalues taken out of the N total ones
       a = first Jacobi parameter (technically, see above, it's a-1)
       b = second Jacobi parameter (technically, see above, it's b-1)
       num_samples = number of independent samples

output: Lambda = num_samples x num_taken (usually N = all), each of its columns an independent JUE sample
"""
function JUE_eval_samples(N::Int64, num_taken::Int64, a::Float64, b::Float64, num_samples::Int64)::Array{Float64,2}
    C = zeros(2*N-1, num_samples) 
    A = zeros(N, num_samples)
    B = zeros(N-1, num_samples)
    unos = ones(num_samples)

    # each column of C corresponds to idependent vectors of independent entries c[n] as in the ref. paper 
    for n in 1:(N-1) 
        C[2*n-1, :] = rand(Beta((N-n) + a, (N-n) + b), num_samples)
        C[2*n, :] = rand(Beta((N-n), (N-n-1) + a + b), num_samples)
    end
    C[2*N-1, :] = rand(Beta(a, b), num_samples)

    # TODO: optimize code below
    A[1, :] = C[1, :]
    for n in 2:N 
        A[n, :] = (unos - C[2*n-3, :]) .* C[2*n-2, :] + (unos - C[2*n-2, :]) .* C[2*n-1, :]
    end

    B[1, :] = C[1, :] .* (unos - C[1, :]) .* C[2, :]
    for n in 1:(N-1)
        B[n, :] = (unos - C[2*n-2, :]) .* C[2*n-1, :] .* (unos - C[2*n-1, :]) .* C[2*n, :]
    end

    # each row of Lambda will be one vector of num_taken GUE points ordered ascendingly 
    Lambda = zeros(num_samples, num_taken) 

    for k in 1:num_samples
        J = SymTridiagonal(A[:, k], sqrt.(B[:, k])) # construct the J(acobi) tridiagonal matrix
        Lambda[k, :] = eigvals(J, range(1, length=num_taken)) # TODO: bound checks on num_taken
    end

    return Lambda
end