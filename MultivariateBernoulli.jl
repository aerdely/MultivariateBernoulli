### Multivariate Bernoulli Distribution
### Author: Dr. Arturo Erdely
### Version: 2024-09-29
### Requires: `Distributions` package

using Distributions 


## Auxiliary functions

"""
    intbin(i::Int, n::Int)

Converts the absolute value of an integer number `i` into its `n` positions binary representation.

## Example 
```
i = 6; n = 5;
v = intbin(i, n);
println(v)
# just checking:
sum(v .* (2 .^ collect((n-1):-1:0)))
```
"""
function intbin(i::Int, n::Int)
    i = abs(i)
    n = max(1, abs(n))
    b = fill(-1, n)
    if i ≤ 2^n - 1
        for k ∈ 1:n
            d, r = divrem(i, 2)
            b[n - k + 1] = r
            i = d            
        end
    else
        @error "Impossible a size $n binary representation for $i"
    end
    return b
end


"""
    binint(b::Vector{<:Integer}) 

Convert a vector `b` of zeros and ones representing a binary number to a positive integer.

## Example 
```
v = [1,1,0,1]
binint(v)
# same as:
sum(v * (2 .^ [3,2,1,0]))
```
"""
function binint(b::Vector{<:Integer})
    # convert: binary --> integer
    b = abs.(b)
    n = length(b)
    i = 0
    for k ∈ 1:n
        i += (2^(n-k))*b[k]
    end
    return i
end


"""
    powerset(x::Vector{T}) where T

Power set of a vector `x`.

## Example 
```
powerset(['a', 'b', 'c', 'd'])
```
"""
function powerset(x::Vector{T}) where T
    result = Vector{T}[[]]
    for elem in x, j in eachindex(result)
        push!(result, [result[j] ; elem])
    end
    result
end


"""
    subconj(n::Int) 

Reverse sort by length of the powerset of {1,…,n} exluding the empty set.

> Dependence: `powerset` function

## Example 
```
subconj(4)
```
"""
function subconj(n::Int)
    # revsort by length of power set of 1:n
    # depends: `powerset` function
    ℘ = powerset(collect(1:n))[2:end] # excluding the empty set
    ℓ = length.(℘)
    A = ℘[sortperm(ℓ, rev = true)]
    return A
end


## Main functions

"""
    MBerDep(pp::Vector{<:Real}; warn = true)

Calculates dependence parameters and measures for a `n`-dimensional Multivariate Bernoulli
joint distribution with vector of parameters `pp` where `length(pp) == 2^n` for some
positive integer `n ≥ 2`, `all(0.0 .≤ pp .≤ 1.0) == true`, and `sum(pp) == 1.0`.
If `sum(pp) ≠ 1.0` then normalization is applied to transform `pp` to force the sum be
equal to 1. If `warn = true` (default) then a message of normalization is displayed in case it 
is required, otherwise turn it off specifying `warn = false`. Returns a named tuple:
- `nrvs` = dimension of the random vector 
- `binprob` = a tuple where `binprob.idx` are the possible values of the random vector, `binprob.value` their probabilities, and `binprob.dic` a dictionary of such pairs of values.
- `dparam` = a tuple where `dparam.idx` are the indexes of the random variables of each parameter, `dparam.value` the value of such parameter, and `dparam.dic` a dictionary of such pairs of values.
- `dmeas` = a tuple where `dmeas.idx` are the indexes of the random variables of each dependence measure, `dmeas.value` the value of such measures, and `dmeas.dic` a dictionary of such pairs of values.

> Dependecies: `intbin` and `subconj` functions

## Example 
```
X = MBerDep([0.10,0.03,0.20,0.05,0.18,0.07,0.25,0.12]);
X.nrvs
[X.binprob.idx X.binprob.value]
[X.dparam.idx X.dparam.value]
[X.dmeas.idx X.dmeas.value]
```
"""
function MBerDep(pp::Vector{<:Real}; warn = true)
    k = length(pp)
    if !ispow2(k)
        @error "Vector length must be a power of 2"
        return nothing
    end
    if minimum(pp) < 0.0 || maximum(pp) > 1.0
        @error "Non valid probability value(s)"
        return nothing 
    end
    if sum(pp) ≠ 1.0
        p = pp ./ sum(pp)
        if warn
            @warn "Sum of vector ≠ 1.0, normalization was applied"
        end
    else
        p = pp
    end
    n = convert(Int, round(log(k)/log(2))) # 2^n == length(pp)
    b = Vector{Int64}[]
    Mb = zeros(Int, 2^n, n)
    for r ∈ 1:(2^n) # binary codification of parameters
        bin = intbin(r-1, n)
        push!(b, bin)
        Mb[r, :] = bin
    end
    Dbp = Dict(b .=> p) # dictionary
    A = subconj(n) # subsets of Bernoulli rvs
    nA = length(A)
    θ = zeros(nA)
    θ[1] = p[1]
    for r ∈ 2:nA # calculating dependence parameters θA
        idx = collect(1:nA)
        for j ∈ A[r]
            idx = idx ∩ findall(Mb[:, j] .== 0)
        end
        for i ∈ idx
            θ[r] += Dbp[Mb[i, :]]
        end
    end
    DAθ = Dict(A .=> θ) # dictionary
    idxA2 = findall(length.(A) .≥ 2) # exclude singletons
    A2 = A[idxA2]
    nA2 = length(A2)
    μ = zeros(nA2)
    for r ∈ 1:nA2 # calculating dependence measures
        θi = Float64[]
        for rr ∈ A2[r]
            push!(θi, DAθ[[rr]])
        end
        Πθ = prod(θi)
        θA = DAθ[A2[r]]
        if θA ≥ Πθ
            csup = minimum(θi)
            μ[r] = (θA - Πθ) / (csup - Πθ)
        else
            cinf = max(sum(θi) - length(θi) + 1, 0)
            μ[r] = (θA - Πθ) / (Πθ - cinf)
        end
    end
    DA2μ = Dict(A2 .=> μ) # dictionary 
    return (nrvs = n, binprob = (idx = b, value = p, dic = Dbp), 
            dparam = (idx = A, value = θ, dic = DAθ), 
            dmeas = (idx = A2, value = μ, dic = DA2μ)
    )
end


"""
    MBerMargin(pp::Vector{<:Real}, idxm::Vector{<:Integer}; warn = false)

Calculates marginal dependence parameters and measures for a `n`-dimensional Multivariate Bernoulli
joint distribution with vector of parameters `pp` where `length(pp) == 2^n` for some
positive integer `n ≥ 2`, `all(0.0 .≤ pp .≤ 1.0) == true`, `sum(pp) == 1.0`, and vector `idxm`
with subindexes of the marginal variables, where `idxm` must be a non-empty ordered subset of `1:n`.
If `sum(pp) ≠ 1.0` then normalization is applied to transform `pp` to force the sum be
equal to 1. If `warn = true` (default) then a message of normalization is displayed in case it 
is required, otherwise turn it off specifying `warn = false`. Returns a named tuple:
- `nrvs` = dimension of the marginal random vector 
- `binprob` = a tuple where `binprob.idx` are the possible values of the marginal random vector, `binprob.value` their probabilities, and `binprob.dic` a dictionary of such pairs of values.
- `dparam` = a tuple where `dparam.idx` are the indexes of the random variables of each parameter, `dparam.value` the value of such parameter, and `dparam.dic` a dictionary of such pairs of values.
- `dmeas` = a tuple where `dmeas.idx` are the indexes of the marginal random variables of each dependence measure, `dmeas.value` the value of such measures, and `dmeas.dic` a dictionary of such pairs of values.

> Dependencies: `MBerDep` and `intbin` functions

## Example:
```
X = MBerDep([0.10,0.03,0.20,0.05,0.18,0.07,0.25,0.12], warn = false);
[X.binprob.idx X.binprob.value]
X13 = MBerMargin(X.binprob.value, [1,3], warn = false);
[X13.binprob.idx X13.binprob.value]
```
"""
function MBerMargin(pp::Vector{<:Real}, idxm::Vector{<:Integer}; warn = false)
    X = MBerDep(pp, warn = warn)
    n = length(idxm)
    if !(allunique(idxm) && issorted(idxm) && idxm ⊆ 1:X.nrvs && 1 ≤ n ≤ X.nrvs)
        @warn "`idxm` must be a non-empty ordered subset of 1:$(X.nrvs)"
        return nothing
    end
    if X.nrvs == n 
        return X 
    end
    pm = zeros(2^n)
    for i ∈ 1:(2^n)
        b = intbin(i-1, n)
        id = findall(x -> x[idxm] == b, X.binprob.idx)
        pm[i] = sum(X.binprob.value[id])
    end
    return MBerDep(pm, warn = false)
end


"""
    MBerCond(pp::Vector{<:Real}, ipred::Vector{<:Integer}, icond::Vector{<:Integer}, condval::Vector{<:Integer}; warn = false) 

Calculates conditional probabilities, dependence parameters and measures of a `n`-dimensional
Multivariate Bernoulli joint distribution with vector of parameters `pp` (see `MBerDep` function)
for random subvector of variable ordered subindexes `ipred` conditional on the values `condval` given by 
the variables with ordered indexes `icond`. Returns a tuple with the characteristics explained for
the `MBerDep` function.

> Dependencies: `MBerDep`, `MBerMargin`, and `intbin` functions

## Example 
```
# joint distribution (X1,X2,X3)
X123 = MBerDep([0.10,0.03,0.20,0.05,0.18,0.07,0.25,0.12]);
X123.nrvs
[X123.binprob.idx X123.binprob.value]
[X123.dparam.idx X123.dparam.value]
[X123.dmeas.idx X123.dmeas.value]
# X2 | X1=1, X3=0
X2g13 = MBerCond(X123.binprob.value, [2], [1,3], [1,0]);
[X2g13.binprob.idx X2g13.binprob.value]
# X1,X3 | X2=1 
X13g2 = MBerCond(X123.binprob.value, [1,3], [2], [1]);
[X13g2.binprob.idx X13g2.binprob.value]
[X13g2.dparam.idx X13g2.dparam.value]
[X13g2.dmeas.idx X13g2.dmeas.value]
# X3 | X1=0 
X3g1 = MBerCond(X123.binprob.value, [3], [1], [0]);
[X3g1.binprob.idx X3g1.binprob.value]
```
"""
function MBerCond(pp::Vector{<:Real}, ipred::Vector{<:Integer}, icond::Vector{<:Integer}, condval::Vector{<:Integer}; warn = false)
    if !all(0 .≤ condval .≤ 1)
        @warn "Required: each `condval` ∈ {0,1}"
        return nothing
    end
    if length(icond) ≠ length(condval)
        @warn "vectors `ipred` and `icond` must have the same length"
    end
    if length(ipred ∩ icond) > 0 || !(issorted(ipred) && issorted(icond))
        @warn "vectors `ipred` and `icond` must be ordered disjoint sets"
        return nothing
    end
    ijoint = sort(ipred ∪ icond)
    njoint = length(ijoint)
    iipred = Int[]
    iicond = Int[]
    for k ∈ 1:njoint
        if ijoint[k] ∈ ipred 
            push!(iipred, k)
        else
            push!(iicond, k)
        end
    end
    Xjoint = MBerMargin(pp, ijoint, warn = warn)
    Xcond = MBerMargin(pp, icond, warn = warn)
    npred = length(iipred)
    ppred = zeros(2^npred)
    for i ∈ 1:(2^npred)
        bpred = intbin(i-1, npred)
        bjoint = zeros(Int, njoint)
        bjoint[iipred] = bpred 
        bjoint[iicond] = condval 
        ppred[i] = Xjoint.binprob.dic[bjoint] / Xcond.binprob.dic[condval]
    end
    return MBerDep(ppred, warn = warn)
end


"""
    MBerSim(pp::Vector{<:Real}, m::Integer; warn = true)

Simulates a random sample of size `m` from a `n`-dimensional Multivariate Bernoulli
joint distribution with vector of parameters `pp` where `length(pp) == 2^n` for some
positive integer `n ≥ 2`, `all(0.0 .≤ pp .≤ 1.0) == true`, and `sum(pp) == 1.0`.
If `sum(pp) ≠ 1.0` then normalization is applied to transform `pp` to force the sum be
equal to 1. If `warn = true` (default) then a message of normalization is displayed in case it 
is required, otherwise turn it off specifying `warn = false`. Returns a `m×n` matrix with one 
simulation per row.

> Dependence: `Categorical` function in `Distributions` package

## Example 
```
MBerSim([0.10,0.03,0.20,0.05,0.18,0.07,0.25,0.12], 15)
```
"""
function MBerSim(pp::Vector{<:Real}, m::Integer; warn = true)
    k = length(pp)
    if !ispow2(k)
        @error "Vector length must be a power of 2"
        return nothing
    end
    if minimum(pp) < 0.0 || maximum(pp) > 1.0
        @error "Non valid probability value(s)"
        return nothing 
    end
    if sum(pp) ≠ 1.0
        p = pp ./ sum(pp)
        if warn
            @warn "Sum of vector ≠ 1.0, normalization was applied"
        end
    else
        p = pp
    end
    n = convert(Int, round(log(k)/log(2))) # 2^n == length(pp)
    sim = rand(Categorical(p), m) .- 1
    matsim = zeros(Bool, m, n)
    for i ∈ 1:m
        matsim[i, :] = intbin.(sim[i], n)
    end
    return matsim
end


"""
    MBerBayes(sample::Matrix{<:Integer}, prior = 1/2) 

Posterior Dirichlet model for the parameters of a `n`-dimensional Multivariate Bernoulli
distribution given an observed random sample of size `m`. The `sample` must be a `m×n` matrix of zeros 
and ones, where each row is an observation of the `n`-dimensional random vector. Returns 
a `Dirichlet` model as defined in the `Distributions` package. The `prior` constant multiplies
a vector of ones to be used as a prior distribution.

> Dependencies: `Dirichlet` function in `Distributions` package, `binint` function

## Example 
```
sample = [0 0 1; 0 1 0; 1 1 0; 0 0 0; 0 1 0; 0 0 1]
post = MBerBayes(sample)
```
"""
function MBerBayes(sample::Matrix{<:Integer}, prior = 1/2)
    if !all(x -> x ∈ [0,1], sample)
        @warn "Not a valid sample from a Multivariate Bernoulli"
        return nothing
    end
    m, n = size(sample)
    α = prior * ones(Int, 2^n)
    nr = zeros(Int, 2^n)
    for i ∈ 1:m
        nr[1 + binint(sample[i, :])] += 1
    end
    return Dirichlet(α + nr)
end


"""
    MBerEst(sample::Matrix{<:Integer}; prior = 1/2)

Posterior point estimation for probabilities `probs`, dependence parameters `dparam`, and dependence
measures `dmeas` given a `sample` from a Multivariate Bernoulli distribution, and a `prior`
value that multiplies a vector of ones to be used as a prior distribution. Returns a named tuple:
- `nrvs` = dimension of the random vector 
- `binprob` = a tuple where `binprob.idx` are the possible values of the random vector, `binprob.value` their mean posterior estimations, and `binprob.dic` a dictionary of such pairs of values.
- `dparam` = a tuple where `dparam.idx` are the indexes of the random variables of each parameter, `dparam.value` the value of such parameter posterior estimation, and `dparam.dic` a dictionary of such pairs of values.
- `dmeas` = a tuple where `dmeas.idx` are the indexes of the random variables of each dependence measure, `dmeas.value` the value of such measures posterior point estimations, and `dmeas.dic` a dictionary of such pairs of values.

> Dependencies: `MBerDep` and `MBerBayes` functions

## Example
```
sample = [0 0 1; 0 1 0; 1 1 0; 0 0 0; 0 1 0; 0 0 1]
estim = MBerEst(sample);
estim.nrvs
[estim.binprob.idx estim.binprob.value]
[estim.dparam.idx estim.dparam.value]
[estim.dmeas.idx estim.dmeas.value]
```
"""
function MBerEst(sample::Matrix{<:Integer}; prior = 1/2)
    if !all(x -> x ∈ [0,1], sample)
        @warn "Not a valid sample from a Multivariate Bernoulli"
        return nothing
    end
    post = MBerBayes(sample, prior)
    pp = post.alpha ./ post.alpha0
    return MBerDep(pp, warn = false)
end


"""
    MBerInf(sample::Matrix{<:Integer}; prior = 1/2, nsim = 10_000, probint = 0.95) 

Posterior point and `probint` interval estimation for probabilities `probs`, 
dependence parameters `dparam`, and dependence measures `dmeas` given a `sample` from
a Multivariate Bernoulli distribution, and a `prior` value that multiplies a vector
of ones to be used as a prior distribution. Returns a tuple with results:
- `probint` = probability of intervals containing the parameter value.
- `prior` = constant multiplying the prior vector of ones.
- `probs` = 4-column matrix where column 1 has the possible values of the random vector, columns 2 and 4 are the extremes of a `probint` interval using the `(1-probint)/2` and `(1+probint)/2` quantiles, and column 3 is the median posterior point estimation, for the probability parameters.
- `dparam` = 4-column matrix where column 1 has the variables' indexes, columns 2 and 4 are the extremes of a `probint` interval using the `(1-probint)/2` and `(1+probint)/2` quantiles, and column 3 is the median posterior point estimation, for the dependence parameters.
- `dmeas` = 4-column matrix where column 1 has the variables' indexes, columns 2 and 4 are the extremes of a `probint` interval using the `(1-probint)/2` and `(1+probint)/2` quantiles, and column 3 is the median posterior point estimation, for the dependence measures.

> Dependencies: functions `MBerDep` `MBerBayes`

> **Warning**: For dimension 10 or higher calculations may take more than an hour. In such case you may consider using `MBerEst` just for point estimations.

## Example
```
sample = [0 0 1; 0 1 0; 1 1 0; 0 0 0; 0 1 0; 0 0 1]
inference = MBerInf(sample);
inference.probint 
inference.probs
inference.dparam 
inference.dmeas 
```
"""
function MBerInf(sample::Matrix{<:Integer}; prior = 1/2, nsim = 10_000, probint = 0.95)
    if !all(x -> x ∈ [0,1], sample)
        @warn "Not a valid sample from a Multivariate Bernoulli"
        return nothing
    end
    if !(0 < probint < 1)
        @warn "Not a valid probability for interval estimations"
        return nothing
    end
    γ = 1 - probint 
    post = MBerBayes(sample, prior)
    simPost = rand(post, nsim)
    # probabilities
    postInf = zeros(size(simPost)[1], 3);
    for i ∈ 1:size(simPost)[1] # 100(1-γ)% interval and point estimation
        postInf[i, :] = quantile(simPost[i, :], [γ/2, 0.5, 1-γ/2])
    end
    X = MBerDep(postInf[:, 2], warn = false)
    p = [X.binprob.idx postInf]
    # dependence parameters and measures
    dparam = zeros(length(X.dparam.idx), size(simPost)[2])
    dmeas = zeros(length(X.dmeas.idx), size(simPost)[2])
    for j ∈ 1:size(simPost)[2]
        Xj = MBerDep(simPost[:, j], warn = false)
        dparam[:, j] = Xj.dparam.value
        dmeas[:, j] = Xj.dmeas.value
    end
    dparamInf = zeros(length(X.dparam.idx), 3)
    dmeasInf = zeros(length(X.dmeas.idx), 3)
    for i ∈ 1:size(dparamInf)[1] # 100(1-γ)% interval and point estimation
        dparamInf[i, :] = quantile(dparam[i, :], [γ/2, 0.5, 1-γ/2])
    end
    for i ∈ 1:size(dmeasInf)[1] # 100(1-γ)% interval and point estimation
        dmeasInf[i, :] = quantile(dmeas[i, :], [γ/2, 0.5, 1-γ/2])
    end
    d = [X.dparam.idx dparamInf]
    μ = [X.dmeas.idx dmeasInf]
    # results
    return (probs = p, dparam = d, dmeas = μ, probint = probint, prior = prior)
end



## Menu of functions

"""
    MBerMenu() 

Displays de names of the main functions in `MultBernoulli.jl`
"""
function MBerMenu()
    @info "MBerDep  MBerMargin  MBerCond  MBerSim  MBerBayes  MBerEst  MBerInf  MBerMenu  binint  intbin  powerset  subconj"
    println("--> Depends on `Distributions` package")
end


println("Execute `MBerMenu()` to display the list of loaded functions:")
MBerMenu()
