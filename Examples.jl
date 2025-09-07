### Examples for Multivariate Bernoulli Distribution
### Author: Dr. Arturo Erdely
### Version: 2025-09-06

begin
    using Distributions, Plots, LaTeXStrings, Random, CSV, DataFrames
    include("MultivariateBernoulli.jl")
end;


## Bivariate Bernoulli figures

begin # compatible θ12 values given (θ1,θ2)
    plot([0,1],[0,1], color = :gray, lw = 0.5, legend = false, size = (417,400), 
         grid = false, xticks = [0,1], yticks = [0,1], xlabel = L"\theta_1", ylabel = L"\theta_2"
    )
    plot!([0,1], [1,0], color = :gray, lw = 0.5)
    plot!([0,1], [0,0], color = :gray, lw = 0.5)
    plot!([1,1], [0,1], color = :gray, lw = 0.5)
    plot!([0,1], [1,1], color = :gray, lw = 0.5)
    plot!([0,0], [0,1], color = :gray, lw = 0.5)
    θ1,θ2 = 0.78,0.58
    scatter!([θ1],[θ2], ms = 6, mc = :black)
    plot!([θ1+θ2-1,1],[1, θ1+θ2-1], color = :black, lw = 1.5)
    m = min(θ1,θ2)
    plot!([m, 1], [m, m], color = :black, lw = 1.5)
    plot!([m, m], [m, 1], color = :black, lw = 1.5)
    plot!([m,m], [0,m], color = :violet, lw = 1.5)
    w = θ1+θ2-1
    plot!([w,w], [0,1], color = :violet, lw = 1.5)
    plot!([w,m], [0,0], color = :red, lw = 7)
    annotate!([(θ1+0.1,θ2+0.02,(L"(\theta_1^{*},\theta_2^{*})", :bottom, 11))])
    P1 = annotate!([((w+m)/2, 0.02,(L"\theta_{12}", :bottom, :red, 12))])
    println()
end

begin # compatible values (θ1,θ2) given θ12
    plot([0,1],[0,1], color = :white, lw = 0.5, legend = false, size = (417,400), 
         grid = false, xticks = [0,1], yticks = [0,1], xlabel = L"\theta_1", ylabel = L"\theta_2"
    )
    # plot!([0,1], [1,0], color = :gray, lw = 0.5)
    plot!([0,1], [0,0], color = :gray, lw = 0.5)
    plot!([1,1], [0,1], color = :gray, lw = 0.5)
    plot!([0,1], [1,1], color = :gray, lw = 0.5)
    plot!([0,0], [0,1], color = :gray, lw = 0.5)
    θ12 = 0.3
    plot!([θ12,θ12],[0,θ12], color = :black, lw = 1.5)
    plot!([0,θ12],[θ12,θ12], color = :black, lw = 1.5)
    scatter!([θ12],[0], mc = :black, ms = 5)
    annotate!([(θ12+0.08, 0.01,(L"\theta_{12}^{*}", :bottom, 12))])
    plot!([θ12, 1], [θ12, θ12], fillrange = [1, θ12], fillcolor = :red, fillalpha = 1)
    plot!([θ12,1], [θ12, θ12], color = :red, lw = 1.5)
    plot!([θ12, θ12], [θ12,1], color = :red, lw = 1.5)
    P2 = plot!([θ12, 1], [1, θ12], color = :red, lw = 1.5)
    println()
end

begin # Figure 1
    println("Figure 1 plots")
    plot(P1, P2, layout = (1,2), legend = false, size = (600,300))
    savefig("Figure1.png")
    savefig("Figure1.pdf")
    println()
end


## Example 1

begin
    g0(θ) = θ*(2θ-1)
    g1(θ) = θ^2
    g2(θ) = 1 - g0(θ)
    g3(θ) = 1 - 3*θ + 3*(θ^2)
    θ = range(0.0, 1.0, length = 1_001)
    plot(θ, θ.^3, label = "indep", lw = 1)
    plot!(θ, g0.(θ), label = "g0", lw = 2)
    plot!(θ, g1.(θ), label = "g1", lw = 2)
    plot!(θ, g2.(θ), label = "g2", lw = 2)
    plot!(θ, g3.(θ), label = "g3", lw = 2)
    # vline!([(1+√13)/6], label = "", color = :gray)
    vline!([4/5], label = "", color = :gray)
    vline!([(1+√5)/4], label = "", color = :gray)
    plot!(θ, θ, label = "", color = :gray)
    plot!(θ, 0.0 .*θ, label = "", color = :gray)
    θ1 = range(0.0, 2/3, length = 100)
    plot!(θ1, 0.0 .*θ1, label = "", color = :gray)
    θ2 = collect(range(2/3,1.0,length = 1000))
    plot!(θ2, 3.0.*θ2.-2, label = "", color = :gray)
    println()
end

begin
    println("Figure 2 plot")
    θ = range(0, 1, length = 1_000)
    θ1 = range(0, 1/2, length = 1_000)
    θ2 = range(1/2, 4/5, length = 1_000)
    θ3 = range(4/5, (1+√5)/4, length = 1_000)
    θw = range(2/3, 1, length = 1_000)
    plot(θ, θ, color = :gray, lw = 2, label = "Fréchet-Hoeffding bounds", xticks = [0,1], yticks = [0,1], 
         xlabel = L"\theta", ylabel = L"\theta_{123}", size = (400,400)
    )
    plot!([0,2/3], [0,0], color = :gray, lw = 2, label = "")
    plot!(θw, 3 .* θw .- 2, color = :gray, lw = 2, label = "")
    plot!(θ, θ .^ 3, color = :blue, lw = 2, label = "independence")
    θr = vcat(θ1,θ2,θ3)
    θ123sup = zeros(length(θr))
    θ123inf = zeros(length(θr))
    θ123inf[1001:end] = g0.(θr[1001:end])
    θ123sup[1:1000] = g1.(θr[1:1000])
    θ123sup[1001:2000] = g3.(θr[1001:2000])
    θ123sup[2001:3000] = g2.(θr[2001:3000])
    plot!(θr, θ123inf, fillrange = θ123sup, fillcolor = :red, label = "compatibility region")
    plot!(θr, θ123inf, color = :red, label = "")
    plot!(θ, θ .^ 3, color = :blue, lw = 2, label = "")
    savefig("Figure2.png")
    savefig("Figure2.pdf")
    println()
end



### Multivariate Bernoulli examples

## Section 6.1: A theoretical example

begin
    println("Example Section 6.1")
    pp = [0.15,0.21,0.21,0.03,0.21,0.03,0.03,0.13];
    X = MBerDep(pp);
    println("Parameters:") # Table 2
    display([X.dparam.idx X.dparam.value])
    println("Dependencies:") # Table 2
    display([X.dmeas.idx X.dmeas.value])
    println()
end



## Section 6.2: A simulation example

begin                                                           
    println("Example 6.2")
    p000 = 0.1                                                  
    p001, p010, p100 = 0.2, 0.1, 0.05                           
    p011, p101, p110 = 0.2, 0.15, 0.1                           
    p111 = 1 - sum([p000, p001, p010, p011, p100, p101, p110])  
    pp = [p000, p001, p010, p011, p100, p101, p110, p111]
    X = MBerDep(pp)
end;

begin
    Random.seed!(1234) # for reproducibility
    simX = MBerSim(pp, 3_000);
    infer = MBerInf(simX, prior = 1/2, nsim = 100_000, probint = 0.99);
    println("Estimated probabilities:") # Table 3
    display(infer.probs) # estimated probabilities
    println("Theoretical probabilities:") # Table 3
    display([X.binprob.idx X.binprob.value]) # theoretical probabilities
    println("Estimated parameters:") # Table 4
    display(infer.dparam) # estimated dependence parameters
    println("Theoretical parameters:") # Table 4
    display([X.dparam.idx X.dparam.value]) # theoretical dependence parameters
    println("Estimated dependencies:") # Table 5
    display(infer.dmeas) # estimated dependence measures
    println("Theoretical dependencies:") # Table 5
    display([X.dmeas.idx X.dmeas.value]) # theoretical dependence measures
    println()
end

# coverage rates
@time begin # WARNING: 10_000 simulations may take around 3.5 hours
    nsim = 10_000
    nobs = 3_000
    cover_probs = zeros(Bool, nsim, 8)
    cover_param = zeros(Bool, nsim, 7)
    cover_meas = zeros(Bool, nsim, 4)
    for i ∈ 1:nsim
        simX = MBerSim(pp, nobs)
        infer = MBerInf(simX, prior = 1/2, nsim = 100_000, probint = 0.99)
        for j ∈ 1:8
            cover_probs[i,j] = (infer.probs[j,2] ≤ pp[j] ≤ infer.probs[j,4])
        end
        for j ∈ 1:7
            cover_param[i,j] = (infer.dparam[j,2] ≤ X.dparam.value[j] ≤ infer.dparam[j,4])
        end
        for j ∈ 1:4
            cover_meas[i,j] = (infer.dmeas[j,2] ≤ X.dmeas.value[j] ≤ infer.dmeas[j,4])
        end
    end
    println("Coverage rates (99% intervals):")
    println("Probabilities: ", mean(cover_probs, dims = 1))
    println("Parameters:    ", mean(cover_param, dims = 1))
    println("Dependencies:  ", mean(cover_meas, dims = 1))
    println()
end



## Section 6.3: Bank churn data 

# Read data as a dataframe
begin
    println("Example Section 6.4: Bank churn data")
    df = CSV.read("churnbindata.csv", DataFrame)
    show(describe(df), allrows = true)
end

# Convert dataframe to a matrix
begin 
    binmat = zeros(Int, size(df))
    for c ∈ 1:ncol(df)
        binmat[:, c] = df[:, c]
    end
    binmat
end

# Point and 95% probability interval estimations 
@time inference4v = MBerInf(binmat, nsim = 1_000_000); # less than 1 minute
inference4v.dmeas

# Discard variable 'crcard' (column 2)
binmat = binmat[:, [1,3,4]]
@time inference3v = MBerInf(binmat, nsim = 1_000_000); # less than 1 minute

# Point and 95% probability interval estimations for...
inference3v.probs # ... probabilities
inference3v.dparam # ... parameters
inference3v.dmeas # ... dependencies

# Probabilities for exiting the bank
begin
    condProb = DataFrame(gender = [missing,1,0,missing,missing,1,1,0,0],
                         active = [missing,missing,missing,1,0,1,0,1,0],
                         Pexit = zeros(9)
    )
    condProb.Pexit[1] = 1 - inference3v.dparam[7,3]
    X123 = MBerDep(Float64.(inference3v.probs[:, 3])) # warning due to rounding
    condProb.Pexit[2] = MBerCond(X123.binprob.value, [3], [1], [1]).binprob.value[2]
    condProb.Pexit[3] = MBerCond(X123.binprob.value, [3], [1], [0]).binprob.value[2]
    condProb.Pexit[4] = MBerCond(X123.binprob.value, [3], [2], [1]).binprob.value[2]
    condProb.Pexit[5] = MBerCond(X123.binprob.value, [3], [2], [0]).binprob.value[2]  
    condProb.Pexit[6] = MBerCond(X123.binprob.value, [3], [1,2], [1,1]).binprob.value[2]
    condProb.Pexit[7] = MBerCond(X123.binprob.value, [3], [1,2], [1,0]).binprob.value[2]
    condProb.Pexit[8] = MBerCond(X123.binprob.value, [3], [1,2], [0,1]).binprob.value[2]
    condProb.Pexit[9] = MBerCond(X123.binprob.value, [3], [1,2], [0,0]).binprob.value[2]
    display(condProb)
end

# Prediction rules
function predNonCond(nsim, condProb)
    B = Bernoulli(condProb.Pexit[1])
    accuracy = 0.0
    for k ∈ 1:nsim
        accuracy += mean(rand(B, 10_000) .== binmat[:, 3]) / nsim
    end
    return accuracy
end
function predGivenGender(nsim, condProb)
    ival1 = findall(binmat[:, 1] .== 1)
    ival0 = findall(binmat[:, 1] .== 0)
    n1 = length(ival1)
    n0 = length(ival0)
    vsim = fill(-99999999, 10_000)
    accuracy = 0.0
    for k ∈ 1:nsim
        vsim[ival1] = rand(Bernoulli(condProb.Pexit[2]), n1)
        vsim[ival0] = rand(Bernoulli(condProb.Pexit[3]), n0)
        accuracy += mean(vsim .== binmat[:, 3]) / nsim
    end
    return accuracy
end
function predGivenActive(nsim, condProb)
    ival1 = findall(binmat[:, 2] .== 1)
    ival0 = findall(binmat[:, 2] .== 0)
    n1 = length(ival1)
    n0 = length(ival0)
    vsim = fill(-99999999, 10_000)
    accuracy = 0.0
    for k ∈ 1:nsim
        vsim[ival1] = rand(Bernoulli(condProb.Pexit[4]), n1)
        vsim[ival0] = rand(Bernoulli(condProb.Pexit[5]), n0)
        accuracy += mean(vsim .== binmat[:, 3]) / nsim
    end
    return accuracy
end
function predGivenGenderActive(nsim, condProb)
    X1X2 = Vector{Int}[]
    for j ∈ 1:10_000
        push!(X1X2, binmat[j, [1,2]])
    end
    ival11 = findall(X1X2 .== [[1,1]])
    ival10 = findall(X1X2 .== [[1,0]])
    ival01 = findall(X1X2 .== [[0,1]])
    ival00 = findall(X1X2 .== [[0,0]])
    n11 = length(ival11)
    n10 = length(ival10)
    n01 = length(ival01)
    n00 = length(ival00)
    vsim = fill(-99999999, 10_000)
    accuracy = 0.0
    for k ∈ 1:nsim
        vsim[ival11] = rand(Bernoulli(condProb.Pexit[6]), n11)
        vsim[ival10] = rand(Bernoulli(condProb.Pexit[7]), n10)
        vsim[ival01] = rand(Bernoulli(condProb.Pexit[8]), n01)
        vsim[ival00] = rand(Bernoulli(condProb.Pexit[9]), n00)
        accuracy += mean(vsim .== binmat[:, 3]) / nsim
    end
    return accuracy
end

# Accuracy for each rule 
nsim = 1_000_000 # 4 minutes approx 
@time predict = DataFrame(rule = ["Non-conditional","Given gender", "Given active",
                                  "Given gender & active"],
                      accuracy = [predNonCond(nsim, condProb), 
                                  predGivenGender(nsim, condProb),
                                  predGivenActive(nsim, condProb), 
                                  predGivenGenderActive(nsim, condProb)]
)

# Accuracy given gender and active status
# but setting 3-variate dependence to zero 
begin
    θ1 = inference3v.dparam[5,3]
    θ2 = inference3v.dparam[6,3]
    θ3 = inference3v.dparam[7,3]
    θ123 = θ1 * θ2 * θ3 
    μ12 = inference3v.dmeas[2,3]
    μ13 = inference3v.dmeas[3,3]
    μ23 = inference3v.dmeas[4,3]
    θ12 = (min(θ1,θ2) - θ1*θ2)*μ12 + θ1*θ2
    θ13 = (θ1*θ3 - max(θ1+θ3-1,0))*μ13 + θ1*θ3
    θ23 = (θ2*θ3 - max(θ2+θ3-1,0))*μ23 + θ2*θ3
    p000 = θ123 
    p001 = θ12 - θ123
    p010 = θ13 - θ123
    p100 = θ23 - θ123
    p011 = θ1 - θ12 - θ13 + θ123
    p101 = θ2 - θ12 - θ23 + θ123
    p110 = θ3 - θ13 - θ23 + θ123
    p111 = 1 - (p000 + p001 + p010 + p011 + p100 + p101 + p110)
    pp = [p000, p001, p010, p011, p100, p101, p110, p111]
    Xzero3dep = MBerDep(pp)
end;
# dependence measures without and with zero 3-variate dependence
[inference3v.dmeas[:, [1,3]] Xzero3dep.dmeas.value]

# Probabilities for exiting the bank with zero 3-variate dependence
begin
    condProb3 = DataFrame(gender = [missing,1,0,missing,missing,1,1,0,0],
                          active = [missing,missing,missing,1,0,1,0,1,0],
                          Pexit = zeros(9)
    )
    condProb3.Pexit[1] = 1 - Xzero3dep.dparam.value[7]
    condProb3.Pexit[2] = MBerCond(Xzero3dep.binprob.value, [3], [1], [1]).binprob.value[2]
    condProb3.Pexit[3] = MBerCond(Xzero3dep.binprob.value, [3], [1], [0]).binprob.value[2]
    condProb3.Pexit[4] = MBerCond(Xzero3dep.binprob.value, [3], [2], [1]).binprob.value[2]
    condProb3.Pexit[5] = MBerCond(Xzero3dep.binprob.value, [3], [2], [0]).binprob.value[2]  
    condProb3.Pexit[6] = MBerCond(Xzero3dep.binprob.value, [3], [1,2], [1,1]).binprob.value[2]
    condProb3.Pexit[7] = MBerCond(Xzero3dep.binprob.value, [3], [1,2], [1,0]).binprob.value[2]
    condProb3.Pexit[8] = MBerCond(Xzero3dep.binprob.value, [3], [1,2], [0,1]).binprob.value[2]
    condProb3.Pexit[9] = MBerCond(Xzero3dep.binprob.value, [3], [1,2], [0,0]).binprob.value[2]
    display(condProb3)
end

# Accuracy for each rule setting 3-variate dependence to zero
nsim = 1_000_000 # 4 minutes approx 
@time predict3 = DataFrame(rule = ["Non-conditional","Given gender", "Given active",
                                  "Given gender & active"],
                      accuracy = [predNonCond(nsim, condProb3), 
                                  predGivenGender(nsim, condProb3),
                                  predGivenActive(nsim, condProb3), 
                                  predGivenGenderActive(nsim, condProb3)]
);

# Comparison of accuracies with and without 3-variate dependence
println("Accuracy with 3-variate dependence:")
predict
println("Accuracy with zero 3-variate dependence:")
predict3



## Section 6.4: COVID-19 data in Mexico

begin
    println("Example Section 6.3: COVID-19 data")
    df = CSV.read("covid2020.csv", DataFrame)
    show(describe(df), allrows = true)
    data = zeros(Int, size(df))
    for c ∈ 1:ncol(df)
        data[:, c] = df[:, c]
    end
    println()
    display(data)
    println()
end

begin
    println("Processing data... (5 minutes approx)")
    @time estim = MBerEst(data); # 5 minutes approx
    iord = sortperm(estim.dmeas.value, rev = true)
    μidx = estim.dmeas.idx[iord]
    μ = estim.dmeas.value[iord]
    iMue = findall(x -> 15 ∈ x, μidx)
    println("Highest multivariate dependencies:") # Table 6
    display([μidx[iMue] μ[iMue]][1:20, :])
    println()
end

begin
    pd = mean(data[:, 15]) # P(death)
    # P(death|hospital, COPD, immunosup, cardio, CKD, age 65+)
    estim2 = MBerEst(data[:, [2,4,5,7,9,14,15]])
    pdcond = MBerCond(estim2.binprob.value, [7], [1,2,3,4,5,6], [1,1,1,1,1,1]).binprob.dic[[1]]
    println("P(death) = ", pd)
    println("P(death|hospital, COPD, immunosup, cardio, CKD, age 65+) = ", pdcond)
end

println()
println("Execute MBerMenu() to display the list of loaded functions from file MultivariateBernoulli.jl")
