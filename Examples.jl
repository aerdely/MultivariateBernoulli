### Examples for Multivariate Bernoulli Distribution
### Author: Dr. Arturo Erdely
### Version: 2024-09-29

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

## Example 1 revisited

begin
    println("Example 1 revisited")
    pp = [0.15,0.21,0.21,0.03,0.21,0.03,0.03,0.13];
    X = MBerDep(pp);
    println("Parameters:")
    display([X.dparam.idx X.dparam.value])
    println("Dependencies:")
    display([X.dmeas.idx X.dmeas.value])
    println()
end


## Example 3 

begin                                                           
    println("Example 3")
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
    println("Estimated probabilities:")
    display(infer.probs) # estimated probabilities
    println("Theoretical probabilities:")
    display([X.binprob.idx X.binprob.value]) # theoretical probabilities
    println("Estimated parameters:")
    display(infer.dparam) # estimated dependence parameters
    println("Theoretical parameters:")
    display([X.dparam.idx X.dparam.value]) # theoretical dependence parameters
    println("Estimated dependencies:")
    display(infer.dmeas) # estimated dependence measures
    println("Theoretical dependencies:")
    display([X.dmeas.idx X.dmeas.value]) # theoretical dependence measures
    println()
end


## Example 4: COVID-19 data

begin
    println("Example 4: COVID-19 data")
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
    println("Highest multivariate dependencies:")
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
