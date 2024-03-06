using LinearAlgebra, SkewLinearAlgebra, Plots, BenchmarkTools, LaTeXStrings
import .Manifold as m 

n = 40; p = 20

betas = 0.6:0.1:1
iters = Array(0:5)
times = zeros(length(iters), length(betas))
δ = 0.15
N = 10
for k ∈ 1:N
    U₁ = m.randorthonormal(Float64, n, p)
    U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n) * δ)) * U₁.Q, false)
    S₁ = m.StiefelVector(U₁, 0)
    S₂ = m.StiefelVector(U₂, 0)
    display(norm(U₁.Q-U₂.Q)/(2sqrt(p)))
    for j ∈ eachindex(betas)
        β = betas[j]
        for i ∈ eachindex(iters)
            iter = iters[i]
            temp = minimum(@elapsed m.logβ(S₁, S₂, β, iter) for _ ∈ 1:10)
            times[i,j] = (times[i, j] * (k-1) + temp)/k
        end
    end
end

P = plot(framestyle=:box, legend=:topleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
margin = 0.3Plots.cm,gridlinewidth=1, gridalpha=0.3,minorgrid=true, minorgridalpha=0.05, titlefontfamily="Computer Modern", titlefontsize =18,ylims=(0.3,1.7))
styles= [:solid, :dash, :dot, :dashdot, :dashdotdot]
for j ∈ eachindex(betas)
    β = betas[j]
    plot!(iters, times[:, j]./times[1, j], label = L"\beta = " * string(β), linewidth=2, linestyle=styles[j])
end
xlabel!("Number of sub-iterations")
ylabel!("Running time % w.r.t PF")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_4_subiterationtest.pdf")
savefig(P, path)
display(P)


