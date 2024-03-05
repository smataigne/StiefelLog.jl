using LinearAlgebra, SkewLinearAlgebra, MatrixEquations, Plots, LaTeXStrings, Colors
import .Manifold as m

n = 60; p = 30


#Plot for different betas
colors = distinguishable_colors(6, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
deltas = 0.01:0.01:0.15
betas = 0.6:0.1:1
N = length(deltas)
Nb = length(betas)
dist = zeros(N)
ratios = zeros(N, Nb)
styles =[:solid, :dash, :dot, :dashdot, :dashdotdot,:solid]
for t∈1:50
    for i∈1:N
        δ  = deltas[i]
        U₁ = m.randorthonormal(Float64, n, p)
        U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n)*δ))*U₁.Q, false)#m.randorthonormal(Float64, n, p)
        dist[i] = norm(U₁-U₂)/sqrt(p)/2
        S₁ = m.StiefelVector(U₁, 0)
        S₂ = m.StiefelVector(U₂, 0)
        for j∈1:Nb
            β = betas[j]
            norms, Niter = m.logβmomentum(S₁, S₂, β)[3:4]
            norms, Niter2 = m.logβ(S₁, S₂, β, 0)[3:4]
            ratios[i,j] = (ratios[i,j] *(t-1) + Niter2 / Niter)/t
        end
    end
end

P= plot(framestyle=:box, legend=(0.6, 0.93),font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
    margin = 0.3Plots.cm,gridlinewidth=1, gridalpha=0.3,minorgrid=true, minorgridalpha=0.05)

for  i∈1:Nb
    plot!(dist, ratios[:,i], label = L"\beta="*string(betas[i]), color=colors[i], linewidth=2, linestyle= styles[i])
end 

ylabel!("Improvement factor")
xlabel!(L"\Vert U-\widetilde{U}\ \Vert_\mathrm{F}/2\sqrt{p}")
display(P)
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_3_ratio_60_30.pdf")
savefig(P, path)