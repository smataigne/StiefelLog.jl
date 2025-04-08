using LinearAlgebra, SkewLinearAlgebra, MatrixEquations, Plots, LaTeXStrings, Colors
import .Manifold as m

n = 80; p = 30
δ = 0.05
U₁ = m.randorthonormal(Float64, n, p)
U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n)*δ))*U₁.Q, false)#m.randorthonormal(Float64, n, p)
display(norm(U₁-U₂)/sqrt(p)/2)
S₁ = m.StiefelVector(U₁, 0)
S₂ = m.StiefelVector(U₂, 0)


m.buildcomplement(S₁)
m.buildcomplement(S₂)
betas = Vector(0.5:0.1:1)
N = length(betas)

#Plot for different betas
colors = distinguishable_colors(8, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
styles =[:solid, :dash, :dot, :dashdot, :dashdotdot,:solid]
for k ∈ [0, 1, 2]
    P= plot(framestyle=:box, legend=:topright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
    margin = 0.3Plots.cm, yaxis=:log,gridlinewidth=1, gridalpha=0.3,minorgrid=true, minorgridalpha=0.02)
    count = 1
    for (i, β) ∈ enumerate(betas)
        norms, Niter = m.logβ(S₁, S₂, β, k)[3:4]
        plot!(1:2, (norms[1:2,2].+norms[1:2,1])./(norms[1,2]+norms[1,1]), linestyle=styles[i], color = colors[count],linewidth=1, label = L"\beta=" * string(β))
        plot!(1:Niter, (norms[1:Niter,2].+norms[1:Niter,1])./(norms[1,2]+norms[1,1]), linestyle=styles[i], color = colors[count],linewidth=2, label = false)
        count +=1
    end
    ylabel!(L"(\Vert C_k\Vert_{\mathrm{F}}+\Vert\widehat{A}_k-A_k\Vert_{\mathrm{F}})^*")
    xlabel!("Iteration")
    display(P)
    script_dir = @__DIR__
    path = joinpath(script_dir, "./figures/figure_5_generalalgorithm_subiter_"*string(k)*".pdf")
    savefig(P, path)
end