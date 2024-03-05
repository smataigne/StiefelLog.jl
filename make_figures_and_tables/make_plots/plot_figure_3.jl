using LinearAlgebra, SkewLinearAlgebra, MatrixEquations, Plots, LaTeXStrings, Colors
import .Manifold as m

n = 120; p = 50
δ = 0.05
U₁ = m.randorthonormal(Float64, n, p)
U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n)*δ))*U₁.Q, false)#m.randorthonormal(Float64, n, p)
display(norm(U₁-U₂)/sqrt(p)/2)
S₁ = m.StiefelVector(U₁, 0)
S₂ = m.StiefelVector(U₂, 0)


betas = Vector(0.5:0.1:1)
colors = distinguishable_colors(8, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
styles =[:solid, :dash, :dot, :dashdot, :dashdotdot,:solid]

P= plot(framestyle=:box, legend=:topright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
margin = 0.3Plots.cm, yaxis=:log,gridlinewidth=1, gridalpha=0.3,minorgrid=true, minorgridalpha=0.05)
count = 1
for i∈eachindex(betas)
    global count
    β = betas[i]
    if β != 0.5
        norms, Niter = m.logβmomentum(S₁, S₂, β)[3:4]
        plot!(1:Niter, (norms[1:Niter,2].+norms[1:Niter,1])./(norms[1,2]+norms[1,1]), color = colors[count],linewidth=2, label = false)
    end
    norms, Niter = m.logβ(S₁, S₂, β, 0)[3:4]
    plot!(1:Niter, (norms[1:Niter,2].+norms[1:Niter,1])./(norms[1,2]+norms[1,1]), linestyle=:dash, color = colors[count],linewidth=2, label = L"\beta=" * string(β))
    count +=1
end
ylabel!(L"(\Vert C_k\Vert_{\mathrm{F}}+\Vert\widehat{A}_k-A_k\Vert_{\mathrm{F}})^*")
xlabel!("Iteration")
display(P)
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_3_acceleration_"*string(floor(Int,100*norm(U₁-U₂)/sqrt(p)/2))*".pdf")
savefig(P, path)
