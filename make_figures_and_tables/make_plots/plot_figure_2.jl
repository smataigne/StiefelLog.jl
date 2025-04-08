using LinearAlgebra, Plots, LaTeXStrings

function τ(β)
    return abs(1 - 2β)
end

function η(δ, β)
    return τ(β) * ((1 + δ * τ(β)+ 2/3 * (δ *τ(β) )^2)+ 2 * (τ(β) * δ)^2 * log(1/(1 - 2 * δ * τ(β))))
end

function α(δ, β)
    a = δ * (1 + τ(β))
    return a^4 / (1 - a)
end

function ξ(δ, β)
    return 2β - η(δ, β) * (1 + 2β*δ + (4β*β/ 3 + 1/ 6)*δ*δ + (δ^3)/(6-δ*δ))
end

function κ(δ, β)
    return max(6 / (6 - δ * δ), η(δ, β) * δ * δ / (ξ(δ, β) * (6-δ*δ)))
end

function checkcondition(δ, β)
    c5 = (2 * δ * τ(β) < 1 ? 1 : 0)
    iszero(c5) && return 0
    c1 = (δ * (1 + τ(β)) < 1 ? 1 : 0)
    c2 = (ξ(δ, β) > 0 ? 1 : 0)
    c3 = (1 - η(δ, β) * α(δ, β) / ξ(δ, β) > 0 ? 1 : 0)
    t = η(δ, β) / ξ(δ, β)
    c4 = ( (t * δ^4)/(6*(6-δ*δ)) + (1 + t) * κ(δ, β) * α(δ, β)/(1- t * α(δ, β)) < 1 ? 1 : 0)
    return c1 * c2 * c3 * c4 * c5
end

function bisection(β)
    δₗ = 0; δₕ = 1
    δₘ = (δₗ + δₕ) / 2
    ϵ = 0.001
    while δₕ - δₗ > ϵ
        if checkcondition(δₘ, β) == 0
            δₕ = δₘ
        else
            δₗ = δₘ 
        end
        δₘ = (δₗ + δₕ) / 2
    end
    return δₘ
end

N = 1000
betas = LinRange(0, 2, N)
deltas = zeros(N)

for i ∈ 1:N
    deltas[i] = bisection(betas[i])
end
P= plot(framestyle=:box, legend=:topright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
margin = 0.3Plots.cm,gridlinewidth=1, gridalpha=0.3,minorgrid=false, minorgridalpha=0.05, titlefontfamily="Computer Modern", titlefontsize =18)

plot!(betas, deltas, linestyle=:dash, label="Admissible δ", linewidth = 2, fillrange =zeros(N), fillalpha =0.3)
xlabel!(L"\beta")
ylabel!(L"\delta")
script_dir = @__DIR__
path = joinpath(script_dir, "./figures/figure_2_condition_delta.pdf")
savefig(P, path)
display(P)





    