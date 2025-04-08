using LinearAlgebra, Plots, Distributions, LaTeXStrings, XLSX
include("Accelerated_gradient_for_regression.jl")

ps = [2, 4, 8, 16, 32]                  #values of p that are considered
radius = zeros(length(ps), 6)
colors = distinguishable_colors(6, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
styles =[:solid, :dash, :dot, :dashdot, :dashdotdot,:solid]

script_dir = @__DIR__
data_dir = "../logistic_data/"
figures_dir = "../logistic_figures/"

for j∈eachindex(ps)
    p = ps[j]
    file_name = "Logistic_data_st"*string(p)*"_1000.xlsx"
    file_path = joinpath(script_dir,data_dir, file_name)
    dtable = XLSX.readtable(file_path, "new_sheet")

    n, p, N, s = dtable.data[1][1], dtable.data[1][2],dtable.data[1][3], dtable.data[1][4]
    betas = dtable.data[2][1:s]
    x = dtable.data[3]
    s = length(betas)
    X = zeros(N, 2)
    X[:, 1] .= 1
    X[:, 2]  = x
    Y = zeros(N, s)
    for i ∈ 1:s
        Y[:,i] = dtable.data[3+i]
    end

    θ₀ = zeros(2)
    L = (norm(X)^2) / 4

    P2= plot(framestyle=:box, legend=:bottomleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
    legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
    margin = 0.3Plots.cm,gridlinewidth=1, gridalpha=0.3,minorgrid=true, minorgridalpha=0.05, titlefontfamily="Computer Modern", titlefontsize =18, ylims=(0,1))

    t = LinRange(0,1, 100)
    max_iter = 10000
    for i ∈ 1:s
        global radius
        β = betas[i]
        if β > 0.5
            θ, norms, iter = AccGradientMethod(θ₀, X, Y[:, i], L, max_iter)
            s = 1 ./(1 .+ exp.(-θ[1] .- θ[2] .* t))
            R2 = 1 - sum((Y[:, i]  - 1 ./(1 .+ exp.(-θ[1] .- θ[2] .* X[:, 2]))).^2)/sum((Y[:, i] .- mean(Y[:, i])).^2)
            display("R2 for (β, p)=(" * string(β) * ","*string(p) * ") : " * string(round(R2, digits =2)))
            plot!(t, s, label = false, linewidth = 2, linestyle = styles[i], color = colors[i])
            plot!(t[1:2], s[1:2], label = L"\beta="*string(β), linewidth = 1, linestyle = styles[i], color = colors[i])
            h = 1
            while h<101 && s[h]>0.9
                h+=1
            end
            h-=1
            h = max(h, 1)
            radius[j,i] = t[h]
        end
    end
    figure_path = joinpath(script_dir,figures_dir, "success_rate_st" * string(p) * "_1000_new.pdf")
    xlabel!(L"\Vert U-\widetilde{U} \ \Vert_\mathrm{F} / 2\sqrt{p}")
    ylabel!("Probability of success")
    savefig(P2, figure_path)
    display(P2)
end


betas = Vector(0.5:0.1:1)
P = plot(framestyle=:box, legend=:bottomleft,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern",
legendfontsize=12,guidefontfamily = "Computer Modern",yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13,
margin = 0.3Plots.cm,gridlinewidth=1, gridalpha=0.3,minorgrid=true, minorgridalpha=0.05, titlefontfamily="Computer Modern", titlefontsize =18, ylims=(0,1))
for i ∈ eachindex(betas)
    β = betas[i]
    if β > 0.5
        global radius
        plot!(log2.(ps), radius[:,i], color=colors[i], linestyle = styles[i], label = false,linewidth=2)
        plot!(log2.(ps[1:2]), radius[1:2,i], color=colors[i], linestyle = styles[i], label = L"\beta="*string(β),linewidth=1)
    end
end

ylabel!(L"r/2\sqrt{p}")
xlabel!(L"\log_2(p)")
figure_path = joinpath(script_dir,figures_dir, "success_rates_new.pdf")
savefig(P, figure_path)
display(P)





