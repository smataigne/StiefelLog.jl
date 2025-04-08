using LinearAlgebra, Plots

n = 10
A = randn(n, n); A = (A - A')/2
B = randn(n, n); B = (B - B')/2
A ./= norm(A); B ./= norm(B)
ts = [1, 0.1, 0.01, 0.001, 0.0001]
s = length(ts)
prec1 = zeros(s)
prec2 = zeros(s)
δ = 0.5
for (i, t) ∈ enumerate(ts)
    true1 = t * (B - A)#exp(-t * A / 2) * log(exp(-t* A) * exp(t * B)) * exp(t * A / 2)
    est1 = log(exp(-t * A) * exp(t * B)) #exp(-t * A / 2) * (t * (B - A))* exp(-t * A / 2)
    est2 = log(exp(-t * A * (1- δ)) * exp(t * B) * exp(-t * A  * δ))  #exp(-t * A) * (t * (B - A))* exp(t * A)
    prec1[i] = norm(est1 - true1)
    prec2[i] = norm(est2 - true1)
end
P = plot()
plot!(ts, prec1, label ="1", xscale=:log, yscale=:log)
plot!(ts, prec2, label ="2")
display(P)