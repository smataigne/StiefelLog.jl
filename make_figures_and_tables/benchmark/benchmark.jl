using LinearAlgebra, SkewLinearAlgebra
import .Manifold as m

n = 100; p = 50
betas = Vector(0.3:0.1:1)
m2 = length(betas)
times = zeros(6, m2)
lim = 10
δ = 0.1
for k ∈ 1:10
    U₁ = m.randorthonormal(Float64, n, p)
    U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n)*δ)) * U₁.Q, false)#m.randorthonormal(Float64, n, p)
    S₁ = m.StiefelVector(U₁, 0)
    S₂ = m.StiefelVector(U₂, 0)
    
    for j ∈ 1:m2
        β = betas[j]
        temp = minimum(@elapsed m.pshooting(S₁, S₂, β, 3) for _ ∈ 1:lim)
        times[1,j] = (times[2, j] * (k-1) + temp)/k
        temp = minimum(@elapsed m.pshooting(S₁, S₂, β, 5) for _ ∈ 1:lim)
        times[2,j] = (times[3, j] * (k-1) + temp)/k
        temp    = minimum(@elapsed m.logβmomentum(S₁, S₂, β) for _ ∈ 1:lim)
        times[3,j] = (times[4, j] * (k-1) + temp)/k
        temp    = minimum(@elapsed m.logβ(S₁, S₂, β, 0) for _ ∈ 1:lim)
        times[4,j] = (times[5, j] * (k-1) + temp)/k
        temp    = minimum(@elapsed m.logβ(S₁, S₂, β, 1) for _ ∈ 1:lim)
        times[5,j] = (times[6, j] * (k-1) + temp)/k
        temp    = minimum(@elapsed m.logβ(S₁, S₂, β, 2) for _ ∈ 1:lim)
        times[6,j] = (times[7, j] * (k-1) + temp)/k
    end
end
display(times)
