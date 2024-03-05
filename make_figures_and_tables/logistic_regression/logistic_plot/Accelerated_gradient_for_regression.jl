using LinearAlgebra

function f(θ::AbstractVector, X::AbstractMatrix, y::AbstractVector)
    z = -X'θ
    return sum(log.(1 .+ exp.(z)) .- ((1 .- y) .* z))
end

function ∇f(θ::AbstractVector, X::AbstractMatrix, y::AbstractVector)
    z = -X * θ
    return X' * (1 ./(1 .+ exp.(z)) .- y)
end

#Unused gradient method
function GradientMethod(θ₀::AbstractVector, X::AbstractMatrix, y::AbstractVector, L::Real, max_iter::Integer)
    θ = copy(θ₀); ϵ = 1e-8
    iter = 1
    norms = zeros(max_iter)  
    g = ∇f(θ, X, y)
    norms[1] = norm(g)
    while  norm(g) > ϵ && iter < max_iter
        θ -= g ./ L
        g = ∇f(θ, X, y)
        iter +=1
        norms[iter] = norm(g)
    end
    return θ, norms, iter
end

function AccGradientMethod(θ₀::AbstractVector, X::AbstractMatrix, y::AbstractVector, L::Real, max_iter::Integer)
    θ = copy(θ₀); ϵ = 1e-8
    iter = 1
    norms = zeros(max_iter)  
    g = ∇f(θ, X, y)
    norms[1] = norm(g)
    t = 1
    η = copy(θ)
    while  norm(g) > ϵ && iter < max_iter
        θ₂ = copy(θ)  
        θ = η - ∇f(η, X, y) ./ L
        tₒ = t
        t = (1 + sqrt(1 + 4 * tₒ * tₒ)) / 2
        η = θ + (tₒ - 1) / t * (θ - θ₂)
        g = ∇f(θ, X, y)
        iter +=1
        norms[iter] = norm(g)
    end
    return θ, norms, iter
end