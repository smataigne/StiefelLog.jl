#=
Algorithm from
D. Bryner, "Endpoint Geodesics on the Stiefel 
Manifold Embedded in Euclidean Space", SIAM ,
https://doi.org/10.1137/16M1103099, 2017.
=#

@views function myexp(Δ::AbstractArray, U::AbstractArray)
    n, p = size(U)
    U₂ = similar(U, n, p)
    U_end = similar(U, n, p)
    temp = similar(U, n, n)
    A = similar(U, p, p)
    mul!(temp, Δ, U', 2, 0)
    mul!(U₂, exp(skewhermitian!(temp)), U)
    mul!(A, U', Δ, -1, 0)
    mul!(U_end, U₂, exp(skewhermitian!(A)))
    orthogonalize!(U_end)
    return U_end
end

@views function subroutine_1!(U::AbstractArray)
    Svd = svd(U)
    mul!(U, Svd.U, Svd.Vt)
    return U
end

@views function subroutine_3!(Δ::AbstractArray, U::AbstractArray)
    l1 = norm(Δ)
    projection!(Δ, U)
    l2 = norm(Δ)
    if l2 > eps(eltype(U))
        Δ .*= (l1/l2)
    else
        Δ .= 0
    end
    return Δ
end

@views function path_initialization!(U₁::AbstractMatrix, U₂::AbstractMatrix, T::AbstractVector, γ::AbstractArray)
    for i ∈ eachindex(T)
        @. γ[:, :, i] = (1-T[i]) * U₁ + T[i] * U₂
        subroutine_1!(γ[:, :, i])
    end
    return γ
end

@views function compute_δγ!(δγ::AbstractArray, γ::AbstractArray, T::AbstractVector)
    N = length(T) - 1 
    for i ∈ eachindex(T)
        if i == 1
            @. δγ[:, :, i] = (γ[:, :, i+1] - γ[:, :, i])/(T[i+1]-T[i])
        elseif i == N + 1
            @. δγ[:, :, i] = (γ[:, :, i] - γ[:, :, i-1])/(T[i]-T[i-1])
        else
            @. δγ[:, :, i] = (γ[:, :, i+1] - γ[:, :, i-1])/(T[i+1]-T[i-1])
        end
        projection!(δγ[:, :, i], γ[:, :, i])
    end
    return δγ
end

@views function compute_covint!(δγ::AbstractArray, γ::AbstractArray, u::AbstractArray)
    n, p, N = size(γ) 
    u[:, :, 1] .= 0
    u_trans = similar(u, n, p)
    for i ∈ 2:N
        u_trans .= u[:, :, i-1]
        subroutine_3!(u_trans, γ[:, :, i])
        @. u[:, :, i]= δγ[:, :, i] / (N-1) + u_trans
    end
    return u
end

@views function backward_parallel_trans!(u::AbstractArray, u₂::AbstractArray, γ::AbstractArray)
    N = size(γ, 3)
    u₂[:, :, N] .= u[:, :, N]
    for i = N-1:-1:1
        u₂[:, :, i] .= u₂[:, :, i+1]
        subroutine_3!(u₂[:, :, i], γ[:, :, i])
    end
    return u₂
end

@views function gradient_of_E!(u::AbstractArray, u₂::AbstractArray, ω::AbstractArray, T::AbstractVector)
    N = size(u, 3)
    for i ∈ eachindex(T)
        @. ω[:, :, i] = u[:, :, i] - u₂[:, :, i] * T[i]
    end
    return ω
end

@views function path_update!(γ::AbstractArray, ω::AbstractArray, T::AbstractVector, δ::Real)
    n, p, N = size(γ)
    temp = similar(γ, n, p)
    for i ∈ eachindex(T)
        @. temp = -δ * ω[:, :, i]
        γ[:, :, i] = myexp(temp, γ[:, :, i])
    end
    return γ
end

function path_straightening(U₁::AbstractMatrix, U₂::AbstractMatrix)
    N = ceil(Int, norm(U₁-U₂)) * 30
    n, p = size(U₁)
    τ = 0.01
    δ = 0.5
    iter = 0
    itermax = 100
    #Memory allocation
    T =  LinRange(0,1,N)#cos(2*(N:-1:1) * π/(2*(N+1)))/2 + 1/2
    γ = similar(U₁, n, p, N)
    δγ = similar(U₁, n, p, N)
    ω = similar(U₁, n, p, N)
    u = similar(U₁, n, p, N)
    u₂ = similar(U₁, n, p, N)

    nω = τ + 1 
    path_initialization!(U₁, U₂, T, γ)
    while nω > τ && iter < itermax 
        iter += 1
        compute_δγ!(δγ, γ, T)
        compute_covint!(δγ, γ, u)
        backward_parallel_trans!(u, u₂, γ)
        gradient_of_E!(u, u₂, ω, T)
        nω = norm(ω)
        path_update!(γ, ω, T, δ)
    end
    return δγ[:, :, 1]
end
#=
n = 5
p = 3
U₁ = Matrix(m.randorthonormal(Float64, n, p))
U₂ = Matrix(m.randorthonormal(Float64, n, p))
Δ,  mem = path_straightening(U₁, U₂)
output  = myexp(Δ , U₁)
display(norm(U₂  - output))
display(U₂)
display(output)
plot(1:length(mem), mem, xaxis=:log, yaxis=:log, framestyle =:box,gridlinewidth=1,gridalpha=0.5,minorgrid=true,minorgridalpha=0.2,label = false)
xlabel!("Iteration [/]")
ylabel!("Norm of the gradient of the vector field")
=#
