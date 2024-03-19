"""
This file contains a structure to deal with the Stiefel manifold.
"""

mutable struct StiefelVector{T<:Number, METRIC, X<:Orthonormal{<:T}}
    U::X 
    Θ::Union{Nothing, Orthonormal{T}} #Θ is the orthonormal complement of U
    #METRIC: real number > 0. METRIC=1 for euclidean, METRIC=1/2 for canonical metric
    function StiefelVector{T, METRIC, X}(U) where {T, METRIC, X}
        Base.require_one_based_indexing(U)
        new{T, METRIC, X}(U, nothing)
    end
end

"""
    StiefelVector(U, METRIC)

Type representing a point of the Stiefel Manifold.
Let S be a StiefelVector, S.U must be orthonormal, i.e `U^TU = I`.
S.Θ is the orthonormal complement of S.U and can be added/build to the structure
using `addcomplement` or `buildcomplement`.
"""
function StiefelVector(U::AbstractMatrix, METRIC = 1)
    isorthonormal(U) || throw(ArgumentError("Columns of U must be orthonormal"))
    return StiefelVector{eltype(U), METRIC, Orthonormal{<:eltype(U)}}(Orthonormal(U, false))
end

StiefelVector(U::Orthonormal, METRIC = 1) = StiefelVector{eltype(U), METRIC, typeof(U)}(U)
Base.size(S::StiefelVector) = size(S.U)
Base.size(S::StiefelVector, n) = size(S.U, n)
dims(S::StiefelVector) = size(S)
Base.copy(S::StiefelVector{T, METRIC}) where {T, METRIC} = StiefelVector(copy(S.U), METRIC)

"""
    buildcomplement(S)

Construct the orthonormal complement of S.U using
Gram-Schmidt orthogonalization.
"""
@views function buildcomplement(U::AbstractMatrix{T}) where T
    m, n = size(U)
    n₂ = m - n
    Θ = randn(T, m, n₂)
    for i ∈ 1:n
        for j ∈ 1:n₂
            α = dot(U[:,i], Θ[:, j] )
            Θ[:, j] .-= α.*U[:,i]
        end
    end

    for i ∈ 1:n₂
        nm = norm(Θ[:, i])
        Θ[:, i] ./= nm
        for j ∈ i+1:n₂
            α = dot(Θ[:,i], Θ[:, j] )
            Θ[:, j] .-= α.*Θ[:,i]
        end
    end

    temp = similar(U, m, m)
    temp[:,1:n] .= U
    temp[:,n+1:end] .= Θ
    if det(temp) < 0
        Θ[:,n₂] .*= -1
    end
    return Θ
end

@views function buildcomplement(S::StiefelVector{T}) where T
    S.Θ = Orthonormal(buildcomplement(S.U.Q), false)
    return S.Θ
end

@views function addcomplement(S::StiefelVector, Q::Orthonormal) 
    iszero(transpose(S.U) * Q) || throw(ArgumentError("Q must be the orthonormal complement of S.U"))
    S.Θ = Q
end

"""
    addcomplement(S, Q)

Add the orthonormal complement Q of S.U to S where S is a StiefelVector
"""
addcomplement(S::StiefelVector, Q::AbstractMatrix) = addcomplement(S, Orthonormal(Q))

@views function projection!(W::AbstractMatrix, S::AbstractMatrix)
    temp₁ = similar(W, size(S, 2), size(W, 2))
    mul!(temp₁, transpose(S), W)
    mul!(W, S, temp₁, -1, 1)
    skewhermitian!(temp₁)
    mul!(W, S, temp₁, 1, 1)
    return W
end

@views function projection!(W::AbstractMatrix, S::StiefelVector)
    return projection!(W, S.U.Q)
end

"""
    projection(W, S)
Returns the projection of any matrix W m×n on the
tangent space to the StiefelVector S ∈ St(m, n). 
"""
projection(W::AbstractMatrix, S::StiefelVector) = projection!(copy(W), S::StiefelVector)

mutable struct TangentVector{T<:Number, METRIC, Am<:SkewHermitian{T}, Bm<:Union{Nothing, AbstractMatrix{T}}}
    A::Am
    B::Bm
    function TangentVector{T, METRIC, Am, Bm}(A, B) where {T, METRIC, Am, Bm}
        Base.require_one_based_indexing(A)
        B === nothing || Base.require_one_based_indexing(B)
        new{T, METRIC, Am, Bm}(A, B)
    end
end

TangentVector(A::SkewHermitian{T}, B::Union{Nothing, AbstractMatrix{T}}, METRIC = 1) where T = TangentVector{T, METRIC, typeof(A), typeof(B)}(A, B)
TangentVector(A::StridedMatrix{T}, B::AbstractMatrix{T}, METRIC = 1) where {T} = TangentVector{T, METRIC, SkewHermitian{T},typeof(B)}(skewhermitian(A), B)

Base.:-(Δ::TangentVector{T,METRIC}) where {T,METRIC}= TangentVector(-Δ.A, -Δ.B, METRIC)

for f in (:+, :-)
    @eval begin
        Base.$f(Δ::TangentVector{T,METRIC}, Γ::TangentVector{T,METRIC}) where {T,METRIC}= TangentVector($f(Δ.A,Γ.A), $f(Δ.B,Γ.B), METRIC)
   end
end
# Scaling:
for op in (:*, :/, :\)
    if op in (:*, :/)
        @eval Base.$op(Δ::TangentVector{T,METRIC}, x::Real) where {T,METRIC} = TangentVector($op(Δ.A, x), $op(Δ.B, x), METRIC)
    end
    if op in (:*, :\)
        @eval Base.$op(x::Real, Δ::TangentVector{T,METRIC}) where {T,METRIC} = TangentVector($op(Δ.A, x), $op(Δ.B, x), METRIC)
    end
end

function gettangent(Δ::TangentVector, S::StiefelVector)
    if Δ.B === nothing 
        temp = similar(Δ.A, size(S, 1), size(Δ.A, 2))
        mul!(temp, S.U, Δ.A)
        return temp
    else
        S.Θ === nothing && throw(ArgumentError("Complementary space of S must be specified"))
        temp₁ = similar(Δ.A, size(S, 1), size(Δ.A, 2))
        mul!(temp₁, S.U, Δ.A)
        mul!(temp₁, S.Θ, Δ.B, 1, 1)
        return temp₁
    end
end

function LA.dot(Δ::TangentVector{T, METRIC}, Γ::TangentVector{T, METRIC}) where {T, METRIC}
    Σ = 0
    n = size(Δ.A, 2)
    for i ∈ 1:n
        Σ += dot(Δ.A[:, i], Γ.A[:, i])  
    end
    #display(Σ)
    Σ *= METRIC
    n = size(Δ.B, 2)
    for i ∈ 1:n
        Σ += dot(Δ.B[:, i], Γ.B[:, i])
    end
    return Σ
end

function LA.norm(Δ::TangentVector)
    return sqrt(dot(Δ, Δ))
end 

@views function dotβ(Δ::AbstractMatrix{T}, Γ::AbstractMatrix{T}, S::StiefelVector{T, METRIC}, β::Real) where {T, METRIC}
    n, p = size(Δ)
    temp = similar(Δ, n,  p)
    res = similar(Δ, p,  p)
    mul!(res, S.U.Q', Γ)
    temp .= Δ
    mul!(temp, S.U.Q, res, -(1-β), 1)
    mul!(res, Δ', temp)
    return tr(res)
end

@views function normβ(Δ::AbstractMatrix{T}, S::StiefelVector{T, METRIC}, β::Real) where {T, METRIC}
    return sqrt(dotβ(Δ, Δ, S, β))
end

@views function gramschmidt!(A::AbstractMatrix)
    n = size(A, 2)
    R = zeros(eltype(A), n, n)
    iszero(A) && return Matrix(randorthonormal(eltype(A), size(A, 1), n)), R
    for i ∈ 1:n
        R[i, i] = norm(A[:, i])
        A[:, i] ./= R[i,i]
        for j ∈ i+1:n
            R[i, j] = dot(A[:, i], A[:, j])
            A[:, j] .-= R[i,j].*A[:, i]
        end
    end
    return A, R
end

gramschmidt(A::AbstractMatrix) = gramschmidt!(copy(A))

function dist(A::AbstractMatrix, B::AbstractMatrix)
    E = eigvals(A' * B)
    logE = zeros(length(E))
    for i ∈ eachindex(E)
        logE[i] = acos(real(E[i])/abs(E[i]))
    end
    return norm(logE)/sqrt(2)
end

@views function myangle(c::Real, s::Real)
    """
    Computes the angle of a Givens rotation.
    Input: cosine c and sine s of an angle θ
    Output: θ ∈ (-π, π]
    """
    if s > 1
        @warn "WARNING: sine expected < 1. sine was set to 1"
        s = 1
    end

    if s < -1
        @warn "WARNING: sine expected > -1. sine was set to -1"
        s = -1
    end

    if c < 0
        return (s < 0 ? -π - asin(s) : π - asin(s))
    else
        return asin(s)
    end
end

@views function multiplybylog(V::AbstractMatrix{T}, S::AbstractMatrix{T}) where T
    """
    Performs a sparse matrix-matrix multiplication between V and log(S).
    Input: Square matrix V and block diagonal matrix S from a Real Schur form.
    Output: V * log(S) where log(S) is the principal logarithm of S.
    """
    n = size(V, 1)
    M = zeros(T, n, n)
    i = 1; mem = 0 #mem allows to treat the -1 eigenvalues
    while i <= n
        if i<n && !iszero(S[i+1, i])
            θ = myangle(S[i,i], S[i+1, i])
            M[:, i]   .=  θ * V[:, i+1]
            M[:, i+1] .= -θ * V[:, i]
            i += 2
        elseif 1 + S[i,i] < eps(T)
            if !iszero(mem)
                M[:, i]   .=   π * V[:, mem]
                M[:, mem] .= - π * V[:, i]
                mem = 0
            else
                mem = i
            end
            i += 1
        else
            i += 1
        end
    end
    return M
end

@views function skewlog(Q::AbstractMatrix{T}) where T
    """
    Computes the real skew-symmetric logarithm of a real orthogonal matrix Q.
    Input: Real Orthogonal matrix Q.
    Output: Real Skew-symmetric logarithm of Q.
    """
    S = schur(Q)
    #The log of S.T is sparse, multpiplybylog takes advantage of the sparsity
    return skewhermitian!(multiplybylog(S.Z, S.T) * S.Z')
end




    

