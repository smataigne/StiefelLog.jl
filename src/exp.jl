@views function expβ(S::StiefelVector{T, METRIC}, Δ::TangentVector{T, METRIC}) where {T, METRIC}
    
    Δ₂ = gettangent(Δ, S)
    temp₁ = similar(Δ₂, size(Δ₂, 1), size(S, 1))
    temp₂ = similar(Δ₂, size(S, 1), size(Δ.A, 2))
    Exp = similar(Δ₂, size(S, 1), size(Δ.A, 2))
    mul!(temp₁, Δ₂, S.U')
    E₁ = exp(skewhermitian!(temp₁) .* 2)
    E₂ = exp((1-2*METRIC) * Δ.A)
    mul!(temp₂, S.U, E₂)
    mul!(Exp, E₁, temp₂)
    return Orthonormal(Exp, false)
end

@views function expβ(S::StiefelVector{T, METRIC}, Δ::AbstractMatrix) where {T, METRIC}
    Δ₂ = Δ
    temp₁ = similar(Δ₂, size(Δ₂, 1), size(S, 1))
    temp₂ = similar(Δ₂, size(S, 1), size(S, 2))
    A = similar(Δ₂, size(S, 2), size(S, 2))
    Exp = similar(Δ₂, size(S, 1), size(S, 2))
    mul!(temp₁, Δ₂, transpose(S.U))
    E₁ = exp(skewhermitian!(temp₁.*2))
    mul!(A, transpose(S.U) , Δ, -1, 0)
    E₂ = exp(skewhermitian!(A))
    mul!(temp₂, S.U, E₂)
    mul!(Exp, E₁, temp₂)
    return Orthonormal(Exp, false)
end

@views function expcanonical(S::StiefelVector, Δ::TangentVector)
    Δ₂ = gettangent(Δ, S)
    m, n = size(S, 1), size(S, 2)
    temp₁ = similar(Δ₂, n, n)
    temp₂ = similar(Δ₂, m, n)
    temp₃ = zeros(eltype(Δ₂), 2*n, 2*n)
    temp₄ = similar(Δ₂, m, 2*n)
    mul!(temp₁, S.U.Q', Δ₂)
    mul!(Δ₂, S.U.Q, temp₁, -1, 1)
    Q = orthogonalize(Δ₂)
    R = Q'Δ₂
    temp₃[1:n, 1:n] .= Δ.A
    temp₃[n+1:end,1:n] .= R.*2
    E₁ = exp(skewhermitian!(temp₃))
    temp₄[:, 1:n] .= S.U
    temp₄[:, n+1:end] .= Q
    mul!(temp₂, temp₄, E₁[:, 1:n])
    return Orthonormal(temp₂, false)
end

@views function expcanonical(S::StiefelVector, Δ::AbstractMatrix)
    return expcanonical(S.U.Q, Δ)
end

@views function expcanonical(S::AbstractMatrix, Δ::AbstractMatrix)
    
    n, p = size(S)
    Q = similar(S, n, p)
    B = similar(S, p, p)
    M = similar(S, 2p, 2p)
    E = similar(S, n, 2p)
    Res = similar(S, n, p)
    mul!(M[1:p, 1:p], S', Δ)
    Q .= Δ
    mul!(Q, S, M[1:p, 1:p], -1, 1)
    orthogonalize!(Q)
    mul!(B, Q', Δ)
    M[p+1:end, 1:p]  .= B .* 2
    skewhermitian!(M)
    #M = (M-M')/2
    E[:, 1:p] .= S 
    E[:, p+1:end] .= Q
    mul!(Res, E, exp(M)[:, 1:p])
    #=
    n, p = size(S)
    A = S'Δ
    Res = exp(-3/2 * S*A *S' +Δ * S'-S*Δ') * S=#
    return Orthonormal(Res, false)
end

function Base.exp(S::StiefelVector{T, METRIC}, Δ::TangentVector{T, METRIC}) where {T, METRIC}
    if iszero(METRIC - 1/2)
        return expcanonical(S, Δ)
    end
    return expβ(S, Δ)
end

function Base.exp(S::StiefelVector{T, METRIC}, Δ::AbstractMatrix{T}) where {T, METRIC}
    if iszero(METRIC-1/2)
        return expcanonical(S, Δ)
    end
    return expβ(S, Δ)
end