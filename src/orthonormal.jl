struct Orthonormal{T<:Number, X<:AbstractMatrix{<:T}} <: AbstractMatrix{T}
    Q::X
    function Orthonormal{T, X}(Q) where {T, X}
        Base.require_one_based_indexing(Q)
        new{T, X}(Q)
    end
end
"""
    Orthonormal(Q, check)

Type for orthogonal matrices `Q^TQ = I`.
When it is initialized, check says if orthogonality must be verified.
"""

function Orthonormal(Q::AbstractMatrix, check = true)
    if check
        isorthonormal(Q) || throw(ArgumentError("Columns of Q must be orthonormal"))
    end
    return Orthonormal{eltype(Q), typeof(Q)}(Q)
end

function isorthonormal(Q::AbstractMatrix)
    transpose(Q) * Q ≈ I || return false
    return true
end

Base.@propagate_inbounds Base.getindex(Q::Orthonormal, i::Integer, j::Integer) = Q.Q[i,j]
Base.@propagate_inbounds function Base.setindex!(Q::Orthonormal, v, i::Integer, j::Integer)
    Q.Q[i,j] = v
    return v
end
#Base.similar(Q::Orthonormal, :: Type{T}) where {T} = Orthonormal(similar(parent(Q), T) .= I)
Base.Matrix(Q::Orthonormal) = Matrix(Q.Q)
Base.Array(Q::Orthonormal) = Array(Q.Q)
Base.parent(Q::Orthonormal) = Q.Q
Base.copy(Q::Orthonormal) = Orthonormal(copy(Q.Q), false)

function Base.copyto!(dest::Orthonormal, src::Orthonormal)
    copyto!(dest.Q, src. Q)
end

function Base.copyto!(dest::Orthonormal, src::AbstractMatrix)
    isorthonormal(src) || throw(ArgumentError("matrix to copy must be orthogonal"))
    copyto!(dest.Q, src)
end

Base.copyto!(dest::AbstractMatrix, src::Orthonormal) = copyto!(dest, src.Q)
Base.size(Q::Orthonormal) = size(Q.Q)
Base.size(Q::Orthonormal, n) = size(Q.Q, n)

Base.isreal(Q::Orthonormal) = isreal(Q.Q)
Base.transpose(Q::Orthonormal) = transpose(Q.Q)
Base.adjoint(Q::Orthonormal) = Orthonormal(Q.Q', false)
Base.real(Q::Orthonormal{<:Real}) = Q
Base.real(Q::Orthonormal) = real(Q.Q)
Base.imag(Q::Orthonormal) = imag(Q.Q)

#Base.conj(Q::Orthonormal) = Orthonormal(conj(Q.Q))
#Base.conj!(Q::Orthonormal) = Orthonormal(conj!(Q.Q))
LA.tr(Q::Orthonormal) = tr(Q.Q)

LA.tril!(Q::Orthonormal) = tril!(Q.Q)
LA.tril(Q::Orthonormal)  = tril!(copy(Q))
LA.triu!(Q::Orthonormal) = triu!(Q.Q)
LA.triu(Q::Orthonormal)  = triu!(copy(Q))
LA.tril!(Q::Orthonormal,k::Integer) = tril!(Q.Q,k)
LA.tril(Q::Orthonormal,k::Integer)  = tril!(copy(Q),k)
LA.triu!(Q::Orthonormal,k::Integer) = triu!(Q.Q,k)
LA.triu(Q::Orthonormal,k::Integer)  = triu!(copy(Q),k)

LA.dot(A::Orthonormal, B::Orthonormal) = dot(A.Q, B.Q)
LA.dot(x::AbstractVector, Q::Orthonormal, y::AbstractVector) = dot(x, Q, y)
Base.:-(Q::Orthonormal) = Orthonormal(-Q.Q, false)
for f ∈ (:+, :-)
    @eval begin
        Base.$f(A::Orthonormal, B::Orthonormal) = $f(A.Q, B.Q)
   end
end

## Matvec
LA.mul!(y::StridedVecOrMat, A::Orthonormal, x::StridedVecOrMat, α::Number, β::Number) =
    LA.mul!(y, A.Q, x, α, β)
LA.mul!(y::StridedVecOrMat, A::Orthonormal, x::StridedVecOrMat) =
    LA.mul!(y, A.Q, x)
LA.mul!(y::StridedVecOrMat, A::Orthonormal, B::Orthonormal, α::Number, β::Number) =
    LA.mul!(y, A.Q, B.Q, α, β)
LA.mul!(y::StridedVecOrMat, A::Orthonormal, B::Orthonormal) =
    LA.mul!(y, A.Q, B.Q)
LA.mul!(y::StridedVecOrMat, A::StridedMatrix, B::Orthonormal, α::Number, β::Number) =
    LA.mul!(y, A, B.Q, α, β)
LA.mul!(y::StridedVecOrMat, A::StridedMatrix, B::Orthonormal) =
    LA.mul!(y, A, B.Q)

"""
    getorthonormal(T, m, n)

Returns an orthonormal matrix of size m×n
using modified Gram-Schmidt on a randomly 
generated matrix.
"""
randorthonormal(T::Type, m::Integer, n::Integer) = orthogonalize!(randn(T, m, n))

@views function orthogonalize!(Q::AbstractMatrix{T}) where T
    n = size(Q, 2)
    @inbounds(for i ∈ 1:n
        nm = norm(Q[:, i])
        if nm > 10 * eps(T)
            Q[:, i] ./= nm
            for j ∈ i+1:n
                α = dot(Q[:, i], Q[:, j])
                Q[:, j] .-= α.*Q[:, i]
            end
        end
    end)
    return Orthonormal(Q, false)
end

"""
    orthogonalize(Q)

Build an orthonormal matrix of size m×n, n<=m 
from Q using Gram-Schmidt algorithm.
"""
orthogonalize(Q::AbstractMatrix) = orthogonalize!(copy(Q))

