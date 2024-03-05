#This file serves as a bin. It contains old unused routines.

@views  function logAlg1(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}) where {T,METRIC}
    ϵ = eps(T) * 100
    N = ceil(Int, norm(S₁.U.Q-S₂.U.Q)) * 5
    δt = 1 / (N - 1)
    t = Array(0:δt:1)
    γ = norm(S₁.U.Q - S₂.U.Q)
    Δ = projection(S₂.U.Q, S₁)
    Δ ./= norm(Δ)
    Δ .*= γ
    Δˢ = similar(Δ, size(Δ, 1), size(Δ, 2))
    itermax = 100; iter = 0
    while γ > ϵ  && iter < itermax
        Δˢ .= Matrix(exp(S₁, Δ))
        Δˢ .-= S₂.U
        γ = norm(Δˢ)
        for j = N:-1:1
            Δˢ = projection(Δˢ, StiefelVector(exp(S₁, Δ.*t[j])))
            Δˢ ./= norm(Δˢ)
            Δˢ .*= γ
        end
        Δ .-= Δˢ
        iter +=1
        #display(Δ)
    end
    return Δ
end

@views function logAlg4(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}) where {T,METRIC}
    p = size(S₁)[2]
    M = similar(S₁.U.Q, p, p)
    mul!(M, S₁.U', S₂.U)
    Q, N = gramschmidt(S₂.U - S₁.U * M)
    V = similar(S₁.U.Q, 2p, 2p)
    V[1:p, 1:p] .= M
    V[p+1:end,1:p] .= N
    V[:,p+1:end] .= buildcomplement(V[:,1:p])
    Svd = svd(V[p+1:end,p+1:end])
    V[p+1:end,p+1:end] = V[p+1:end,p+1:end] * Svd.Vt' * Svd.U'
    V[1:p,p+1:end] = V[1:p,p+1:end] * Svd.Vt' * Svd.U'
    if det(V)<0
        V[:,end] .*= (-1)
    end
    temp = similar(V, 2p, p)
    ϵ = 1e-8
    max_iter = 100
    iter = 0
    logV= similar(V, 2p, 2p)
    iterates = zeros(max_iter)
    while iter < max_iter
        logV .= logleast(V)
        iterates[iter+1] = norm(logV/2)
        if norm(logV[p+1:end,p+1:end]) < ϵ
            break
        end
        S = (logV[p+1:end,1:p]*logV[p+1:end,1:p]')/12-I/2
        Φ = exp(skewhermitian(lyapc(S,-logV[p+1:end,p+1:end])))
        mul!(temp, V[:,p+1:end], Φ)
        V[:,p+1:end] .= temp
        iter+=1
    end
    #display(norm(logV/2))
    return S₁.U.Q * skewhermitian!(logV[1:p,1:p]) + Q * logV[p+1:end,1:p],Q, iterates,iter
end

@views function logAlg4SVD(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}) where {T,METRIC}
    n, p = size(S₁)
    temp₁= similar(S₁.U.Q, n, n)
    temp₂= similar(S₂.U.Q, n, n)
    
    Svd = svd(Matrix(S₁.Θ'*S₂.Θ))
    temp₁[:,1:p] .= S₁.U
    mul!(temp₁[:,p+1:end], S₁.Θ, Svd.U)
    temp₂[:,1:p] .= S₂.U
    mul!(temp₂[:,p+1:end], S₂.Θ, Svd.Vt')
    
    #check ∈SO(n)
    if det(temp₁) < 0
        temp₁[:,end] .*= -1
    end
    if det(temp₂) < 0 
        temp₂[:,end] .*= -1
    end
    Q = similar(temp₁, n, n)
    mul!(Q, temp₁', temp₂)
    V = similar(S₁.U.Q, 2p, 2p)
    V[1:p, 1:p] .= Q[1:p, 1:p]
    V[p+1:end,1:p] .= Q[n-p+1:n, 1:p]
    V[1:p,p+1:end] .= Q[1:p, n-p+1:n]
    V[p+1:end,p+1:end] .= Q[n-p+1:n,n-p+1:n]
    temp = similar(V, 2p, p)
    ϵ = 1e-4
    max_iter = 100
    iter = 0
    logV= similar(V, 2p, 2p)
    iterates = zeros(max_iter)
    while iter < max_iter
        logV = logleast(V)
        iterates[iter+1] = norm(logV/2)
        if norm(logV[p+1:end,p+1:end]) < ϵ
            break
        end
        S = (logV[p+1:end,1:p]*logV[p+1:end,1:p]')/12-I/2
        Φ = exp(skewhermitian(lyapc(S,-logV[p+1:end,p+1:end])))
        mul!(temp, V[:,p+1:end], Φ)
        V[:,p+1:end] .= temp
        iter+=1
    end
    return S₁.U.Q * skewhermitian!(logV[1:p,1:p]) + temp₁[:,n-p+1:end] * logV[p+1:end,1:p],temp₁[:,n-p+1:end], iterates, iter
end
@views function shootingsubproblem_direct(V::AbstractMatrix{T}, A::AbstractMatrix{T}, β::Real, max_iter::Integer) where T
    #shooting method with direct projection of the error on the tangent space to the initial point.
    p = size(V, 1) ÷ 2
    if iszero(max_iter)
        return copy(A[1:p, 1:p]) ./ (2β)
    end
    ϵ = 1e-12; iter = 0;  ν = 1
    D = zeros(p,p)#skewhermitian!(zeros(p, p))
    Δ = zeros(2p, 2p)
    M  = zeros(2p, 2p)
    V₊ = zeros(2p, 2p)
    while ν > ϵ && iter < max_iter
        D .= A[1:p, 1:p]
        D .*= (1-2β)/ 2β 
        EA = exp(skewhermitian!(A))
        ED = exp(skewhermitian!(D))
        mul!(V₊[:, 1:p] , EA[:, 1:p], ED)
        V₊[:, p+1:end] .= EA[:, p+1:end]
        Δ .= V - V₊
        ν = norm(Δ)
        mul!(M , V₊', Δ)
        skewhermitian!(M)
        mul!(Δ , V₊ , M)
        M .= Δ
        skewhermitian!(M)
        M .*= ν / norm(M)
        A[1:p, 1:p] ./= 2β
        A .+= M
        A[1:p, 1:p] .*= 2β
        iter += 1    
    end
    D .= A[1:p, 1:p]
    D ./= 2β 
    return D
end