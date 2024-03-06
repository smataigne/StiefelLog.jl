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
@views function dummypshooting(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}, β::Real, m::Integer) where {T,METRIC}
    t = LinRange(0, 1, m)
    max_iter = 1000
    n, p = size(S₁)
    ϵ = 1e-12

    #Pre-allocating memory for efficiency
    M = zeros(T, p, p)
    V = similar(M, n, p)
    A = similar(M, p, p)
    R = similar(M, p, p)
    Aˢ= similar(M, p, p)
    Rˢ= similar(M, p, p)
    R = similar(M, p, p)
    S = similar(M, 2p, 2p)
    Sym1 = similar(M, p, p)
    Sym2 = similar(M, p, p)
    temp = similar(M, 2p, 2p)
    MN = similar(M, 2p, p, m)
    #End of memory pre-allocation

    mul!(M, S₁.U.Q', S₂.U.Q)
    V .= S₂.U.Q
    mul!(V, S₁.U.Q, M,-1, 1)
    Q, N = gramschmidt!(V)
    ν = sqrt(norm(M-I)^2 + norm(N)^2) #||U-̃U ||_F
    A .= (M-M')/2; R .= N
    γ = ν / sqrt(β * norm(A)^2+ norm(N)^2)
    A .*= γ ; R .*= γ
    iter = 0
    while ν > ϵ && iter<max_iter
        
        temp .= 0
        temp[1:p, 1:p] .= A
        temp[1:p, 1:p] .*= 2β
        temp[p+1:end, 1:p] .= R
        temp[1:p, p+1:end] .= -R'
        S .= temp
        for j ∈ 1:m
            MN[:, :, j] = (exp(t[j]*S)[:,1:p])*exp(t[j]*(1-2β) * A)
        end
        Aˢ .= MN[1:p, :, m] ;  Aˢ.-= M
        Rˢ .= MN[p+1:end, :, m] ; Rˢ.-= N
        ν  = sqrt(norm(Aˢ)^2 + norm(Rˢ)^2)
        for j ∈ m:-1:1
            mul!(Sym2, MN[1:p, :, j]', Aˢ)
            mul!(Sym2, MN[p+1:end, :, j]', Rˢ, 1, 1)
            #take symmetric part S = S2 + S2' / 2
            Sym1 .= Sym2
            Sym1.+= Sym2'
            Sym1 ./= 2
            mul!(Aˢ, MN[1:p, :, j], Sym1, -1, 1)
            mul!(Rˢ, MN[p+1:end, :, j], Sym1, -1, 1)
            l = sqrt(norm(Aˢ)^2+norm(Rˢ)^2)
            if l > ϵ
                Aˢ .*= ν/l
                Rˢ .*= ν/l
            else
                Aˢ .= 0 ; Rˢ .= 0
            end
        end
        A .-= Aˢ
        R .-= Rˢ 
        iter +=1 
    end
    return S₁.U.Q * A + Q * R, Q, iter
end

@views function Base.log(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}) where {T,METRIC}
    if iszero(METRIC)
        return path_straightening(S₁.U, S₂.U) 
    else
        return logAlg4(S₁,S₂)
    end
end