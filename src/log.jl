
@views function subproblem(Q::AbstractMatrix, A₀::AbstractMatrix, β::Real, max_iter::Integer)
    """
    This routine is Algorithm 4.2 from the associated paper.
    The routine performs max_iter iterations to solve the subproblem.
    Input: Q an matrix from SO(2p), A₀ an initial p×p skew-symmetric matrix, β the parameter of the metric.
    Output: a p×p skew-symmetric matrix D.
    """
    if iszero(max_iter)
        return copy(A₀)
    end
    p = size(Q, 1) ÷ 2; iter = 0
    τ = 2β - 1
    if iszero(τ)
        return skewlog(Q)[1:p,1:p]/(2β)
    end

    R = similar(Q, 2p, 2p)
    A = similar(Q, p, p)
    E = similar(Q, p, p)
    D = copy(A₀)
    while  iter < max_iter
        R .= Q
        mul!(R[:, 1:p] , Q[:, 1:p], exp(τ * D))
        A = skewhermitian!(skewlog(R)[1:p, 1:p] /(2β))
        E = exp(τ * A)
        #Use R as memory space
        R[1:p, 1:p] .= A - D
        mul!(R[p+1:end, 1:p], R[1:p, 1:p], E')
        D .= A
        mul!(D, E, R[p+1:end, 1:p], τ, 1)
        iter += 1
    end
    return D #, iter - 1
end

@views function shootingsubproblem(V::AbstractMatrix{T}, A::AbstractMatrix{T}, β::Real, max_iter::Integer) where T
    """
    The subproblem can eventually be solved using a shooting method (not prefered).
    This routine performs max_iter shooting iterations to solve the subproblem.
    Input: V an matrix from SO(2p), A₀ an initial p×p skew-symmetric matrix, β the parameter of the metric.
    Output: a p×p skew-symmetric matrix D.
    """
    p = size(V, 1) ÷ 2
    if iszero(max_iter)
        return copy(A[1:p, 1:p]) ./ (2β)
    end
    m = 4; ϵ = 1e-12; iter = 0;  ν = 1
    t = LinRange(0, 1, m)
    
    D = zeros(p,p)
    Δ = zeros(2p, 2p)
    M  = zeros(2p, 2p)
    V₊ = zeros(2p, 2p)
    temp = zeros(Complex{T}, 2p, 2p)
    temp2 = zeros(Complex{T}, 2p, 2p)
    while ν > ϵ && iter < max_iter
        D .= A[1:p, 1:p]
        D .*= (1-2β)/ 2β 
        EA = eigen(skewhermitian!(A))
        ED = eigen(skewhermitian!(D))
        for i ∈ m:-1:1
            mul!(temp, EA.vectors, Diagonal(exp.(t[i].* EA.values)))
            mul!(temp2, temp, EA.vectors')
            mul!(temp[1:p, 1:p], ED.vectors, Diagonal(exp.(t[i] .* ED.values)))
            mul!(temp[p+1:end, 1:p], temp[1:p, 1:p], ED.vectors')
            M[:, 1:p] .= real.(temp2[:, 1:p])
            M[1:p, p+1:end] .= real.(temp[p+1:end, 1:p])
            mul!(V₊[:, 1:p] , M[:, 1:p], M[1:p, p+1:end])
            V₊[:, p+1:end] .= real(temp2[:, p+1:end])
            if i == m 
                Δ = V - V₊
                ν = norm(Δ)
            end
            mul!(M ,V₊', Δ)
            skewhermitian!(M)
            M .*= ν / norm(M)
            mul!(Δ , V₊ , M)
        end
        A[1:p, 1:p] ./= 2β
        A .+= M
        A[1:p, 1:p] .*= 2β
        iter += 1
    end
    D .= A[1:p, 1:p]
    D ./= 2β 
    return D
end

@views function initial_approximant(V::AbstractMatrix, β::Real)
    """
    Estimates the initial approximation Â₀ using the BCH series expansion. (See Section 4.5 of the paper)
    Input: V an orthogonal matrix of SO(2p), β the parameter of the metric.
    Output: the initial p×p skew-symmetric approximation Â₀ of A₀.
    """
    p = size(V, 1) ÷ 2
    S = similar(V, p, p)
    X = skewlog(V)
    X₂₁ = X[p+1:end, 1:p]
    mul!(S, X₂₁', X₂₁, (2β-1)/ 12, 0)
    S += I/2
    return lyapc(S, -X[1:p, 1:p])
end

@views function initial_basis(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}) where {T, METRIC}
    """
    Given two points on the Stiefel manifold S₁,S₂ ∈ St(n,p). Compute the initial matrix V from Algorithms  3.1 and 4.1.
    Input: S₁,S₂ ∈ St(n,p) (type StiefelVector).
    Output: The matrix V ∈ SO(2p).
    """
    n, p = size(S₁)
    #Allocate memory
    V = randn(2p, 2p)
    Q = similar(S₂.U.Q, n, p)
    temp = similar(S₂.U.Q, p, p)

    Q .= S₂.U.Q
    mul!(V[1:p,1:p], S₁.U.Q', S₂.U.Q)
    mul!(temp,S₁.U.Q', Q)
    mul!(Q, S₁.U.Q, temp, -1, 1)
    orthogonalize!(Q)
    lowrank = false
    for i ∈ 1:p
        #Complete basis Q if not full column rank
        if norm(Q[:, i]) < 0.5 #0.5 to avoid potential numerical issues of <1
            Q[:, i] .= randn(n)  
            lowrank = true
        end
    end
    if lowrank == true
        #Orthogonalize new completed basis
        mul!(V[p+1:end, 1:p], S₁.U.Q', Q)
        mul!(temp, S₁.U.Q', Q)
        mul!(Q, S₁.U.Q, temp, -1, 1)
        orthogonalize!(Q)
    end
    mul!(V[p+1:end, 1:p] , Q', S₂.U.Q)
    orthogonalize!(V) #security orthogonalization, theoretically not needed
    S = svd(V[p+1:end, p+1:end])

    #= Initialization suggested by Zimmermann 2017
    mul!(temp, S.Vt', S.U')
    temp2 = similar(temp, 2p, p)
    temp2 .= V[:, p+1:end]
    mul!(V[:, p+1:end], temp2, temp)
    =#

    #Alternative initialization with better performance
    temp2 = similar(temp, 2p, p)
    temp2 .= V[:, p+1:end]
    mul!(V[:, p+1:end], temp2, S.Vt')
    temp2 .= V[p+1:end, :]'
    mul!(V[p+1:end, :], S.U', temp2')
    Q₂ = copy(Q)
    mul!(Q, Q₂, S.U)
    
    if logabsdet(V)[2] < 0        #Using logabsdet may be more stable than det to ensure V ∈ SO(2p)
        V[:, end] .*= -1
    end
    return V, Q
end

@views function logβ(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}, β::Real, max_sub_iter::Integer) where {T,METRIC}
    """
    Computes the Riemannian logarithm Log(S₁,S₂) on St(n,p) equipped the β-metric.
    Performs pseudo-backward iterations of Algorithm 4.1 (paper).
    Input: S₁,S₂ two vectors on St(n,p), β > 0 the metric parameter 
    and  max_sub_iter the number of sub-iteration of Algorithm  4.2 (function subproblem) to peform.
    Output: The Riemannian logarithm Δ = S₁ * A + Q * B , Q, the residuals and the number of iterations performed.  
    NB: max_sub_iter=0 corresponds to the pure forward iteration.
    """
    #parameters
    ϵ = 1e-12; max_iter = 100; iter = 0
    p = size(S₁, 2)
    τ = 2β - 1
    norms = zeros(max_iter, 2)
    V, Q = initial_basis(S₁, S₂)

    #Pre-allocating memory
    temp = similar(V, 2p, p)
    logV= similar(V, 2p, 2p)
    Vₐ= similar(V, 2p, 2p)
    A = similar(V, p, p)
    Ahat = similar(V, p, p)
    S = similar(V, p, p)
    Φ = similar(V, p, p)
    E = similar(V, p, p)
    #End memory pre-allocation
    
    Ahat = initial_approximant(V, β)
    while iter < max_iter
        E = exp(skewhermitian!(τ * Ahat))
        mul!(Vₐ[:,1:p], V[:,1:p], E)
        Vₐ[:,p+1:end] .= V[:,p+1:end]
        logV = skewlog(Vₐ)
        A .= skewhermitian!(logV[1:p,1:p])
        A ./= 2β
        B = logV[p+1:end,1:p]
        C = logV[p + 1:end, p + 1:end]

        norms[iter+1, :] = [abs(1-2β) * norm(A-Ahat), norm(C)]
        if (norms[iter+1,1] + norms[iter+1,2]) < max(ϵ * max(norms[1,1]+norms[1,2],1), 2p * eps(T))
            iter += 1
            break
        end 
        mul!(S, B, B', 1/12, 0)
        S -= I/2
        Φ = exp(skewhermitian!(lyapc( S, -C)))
        mul!(temp, V[:,p+1:end], Φ)
        V[:,p+1:end] .= temp
        if β != 0.5
            Ahat = subproblem(V, A, β, max_sub_iter)  ##subproblem(V, A, β, max_sub_iter, ϵ/100) 
        else
            Ahat .= 0
        end
        iter+=1
    end
    return S₁.U.Q * A + Q * logV[p+1:end,1:p], Q, norms, iter
end

@views function logβmomentum(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}, β::Real) where {T,METRIC}
    """
    Computes the Riemannian logarithm Log(S₁,S₂) on St(n,p) equipped the β-metric.
    Performs accelerated forward iterations of Algorithm 4.1 (paper).
    Input: S₁,S₂ two vectors on St(n,p), β > 0 the metric parameter 
    and  max_sub_iter the number of sub-iteration of Algorithm  4.2 (function subproblem) to peform.
    Output: The Riemannian logarithm Δ = S₁ * A + Q * B , Q, norms (the residuals) and the number of iterations performed.  
    """
    #parameters
    ϵ = 1e-12; max_iter = 100; iter = 0
    n, p = size(S₁)
    τ = 2β - 1
    norms = zeros(max_iter, 2)
    V, Q = initial_basis(S₁, S₂)

    #Pre-allocating memory
    temp = similar(V, 2p, p)
    logV= similar(V, 2p, 2p)
    Vₐ= similar(V, 2p, 2p)
    A = similar(V, p, p)
    Ahat = similar(V, p, p)
    S = similar(V, p, p)
    Φ = similar(V, p, p)
    E = similar(V, p, p)
    #End memory pre-allocation

    Ahat = initial_approximant(V, β)
    while iter < max_iter
        E = exp(τ * Ahat)
        mul!(Vₐ[:,1:p], V[:,1:p], E)
        Vₐ[:,p+1:end] .= V[:,p+1:end]
        logV = skewlog(Vₐ)
        A .= skewhermitian!(logV[1:p,1:p])
        A ./= 2β
        B = logV[p+1:end,1:p]
        C = logV[p + 1:end, p + 1:end]

        norms[iter+1, :] = [abs(1-2β)*norm(A-Ahat), norm(C)]
        if (norms[iter+1,1]+norms[iter+1,2])< max(ϵ * max(norms[1,1]+norms[1,2],1), 2*p*eps(T))
            iter += 1
            break
        end 
        mul!(S, B, B', 1/12, 0)
        S -= I/2
        Φ = exp(skewhermitian!(lyapc( S, -C)))
        mul!(temp, V[:,p+1:end], Φ)
        V[:,p+1:end] .= temp
        if β == 1/2
            Ahat .= 0
        else
            E2 = exp(τ * A)
            # Use Vₐ as temp memory
            # Ahat = A + τ *  E2*(A - Ahat)*E2' 
            Vₐ[1:p, 1:p] .= τ * (A - Ahat) #- τ^2 * (A * Ahat - Ahat * A)/2
            mul!(Vₐ[p+1:end, 1:p], Vₐ[1:p, 1:p], E2')
            mul!(Ahat, E2, Vₐ[p+1:end, 1:p])
            Ahat .+= A 
        end
        iter+=1
    end
    return S₁.U.Q * A + Q * logV[p+1:end,1:p], Q, norms, iter
end


@views function pshooting(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}, β::Real, m::Integer) where {T,METRIC}
    """
    Computes the Riemannian logarithm Log(S₁,S₂) on St(n,p) equipped the β-metric.
    Implementation of Algorithm 2 from 
    Zimmermann, R., Hüper, K.: Computing the Riemannian Logarithm on the Stiefel Mani-
    fold: Metrics, Methods, and Performance. SIAM Journal on Matrix Analysis and Appli-
    cations 43(2), 953–980 (2022).
    Input: S₁,S₂ two vectors on St(n,p), β > 0 the metric parameter 
    and  m the discretization parameter.
    Output: The Riemannian logarithm Δ = S₁ * A + Q * B , Q and the number of iterations performed.  
    """
    t = LinRange(0, 1, m)
    max_iter = 1000
    n, p = size(S₁)
    ϵ = 1e-12

    #Pre-allocating memory for efficiency
    M = zeros(T, p, p)
    N = similar(M, p, p)
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
    temp2 = zeros(Complex{T}, 2p, 2p)
    temp2bis = similar(temp2, 2p, 2p)
    MN = similar(M, 2p, p)
    temp4 = similar(temp2, p, p)
    temp5 = similar(temp2, p, p)
    #End of memory pre-allocation

    mul!(M, S₁.U.Q', S₂.U.Q)
    V .= S₂.U.Q
    mul!(V, S₁.U.Q, M, -1, 1)
    Q = orthogonalize!(V)
    mul!(N, Q', S₂.U.Q)
    ν = sqrt(β * norm(M-I)^2 + norm(N)^2)
    A .= skewhermitian(M); R .= N
    γ = ν / sqrt(β * norm(A)^2+ norm(N)^2)
    A .*= γ ; R .*= γ
    iter = 0
    while ν > ϵ && iter<max_iter
        
        temp .= 0
        temp[1:p, 1:p] .= A
        temp[1:p, 1:p] .*= 2β
        temp[p+1:end, 1:p] .= R
        temp[1:p, p+1:end] .= -R'
        S .= skewhermitian(temp)
        E₁ = eigen!(S)
        E₂ = eigen((1-2β) * A)
        for j ∈ m:-1:1
            #if j !=1
            D₁ = Diagonal(exp.(t[j].* E₁.values))
            D₂ = Diagonal(exp.(t[j].* E₂.values))
            mul!(temp2bis, E₁.vectors, D₁)
            mul!(temp2, temp2bis, E₁.vectors')
            mul!(temp4, E₂.vectors, D₂)
            mul!(temp5, temp4, E₂.vectors')
            mul!(MN, real.(temp2[:,1:p]), real.(temp5))
            #else
            #    MN[1:p, :] .= Matrix(1.0I, p, p)
            #end
            if j == m
                Aˢ .= MN[1:p, :] ;  Aˢ.-= M
                Rˢ .= MN[p+1:end, :] ; Rˢ.-= N
                ν  = sqrt(β * norm(Aˢ)^2 + norm(Rˢ)^2)
            else
                mul!(Sym2, MN[1:p, :]', Aˢ)
                mul!(Sym2, MN[p+1:end, :]', Rˢ, 1, 1)
                #take symmetric part S = S2 + S2' / 2
                Sym1 .= Sym2
                Sym1.+= Sym2'
                Sym1 ./= 2
                mul!(Aˢ, MN[1:p, :], Sym1, -1, 1)
                mul!(Rˢ, MN[p+1:end, :], Sym1, -1, 1)
                l = sqrt(β * norm(Aˢ)^2+norm(Rˢ)^2)
                if l > ϵ
                    Aˢ .*= ν/l
                    Rˢ .*= ν/l
                else
                    Aˢ .= 0 ; Rˢ .= 0
                end
            end
        end
        A .-= Aˢ
        R .-= Rˢ 
        iter +=1 
    end
    return S₁.U.Q * A + Q * R, Q, iter
end