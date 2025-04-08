using LinearAlgebra, SkewLinearAlgebra, Plots


@views function logcanonical(S₁::StiefelVector{T, METRIC}, S₂::StiefelVector{T, METRIC}, param::Number) where {T,METRIC}
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
    norms = zeros(max_iter)
    V, Q = initial_basis(S₁, S₂)

    #Pre-allocating memory
    temp = similar(V, 2p, p)
    logV= similar(V, 2p, 2p)
    C = similar(V, p, p)
    #Vₐ= similar(V, 2p, 2p)
    A = similar(V, p, p)
    #Ahat = similar(V, p, p)
    S = similar(V, p, p)
    Φ = similar(V, p, p)
    #E = similar(V, p, p)
    #End memory pre-allocation
    
    #Ahat = initial_approximant(V, β)
    while iter < max_iter
        #E = exp(skewhermitian!(τ * Ahat))
        #mul!(Vₐ[:,1:p], V[:,1:p], E)
        #Vₐ[:,p+1:end] .= V[:,p+1:end]
        #logV = skewlog(Vₐ)
        logV = skewlog(V)
        A .= skewhermitian!(logV[1:p,1:p])
        #A ./= 2β
        B = logV[p+1:end,1:p]
        if iter > 0
            Cold = copy(C)
        end
        C .= logV[p + 1:end, p + 1:end]

        norms[iter+1, :] = [norm(C)]
        if (norms[iter+1])< max(ϵ * max(norms[1,1],1), 2*p*eps(T))
            iter += 1
            break
        end 
        mul!(S, B, B', 1/12, 0)
        S -= I/2
        if param == 0
            Φ = exp(-C) #exp(skewhermitian!(lyapc( S, -C)))
        elseif param == 1
            Φ = exp(skewhermitian!(lyapc( S, -C)))
        else
            if iter >= 0 
                c = maximum(C); J = randn(p,p); J = J - J'; J .*= c/norm(J)
                Φ = exp(-J)
            else
                Φ = exp( -C)
            end
        end
        mul!(temp, V[:,p+1:end], Φ)
        V[:,p+1:end] .= temp
        #=
        if β != 0.5
            Ahat = subproblem(V, A, β, max_sub_iter)  ##subproblem(V, A, β, max_sub_iter, ϵ/100) 
        else
            Ahat .= 0
        end
        =#
        iter+=1
    end
    return S₁.U.Q * A + Q * logV[p+1:end,1:p], Q, norms, iter
end
#=
n = 40; p =10
δ = 0.05
U₁ = randorthonormal(Float64, n, p)
U₂ = Orthonormal(exp(skewhermitian!(randn(n,n)*δ))*U₁.Q, false)#m.randorthonormal(Float64, n, p)
display(norm(U₁-U₂)/sqrt(p)/2)
S₁ = StiefelVector(U₁, 0)
S₂ = StiefelVector(U₂, 0)
norms0, iter0 = logcanonical(S₁, S₂, 0)[3:4]
norms1, iter1 = logcanonical(S₁, S₂, 2)[3:4]
P = plot(1:iter0, norms0[1:iter0], yscale=:log)
scatter!(1:iter0, norms0[1:iter0])
plot!(1:iter1, norms1[1:iter1])
scatter!(1:iter1, norms1[1:iter1])
display(P)
=#