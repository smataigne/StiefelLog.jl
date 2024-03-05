using LinearAlgebra, SkewLinearAlgebra, XLSX
import .Manifold as m


n = 64; p = 32           #Choose n and p, n > p
N = 1000                 #Choose the number of samples N

distances = zeros(N)
βs = 0.5:0.1:1           #Values of β that are considered
k = length(βs)
succ = zeros(N, k)       #Remember success of failure, 1 or 0.
for i ∈ 1:N
    U₁ = m.randorthonormal(Float64, n, p)
    if i > (3 * N / 4)   #Sample 25% from antipodal antipodal point
        δ = 8 * (i - 3 * N / 4) / N
        U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n) * δ)) * (-U₁.Q), false)
    else 
        δ = (i / N) * 4 / 3
        U₂ = m.Orthonormal(exp(skewhermitian!(randn(n,n) * δ)) * U₁.Q, false)
    end
    S₁ = m.StiefelVector(U₁, 0)
    S₂ = m.StiefelVector(U₂, 0)
    distances[i] = norm(U₁.Q-U₂.Q)/ (2*sqrt(p))  #Compute relative Frobenius distance
    for j ∈ 1:k
        β = βs[j]
        try 
            Δ, Q, norms, Niter = m.logβ(S₁, S₂, β, 2)
            if Niter < 98                        #Returned successfully
                succ[i, j] = 1
            end
        catch
            succ[i,j] = 0
        end
    end
end

#The data is put in Excel file to be used by Logistic_regression.jl
columns = Vector()
push!(columns, distances)
for i ∈ 1:k
    push!(columns, succ[:,i])
end
columns2 = Vector()
push!(columns2,[n; p; N; length(βs)])
columns3 = Vector()
push!(columns3,Vector(βs))
labels = zeros(k+1)
labels[1] = 0
labels[2:k+1] = βs
XLSX.openxlsx("./logistic_data/Logistic_data_st"*string(p)*"_1000.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "new_sheet")
    XLSX.writetable!(sheet, columns2, ["n,p,N,size"], anchor_cell=XLSX.CellRef("A1"))
    XLSX.writetable!(sheet, columns3, ["βs"], anchor_cell=XLSX.CellRef("B1"))
    XLSX.writetable!(sheet, columns, labels, anchor_cell=XLSX.CellRef("C1"))   
end
