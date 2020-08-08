using LinearAlgebra, Plots

const k = 14.2e2
const mo, mc = 16, 12
const Mo, Mc = 16*1.6605e-27, 12*1.6605e-27
const c = 3e8

A(ω) = begin
    [(k-mo*ω^2) -k 0 ; -k (2k - mc*ω^2) -k ; 0 -k (k-mo*ω^2)]
end

Ω = [0 sqrt(k/Mo) sqrt(k/Mo)*sqrt(1 + 2*(Mo/Mc))]

B = [A(i) for i in  Ω]

λ(ω) = begin
    2π*c/ω
end

Λ = [λ(i) for i in Ω]

for i = 2:3
    println("λ", i, " = ", Λ[i]*1e6, "μm")
end

for i = 1:3
    display(eigvecs(A(Ω[i])))
end