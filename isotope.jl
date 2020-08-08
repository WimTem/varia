using LinearAlgebra, Plots

# Numeric integration
RK4(f, y0, dom, h) = begin
    x = dom[1]:h:dom[2]
    n = length(x)
    m = length(y0)

    y = zeros(m, n)
    y[:, 1] = y0

    h2 = h/2

    for i = 1:1:n-1
        k1 = f(x[i], y[:, i])
        k2 = f(x[i] .+ h2, y[:, i] .+ h2*k1)
        k3 = f(x[i] .+ h2, y[:, i] .+ h2*k2)
        k4 = f(x[i] .+ h, y[:, i] .+ h*k3)

        y[:, i+1] = y[:, i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    end
    return y
end

const k1, k2, k3, k4 = 5.78e-3, 2.87e-5, 2.09e-5, 1.1e-12 # s^-1

f(t, p) = begin
    [-k1*p[1] ; k1*p[1] - k2*p[2] ; k2*p[2] - k3*p[3] ; k3*p[3] - k4*p[4] ; k4*p[4]]
end

dom1 = [0 300]
dom2 = [0 1e4]

IC = [6.033e23 0 0 0 0]

u1 = RK4(f, IC, dom1, 10)
u2 = RK4(f, IC, dom2, 10)

tspan1 = dom1[1]:10:dom1[2]
tspan2 = dom2[1]:10:dom2[2]

p1 = plot(tspan1, u1', yaxis=:log, ylims=(1, 1e24), label=["Te" "I" "Xe" "Cs" "Ba"], xlabel="t [s]", ylabel="c [n/V]", lw=2)

p2 = plot(u2', xaxis=:log, yaxis=:log, ylims=(1, 1e24), label=["Te" "I" "Xe" "Cs" "Ba"], xlabel="t [s]", ylabel="c [n/V]", lw=2)

plot(p1, p2, layout=(2,1))
savefig("images/isotopic_decay.pdf")