using LinearAlgebra, Plots
# x" - [2ωsin(l)]y' + (g/L)x = 0
# y" - [2ωsin(l)]x' + (g/L)y = 0

const ω = 7.29e-5   # rad*s^-1
const g = 9.81      # m*s^-2
const L = 10        # m
const l = 51.05     # Latitude: 51° 02' 60.00" (Gent, Belgie)
const α = 2ω*sin(l)
const β = g/L

IC = [1 0 0 0] # x(0) x'(0) y(0) y'(0)
dom = [0 100]


# x' = u; u' = f(x, y, t, u, v); y' = v; v' = g(x, y, t, u, v)
# p' = [x'; u'; y'; v']

f(t, p) = begin
    [p[2] ; α*p[3] .- β*p[1] ; p[4] ; α*p[1] .- β*p[3]]
end

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

tspan = dom[1]:0.1:dom[2]
u = RK4(f, [1 0 0 0], dom, .1)

plot(u[1,:][1:5:end], u[3,:][1:5:end], linestyle=:dot, legend=:outertopleft, label="Motion of Pendulum")
savefig("images/foucault_pendulum.pdf")

θ = atan.(u[4,:] ./ u[2, :])