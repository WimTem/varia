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
