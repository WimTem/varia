using LinearAlgebra, Plots

#solve: 
# y''' + y' - sin(y) = 0
# y(0) = 0
# y'(0) = 0
# y(1) = 10

# y = x1
# y' = x2
# y" = x3
# y''' = x4 = -x2 + sin(x1)
f(t, y) = begin
    [y[2] ; y[3] ; -y[2] + sin(y[1])]
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

tspan = 0:0.1:1
#Guess initial BC too low
y0 = [0 0 10]

u1 = RK4(f, y0, [0 1], 0.1)

#Check y(1)
u1[1, end]  # Too low

y1 = [0 0 30]

u2 = RK4(f, y1, [0 1], 0.1)

#Check y(1)
u2[1, end] # Too high

function bisection(a, b)
    return 0.5(a+b)
end

function shooting_method(f, y0, dom, h, target, bc_low, bc_high, n_iter=10, tol=1e-3)
    result = []
    for i = 1:n_iter
        u = RK4(f, hcat(y0, bc_low), dom, h)
        result = push!(result, u[1,:])
        v = RK4(f, hcat(y0, bc_high), dom, h)
        result = push!(result, v[1,:])
        
        w = RK4(f, hcat(y0, bisection(bc_low, bc_high)), dom, h)

        if abs(w[1, end] - target) < tol
            println("Convergence after ", i, " iterations")
            println("Guess for y''(0): ", (bc_low+bc_high)/2)
            push!(result, w[1, :])
            #Niet alles plotten
            #ymax = sol
            return result;
        end

        if w[1, end] > target && abs(w[1, end] - target) > tol
            bc_high = bisection(bc_low, bc_high)
        else
            bc_low = bisection(bc_low, bc_high)
        end
        println(round(u[1,end], digits=2), " : ", round(v[1,end], digits=2))
    end
    return println("No convergence after ", n_iter, " iterations.")
end

result = shooting_method(f, [0 0], [0 1], 0.01, 10, 10, 30, 100, 1e-4)
length(result)
plot([0:0.01:1 for i in 1:3], result[1:3], legend=:outertopleft, linestyle=:dash)
plot!(0:0.01:1, result[end], legend=:outertopleft, lw=3, label="Sol")
savefig("images/shooting_method.pdf")