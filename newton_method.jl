using LinearAlgebra, Plots
x = -1:0.1:1
y = -1:0.1:1

f(x, y) = begin
    x^2 + y^2 - 1
end

fx(x, y) = begin
    2x
end

fy(x, y) = begin
    2y
end

h(x, y) = begin
    x^4 - y^4 + x*y
end

hx(x, y) = begin
    4x^3 + y
end

hy(x, y) = begin
    -4y^3 + x
end


function jacobian(fx, fy, gx, gy)
    return [fx(x[1], x[2]) fy(x[1], x[2]) ; hx(x[1], x[2]) hy(x[1], x[2])]
end

function newton(f, h, n_iter, tol=1e-6)
    #Startin value
    x = [1; 1]
    #Vector function F
    J = jacobian(fx, fy, hx, hy)
    for i = 1:n_iter
        x = hcat(x, x[:,end] + J\[f(x[1, end], x[2, end]) ; h(x[1, end], x[2, end])])
        if norm(x[:,end] - x[:, end-1]) <= tol
            println("Convergence after ", i, " iterations.")
            println("Best estimate for x: ", x[:, end])
            plot(1:length(x[1,:]), x[1,:], label="x val")
            p1 = plot!(1:length(x[1,:]), x[2,:], label="y val")
            savefig("/images/newton_method.pdf")
            return p1
        end
    end
    println("No convergence after ", n_iter, " iterations.")
    println("Best estimate for x: ", x[:, end])
    plot(1:length(x[1,:]), x[1,:], label="x val")
    p2 = plot!(1:length(x[1,:]), x[2,:], label="y val")
    savefig("images/newton_method.pdf")
    return p2
end

newton(f, h, 20, 1e-5)
