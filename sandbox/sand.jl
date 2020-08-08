function f(a, dt)

    y = 1.0
    for n = 100:-1:2
        y = 1. - dt/n * a * y
    end
    return y * dt

end