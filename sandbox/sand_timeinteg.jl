function RK4_split(lin, nonlin, x, t, dt)

    # lin(t_new, t_old, x)
    # nonlin(x, t)

    k1 = dt * nonlin(x, t)
    k2 = dt * nonlin(lin(t+dt/2, t, x+k1/2), t+dt/2)
    k3 = dt * nonlin(lin(t+dt/2, t, x)+k2/2, t+dt/2)
    k4 = dt * nonlin(lin(t+dt, t+dt/2, lin(t+dt/2, t, x)+k3), t+dt)

    return lin(t+dt, t+dt/2, lin(t+dt/2, t, x+k1/6)+(k2+k3)/2) + k4/6

end
