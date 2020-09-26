using Plots
using SpectralMethodFFT


function k_lin_k(k_u::KFunc, dt::Real)

    k_K = k_Kgen(k_u.config)
    xrange = k_u.config.xranges[1]
    len = xrange[2] - xrange[1]

    k_factor = exp(im * k_K^3 * dt / (len^3))
    return k_u * k_factor

end


k_nonlin_k(k_u::KFunc) = - k_u ⊗ (k_∂x_k(k_u))


function k_RK4_k(k_u::KFunc, dt::Real)

    k_k1 = dt * k_nonlin_k(k_u)
    k_k2 = dt * k_nonlin_k(k_lin_k(k_u + k_k1/2, dt/2))
    k_k3 = dt * k_nonlin_k(k_lin_k(k_u, dt/2) + k_k2/2)
    k_k4 = dt * k_nonlin_k(k_lin_k(k_lin_k(k_u, dt/2) + k_k3, dt/2))

    return k_lin_k(
        k_lin_k(k_u + k_k1/6, dt/2) + (k_k2 + k_k3)/3, dt/2
    ) + k_k4/6

end


function k_RK4_nstep_k(k_u::KFunc, dt::Real, nstep::Int)

    if nstep == 0
        return k_u
    elseif nstep > 0
        return k_RK4_k(k_RK4_nstep_k(k_u, dt, nstep-1), dt)
    end

end


function output(x_X::XFunc, k_u::KFunc, anim)

    x_u = x_k(k_u)
    plt = plot(x_X, x_u)
    frame(anim, plt)

end


function main()

    # *** configure FFT ***
    ngrids = (128,)
    xranges = ((0., 2.),)
    c = ConfigFFT(ngrids, xranges)

    # *** setting for time integration ***
    t_st = 0.
    t_ed = 0.1.
    nt = 100000
    dt = (t_ed - t_st) / nt
    nt_out = 1000
    out_period = nt ÷ nt_out

    # *** set initial condition ***
    x_X = x_Xgen(c)
    x_u = cos(π * x_X)
    k_u = k_x(x_u)

    # *** time integration ***
    anim = Animation()

    for it = 0:nt_out
        output(x_X, k_u, anim)
        if it < nt_out
            k_u = k_RK4_nstep_k(k_u, dt, out_period)
        end
    end

    gif(anim, "KdV_equation.gif", fps=10)

end

main()
