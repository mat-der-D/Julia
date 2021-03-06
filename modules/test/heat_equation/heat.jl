include("../../generalFFT.jl")
using Plots


function init_cond!(x_u::XFunc{T,1} where T)

    c = x_u.config
    x_X = x_Xgen(c)
    x_u.vals .= sinc(5x_X).vals

end

function time_develop(k_u::KFunc{T,1} where T, dt::Float64)

    config = k_u.config
    k_K = k_Kgen(config)
    xrange = config.xranges[1]
    xlen = xrange[2] - xrange[1]

    return k_u * exp( - (2π*k_K/xlen)^2 * dt)

end


function main()

    ngrids = (128,)
    xranges = ((-10., 10.),)

    c = ConfigFFT(ngrids, xranges)
    x_u = XFunc(undef, c)
    init_cond!(x_u)
    k_u = k_x(x_u)

    Xcoord = c.Xcoords[1]
    # plot(x_X, x_u.vals)

    t_st = 0.
    t_ed = 10.
    nt = 100
    dt = (t_ed - t_st) / nt

    anim = Animation()

    for it = 0:nt
        x_u = x_k(k_u)
        plt = plot(Xcoord, x_u.vals,
                ylims=(0., 0.25), labels="heat")
        frame(anim, plt)
        if it < nt
            k_u = time_develop(k_u, dt)
        end
    end
    gif(anim, "heat_diffusion.gif", fps = 10)

end

main()
