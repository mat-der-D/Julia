include("../../generalFFT.jl")
using Plots



function init_cond!(x_u::XFunc{T,1} where T)

    c = x_u.config
    x_X = x_Xgen(c)
    x_u.vals .= sinc(5x_X).vals

end

function time_develop(k_u::KFunc{T,1} where T, dt::Float64)

    config = k_u.config
    k_K = k_Kgen(k_u.config)
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

    x_X = c.Xcoords[1]
    # plot(x_X, x_u.vals)

    k_u = time_develop(k_u, 50.0)
    x_u = x_k(k_u)

    plot(x_X, x_u.vals)

end

main()
