include("../../generalFFT.jl")
using Plots


function init_cond!(xy_u::XFunc{T,2} where T)

    c = xy_u.config
    xy_X = xy_Xgen(c)
    xy_Y = xy_Ygen(c)
    xy_u.vals .= sinc(sqrt(xy_X^2 + xy_Y^2)).vals

end

function time_develop(kl_u::KFunc{T,2} where T, dt::Float64)

    config = kl_u.config
    kl_K = kl_Kgen(config)
    kl_L = kl_Lgen(config)
    xranges = config.xranges
    xlens = (x -> -(.-(x...))).(xranges)

    return kl_u * exp(- (
        (2π*kl_K/xlens[1])^2 + (2π*kl_L/xlens[2])^2
    )*dt)

end


function main()

    ngrids = (128,128)
    xranges = ((-10., 10.),(-10., 10.))

    c = ConfigFFT(ngrids, xranges)
    xy_u = XFunc(undef, c)
    init_cond!(xy_u)
    kl_u = kl_xy(xy_u)

    x, y = c.Xcoords[1][:,1], c.Xcoords[2][1,:]

    t_st = 0.
    t_ed = 2.
    nt = 100
    dt = (t_ed - t_st) / nt

    anim = Animation()

    for it = 0:nt
        xy_u = xy_kl(kl_u)
        plt = plot(x, y, xy_u,
            zlims=(0., 0.5), st=:wireframe)
        frame(anim, plt)
        if it < nt
            kl_u = time_develop(kl_u, dt)
        end
    end
    gif(anim, "heat_diffusion_2d.gif", fps = 10)

end

main()
