include("../../../FourierTrans/FFT2D_funcClass/FFT2D_funcClass.jl")

# ***** Settings of the system *****
struct HirutaFlow_settings

    # *** grid points ***
    nt::Int64
    # *** range ***
    t_st::Float64
    t_ed::Float64
    dt::Float64
    # *** physical params ***
    Re_inv::Float64	  # 1/Re
    κ::Float64
    Re0_inv::Float64	  # 1/Re0
    U_y::Float64
    n_force::Int64	  # "n" in exterior force term

end


function set_HirutaFlow(nt::Int64,
			t_st::Float64, t_ed::Float64,
			Re_inv::Float64,
			κ::Float64,
			U_y::Float64,
			n_force::Int64)::HirutaFlow_settings

    dt = (t_ed - t_st) / nt
    Re0_inv = κ^2 * Re_inv

    return HirutaFlow_settings(nt,
			       t_st, t_ed,
			       dt,
			       Re_inv, κ, Re0_inv,
			       U_y, n_force)

end

# ***** time tools *****
struct TimeIntegrationTools

    lk_addFactor::lk_Func
    lk_expFactor::lk_Func

end

function set_TimeIntegrationTools(
			s::HirutaFlow_settings,
			c::FFT2D_config)

    lk_K = lk_Func(c.lk_K, c)
    lk_L = lk_Func(c.lk_L, c)
    lk_alp = (s.Re0_inv + im * lk_L * s.U_y
              + s.Re_inv * (lk_K*lk_K + lk_L*lk_L))
    lk_expFactor = exp(- s.dt * lk_alp)

    yx_Y = yx_Func(c.yx_Y, c)
    yx_ncosny = s.n_force * cos(s.n_force * yx_Y)
    lk_ncosny = lk_yx(yx_ncosny)
    lk_addFactor = lk_Func(lk_ncosny / lk_alp, c)

    return TimeIntegrationTools(lk_addFactor, lk_expFactor)

end

function lk_linear_next_lk(lk_func::lk_Func,
	                   tit::TimeIntegrationTools)::lk_Func

    return ( (lk_func - tit.lk_addFactor) * tit.lk_expFactor
    	     + tit.lk_addFactor )

end


function lk_Jacobi_lk_lk(lk_f::lk_Func, lk_g::lk_Func)::lk_Func

    yx_fx = yx_lk(lk_Dx_lk(lk_f))
    yx_fy = yx_lk(lk_Dy_lk(lk_f))
    yx_gx = yx_lk(lk_Dx_lk(lk_g))
    yx_gy = yx_lk(lk_Dy_lk(lk_g))

    return lk_yx(yx_fx * yx_gy - yx_fy * yx_gx)

end

# ***** test program *****
function test_HirutaFlow()

    # *** parameter setting ***
    # -- Space --
    nx = 16384	    # # of grid points
    ny = 64	    # # of grid points
    xmin = 0.0	    # leftmost(lowest) value of x
    xmax = 512pi    # rightmost(highest) value of x
    ymin = 0.0	    # leftmost(lowest) value of y
    ymax = 2pi	    # rightmost(highest) value of y
    # -- Time --
    nt = 100	    # # of time steps
    t_st = 0.	    # start time
    t_ed = 0.01	    # end time
    # -- Physical Params --
    Re_inv = 1./240
    κ = 0. # sqrt(30)
    U_y = 0.5
    n_force = 4

    c = configure_FFT2D(nx, ny,
		        xmin, xmax, ymin, ymax)
    s = set_HirutaFlow(nt,
		       t_st, t_ed,
		       Re_inv, κ, U_y, n_force)

    yx_X = yx_Func(c.yx_X, c)
    yx_Y = yx_Func(c.yx_Y, c)

end