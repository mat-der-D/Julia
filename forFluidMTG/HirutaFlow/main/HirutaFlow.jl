include("../../../FourierTrans/FFT2D_funcClass/FFT2D_funcClass.jl")

# **********************************
# ***** Settings of the system *****
# **********************************
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


# **********************
# ***** time tools *****
# **********************
struct TimeIntegrationTools

    lk_addFactor::lk_Func
    lk_expFactor::lk_Func

end


function set_TimeIntegrationTools(
			s::HirutaFlow_settings,
			c::FFT2D_config;
			dtFactor=1.0)

    dt = s.dt * dtFactor

    lk_K = lk_Func(c.lk_K, c)
    lk_L = lk_Func(c.lk_L, c)
    lk_alp = (s.Re0_inv + 1.0im * lk_L * s.U_y
              + s.Re_inv * (lk_K*lk_K + lk_L*lk_L))
    lk_expFactor = exp(- dt * lk_alp)

    yx_Y = yx_Func(c.yx_Y, c)
    yx_ncosny = s.n_force * cos(s.n_force * yx_Y)
    lk_ncosny = lk_yx(yx_ncosny)
    lk_addFactor = lk_ncosny * dt * (
        1.0 - dt/2.0 * lk_alp * (
	    1.0 - dt/3.0 * lk_alp * (
	        1.0 - dt/4.0 * lk_alp * (
		    1.0 - dt/5.0 * lk_alp * (
		        1.0 - dt/6.0 * lk_alp * (
			    1.0 - dt/7.0 * lk_alp * (
			        1.0 - dt/8.0 * lk_alp * (
				    1.0 - dt/9.0 * lk_alp * (
				        1.0 - dt/10.0 * lk_alp
				    )
				)
			    )
			)
		    )
		)
	    )
	)
    )

    return TimeIntegrationTools(lk_addFactor, lk_expFactor)

end


function lk_linear_next_lk(lk_func::lk_Func,
	                   tit::TimeIntegrationTools)::lk_Func

    return lk_func * tit.lk_expFactor - tit.lk_addFactor

end


function lk_Jacobi_lk_lk(lk_f::lk_Func, lk_g::lk_Func)::lk_Func

    yx_fx = yx_lk(lk_Dx_lk(lk_f))
    yx_fy = yx_lk(lk_Dy_lk(lk_f))
    yx_gx = yx_lk(lk_Dx_lk(lk_g))
    yx_gy = yx_lk(lk_Dy_lk(lk_g))

    return lk_yx(yx_fx * yx_gy - yx_fy * yx_gx)

end


function lk_TimeDevelopHeun_lk(
			lk_ω::lk_Func,
			s::HirutaFlow_settings,
			tit::TimeIntegrationTools)::lk_Func

    dt = s.dt
    lk_ψ = lk_stream_vorticity_lk(lk_ω)

    # STEP1,2
    lk_ω2 = lk_linear_next_lk(
    			lk_ω + dt * lk_Jacobi_lk_lk(lk_ψ, lk_ω),
    			tit)

    # STEP3
    lk_ψ2 = lk_stream_vorticity_lk(lk_ω2)
    lk_ω3 = lk_Jacobi_lk_lk(lk_ψ2, lk_ω2)

    # STEP4,5
    lk_ω5 = lk_linear_next_lk(
			lk_ω + dt/2 * lk_Jacobi_lk_lk(lk_ψ, lk_ω),
			tit)

    return lk_ω5 + dt/2 * lk_ω3

end


function lk_TimeDevelopRK4_lk(
			lk_ω::lk_Func,
			s::HirutaFlow_settings,
			tit::TimeIntegrationTools)::lk_Func

    # Note: initialize tit as dtFactor=0.5

    dt = s.dt
    lk_ω0 = lk_Func(lk_ω.vals, lk_ω.config)

    # Step1
    lk_ψ0 = lk_stream_vorticity_lk(lk_ω0)
    lk_k1 = dt * lk_Jacobi_lk_lk(lk_ψ0, lk_ω0)

    # Step2
    lk_ω1 = lk_linear_next_lk(lk_ω0 + lk_k1 / 2.0, tit)
    lk_ψ1 = lk_stream_vorticity_lk(lk_ω1)
    lk_k2 = dt * lk_Jacobi_lk_lk(lk_ψ1, lk_ω1)

    # Step3
    lk_ω2 = lk_linear_next_lk(lk_ω0, tit) + lk_k2 / 2.0
    lk_ψ2 = lk_stream_vorticity_lk(lk_ω2)
    lk_k3 = dt * lk_Jacobi_lk_lk(lk_ψ2, lk_ω2)

    # Step4
    lk_ω3 = lk_linear_next_lk(lk_linear_next_lk(lk_ω0, tit) + lk_k3, tit)
    lk_ψ3 = lk_stream_vorticity_lk(lk_ω3)
    lk_k4 = dt * lk_Jacobi_lk_lk(lk_ψ3, lk_ω3)

    lk_ωNext = (
        lk_linear_next_lk(
	    lk_linear_next_lk(
	        lk_ω0 + lk_k1 / 6.0,
	        tit
	    )
	    + (lk_k2 + lk_k3) / 3.0,
	    tit
	)
	+ lk_k4 / 6.0
    )

    return lk_ωNext

end






# *****************************
# ***** Basic Fluid tools *****
# *****************************
function lk_stream_vorticity_lk(lk_ω::lk_Func)::lk_Func

    c = lk_ω.config
    lk_K = lk_Func(c.lk_K, c)
    lk_L = lk_Func(c.lk_L, c)
    lk_KKLL = lk_K*lk_K + lk_L*lk_L
    lk_KKLL.vals[1, 1] = 1.

    lk_ψ = lk_Func(lk_ω.vals, lk_ω.config)
    lk_ψ.vals[1, 1] = 0.

    return lk_ψ / lk_KKLL

end


function lk_vorticity_stream_lk(lk_ψ::lk_Func)::lk_Func

    return - ( lk_Dx_lk(lk_Dx_lk(lk_ψ))
    	     + lk_Dy_lk(lk_Dy_lk(lk_ψ)) )

end


function lk_lk_velocity_stream_lk(
			lk_ψ::lk_Func)::Tuple{lk_Func, lk_Func}

    return lk_Dy_lk(lk_ψ), - lk_Dx_lk(lk_ψ)

end


# ************************
# ***** test program *****
# ************************
function test_HirutaFlow()

    # *** parameter setting ***
    # -- Space --
    aspect = 4		 # aspect ratio
    nx = 64*aspect    	 # # of grid points
    ny = 64	      	 # # of grid points
    xmin = 0.0	      	 # minimum value of x
    xmax = 2pi*aspect    # maximum value of x
    ymin = 0.0	    	 # minimum value of y
    ymax = 2pi	    	 # maximum value of y
    # -- Time --
    nt = 25000		# # of time steps
    nt_output = 25000 	# output every nt_output
    t_st = 0.	    	# start time
    t_ed = 50.0      	# end time
    # -- Physical Params --
    Re_inv = 1.0 / 8.0
    κ = 0.0 # sqrt(30)
    U_y = 0.0
    n_force = 4

    c = configure_FFT2D(nx, ny,
		        xmin, xmax, ymin, ymax)
    s = set_HirutaFlow(nt,
		       t_st, t_ed,
		       Re_inv, κ, U_y, n_force)

    # for Heun
    # tit = set_TimeIntegrationTools(s, c)
    # for RK4
    tit = set_TimeIntegrationTools(s, c, dtFactor=0.5)

    yx_X = yx_Func(c.yx_X, c)
    yx_Y = yx_Func(c.yx_Y, c)

    open("testcase.dat", "w") do f
        lk_ω = lk_Func_undef(c)
        set_initial_condition(lk_ω)
	println("t = ", 0 * s.dt)
	println(f, "# dt=", s.dt)
        output_flowField(f, lk_ω)
	lk_ω_old = lk_Func(lk_ω.vals, lk_ω.config)
	
        for it = 1:s.nt
            # lk_ω = lk_TimeDevelopHeun_lk(lk_ω, s, tit)
	    lk_ω = lk_TimeDevelopRK4_lk(lk_ω, s, tit)
	    if it % 500 == 0
	        println("t = ", it * s.dt)
	    end
            if it % nt_output == 0
	        # println("t = ", it * s.dt)
                output_flowField(f, lk_ω)

            end
        end
    end

end

function output_flowField(f::IOStream, lk_ω::lk_Func)

    c = lk_ω.config
    yx_ω = yx_lk(lk_ω)
    for ix = 1:c.nx
        for iy = 1:c.ny
	    println(f, c.yx_X[iy, ix], " ",
	    	       c.yx_Y[iy, ix], " ",
		       yx_ω.vals[iy, ix])
	end
	println(f)
    end
    println(f)

end


function set_initial_condition(lk_ω::lk_Func)
    c = lk_ω.config
    lk_ω.vals = rand(c.ny, c.nx)
end


function test_ncosny()

    # *** parameter setting ***
    # -- Space --
    aspect = 4		 # aspect ratio
    nx = 64*aspect    	 # # of grid points
    ny = 64	      	 # # of grid points
    xmin = 0.0	      	 # minimum value of x
    xmax = 2pi*aspect    # maximum value of x
    ymin = 0.0	    	 # minimum value of y
    ymax = 2pi	    	 # maximum value of y
    # -- Time --
    nt = 25000		# # of time steps
    nt_output = 25000 	# output every nt_output
    t_st = 0.	    	# start time
    t_ed = 50.0      	# end time
    # -- Physical Params --
    Re_inv = 1.0 / 8.0
    κ = 0.0 # sqrt(30)
    U_y = 0.0
    n_force = 4

    c = configure_FFT2D(nx, ny,
		        xmin, xmax, ymin, ymax)
    s = set_HirutaFlow(nt,
		       t_st, t_ed,
		       Re_inv, κ, U_y, n_force)    

    yx_X = yx_Func(c.yx_X, c)
    yx_Y = yx_Func(c.yx_Y, c)
    println("s.n_force = ", s.n_force)
    yx_ncosny = s.n_force * cos(s.n_force * yx_Y)
    yx_ncosy = s.n_force * cos(yx_Y)
    open("test_ncosny.dat", "w") do f
        for ix = 1:c.nx
	    for iy = 1:c.ny
	        println(f, yx_X.vals[iy, ix], " ",
			   yx_Y.vals[iy, ix], " ",
			   yx_ncosny.vals[iy, ix], " ",
			   yx_ncosy.vals[iy, ix])
	    end
	    println(f)
	end
	println(f)
    end

end