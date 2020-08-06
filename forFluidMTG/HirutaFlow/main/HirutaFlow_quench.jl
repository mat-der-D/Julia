include("../../../modules/FFT2D.jl")


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


function set_HirutaFlow(; nt::Int64,
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
			dtFactor::Float64=1.0)

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


function lk_TimeDevelopRK4_lk(
			lk_ω::lk_Func,
			s::HirutaFlow_settings,
			tit::TimeIntegrationTools)::lk_Func

    # Note: initialize tit as dtFactor=0.5

    dt = s.dt
    lk_ω0 = copy(lk_ω)

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
    lk_ω3 = lk_linear_next_lk(
    	        lk_linear_next_lk(lk_ω0, tit)
		+ lk_k3, tit)
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
    lk_ψ = copy(lk_ω)

    # avoid NaN
    lk_KKLL[1, 1] = 1.
    lk_ψ[1, 1] = 0.

    if c.nx % 2 == 0 && c.ny % 2 == 0
    
        half_nx = div(c.nx, 2)
	half_ny = div(c.ny, 2)
	
        lk_KKLL[1, 1 + half_nx] = 1.
	lk_KKLL[1 + half_ny, 1] = 1.
	lk_KKLL[1 + half_ny, 1 + half_nx] = 1.
	
	lk_ψ[1, 1 + half_nx] = 0.
	lk_ψ[1 + half_ny, 1] = 0.
	lk_ψ[1 + half_ny, 1 + half_nx] = 0.
	
    elseif c.nx % 2 == 0 # but not c.ny % 2 == 0

        half_nx = div(c.nx, 2)
	lk_KKLL[1, 1 + half_nx] = 1.
	lk_ψ[1, 1 + half_nx] = 0.

    elseif c.ny % 2 == 0 # but not c.nx % 2 == 0

        half_ny = div(c.ny, 2)
	lk_KKLL[1 + half_ny, 1] = 1.
	lk_ψ[1 + half_ny, 1] = 0.

    end

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
# ***** main program *****
# ************************
function main_HirutaFlow()


    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # @@@ CONFIGURATION
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    # *** configure FFT2D *****
    # -- Space --
    aspect = 4		 # aspect ratio
    nx = 64*aspect    	 # # of grid points
    ny = 64	      	 # # of grid points
    xmin = 0.0	      	 # minimum value of x
    xmax = 2pi*aspect    # maximum value of x
    ymin = 0.0	    	 # minimum value of y
    ymax = 2pi	    	 # maximum value of y

    c = configure_FFT2D(
		nx=nx, ny=ny,
		xmin=xmin, xmax=xmax,
		ymin=ymin, ymax=ymax)

    # *** PARAMETER CASE 1 ***
    # -- Time --
    nt = 5000		# # of time steps
    t_st = 0.	    	# start time
    t_ed = 10.      	# end time
    # -- Physical Params --
    Re_inv = 1.0 / 30.0
    κ = 0.0 # sqrt(30)
    U_y = 0.6
    n_force = 4

    s1 = set_HirutaFlow(nt=nt,
		        t_st=t_st, t_ed=t_ed,
		        Re_inv=Re_inv,
			κ=κ, U_y=U_y, n_force=n_force)
			
    tit1 = set_TimeIntegrationTools(s1, c, dtFactor=0.5)


    # *** PARAMETER CASE 2 ***
    # -- Time --
    nt = 5000		# # of time steps
    t_st = 10.		# start time
    t_ed = 20.      	# end time
    # -- Physical Params --
    Re_inv = 1.0 / 12.0
    κ = κ
    U_y = U_y
    n_force = n_force

    s2 = set_HirutaFlow(nt=nt,
		        t_st=t_st, t_ed=t_ed,
		        Re_inv=Re_inv,
			κ=κ, U_y=U_y, n_force=n_force)
			
    tit2 = set_TimeIntegrationTools(s1, c, dtFactor=0.5)


    # *** output settings ***
    nt_alert  = 50	# println every nt_alert
    nt_output = 50 	# output every nt_output
    nt_st_out = 0	# start ouput at it = nt_st_out


    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # @@@ TIME INTEGRATION
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    open("vorticity.dat", "w") do f
        lk_ω = lk_Func_undef(c)
	set_initial_condition(lk_ω)
	println("t = ", s1.t_st)

	# MEMO: Initial Condition is not output

	# *** PARAMETER CASE 1 PHASE ***
	for it = 1 : s1.nt
	    lk_ω = lk_TimeDevelopRK4_lk(lk_ω, s1, tit1)
	    # --- stdout ---
	    if it % nt_alert == 0
	        t_now = s1.t_st + it * s1.dt
	        println("t = ", t_now)
	    end
	end

	# --- file output ---
	if ( 0 % nt_output == 0 &&
	     0 >= nt_st_out )
	    output_flow_field(f, lk_ω)
	    println("OUTPUT")
        end

	# *** PARAMETER CASE 2 PHASE ***
	for it = 1 : s2.nt
	    lk_ω = lk_TimeDevelopRK4_lk(lk_ω, s2, tit2)

	    # --- stdout ---
	    if it % nt_alert == 0
	        t_now = s2.t_st + it * s2.dt
		println("t = ", t_now)
	    end
	    # --- file output ---
	    if ( it % nt_output == 0 &&
	         it >= nt_st_out )
		 output_flow_field(f, lk_ω)
		 println("OUTPUT")
            end
	end

    end


end


function output_flow_field(f::IOStream, lk_ω::lk_Func)

    c = lk_ω.config
    yx_ω = yx_lk(lk_ω)
    for ix = 1:c.nx
        for iy = 1:c.ny
	    println(f, c.yx_X[iy, ix], " ",
	    	       c.yx_Y[iy, ix], " ",
		       yx_ω[iy, ix])
	end
	println(f)
    end
    println(f)

end


function set_initial_condition(lk_ω::lk_Func)
    c = lk_ω.config
    lk_ω.vals = rand(c.ny, c.nx)
end