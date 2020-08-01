using FFTW

struct FFT2D_config

    # *** grid points ***
    nx::Int64
    ny::Int64
    # *** range ***
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    # *** corrdinate ***
    yx_X::Matrix{Float64}
    yx_Y::Matrix{Float64}
    lk_K::Matrix{Float64}
    lk_L::Matrix{Float64}
    # *** dealiasing flag ***
    dealiasing::Bool
    
end


function configure_FFT2D(nx::Int64, ny::Int64,
                         xmin::Float64, xmax::Float64,
		         ymin::Float64, ymax::Float64;
			 dealiasing=true
			 )::FFT2D_config

    # *** initialize x-coordinate ***
    len_x = xmax - xmin
    yx_X = Matrix{Float64}(undef, ny, nx)
    yx_X[1, :] = [xmin + len_x*ix/nx for ix = 0:nx-1]
    for iy = 2:ny
        yx_X[iy, :] = yx_X[1, :]
    end

    # *** initialize y-coordinate ***
    len_y = ymax - ymin
    yx_Y = Matrix{Float64}(undef, ny, nx)
    yx_Y[:, 1] = [ymin + len_y*iy/ny for iy = 0:ny-1]
    for ix = 2:nx
        yx_Y[:, ix] = yx_Y[:, 1]
    end

    # *** initialize k-coordinate ***
    K = [2pi/len_x * ix
         for ix = 0:div(nx,2)-1]
    if nx % 2 != 0
        push!(K, 0.0)
    end
    append!(K,
	    [2pi/len_x * (ix - div(nx, 2))
	     for ix = 0:div(nx,2)-1])

    lk_K = Matrix{Complex{Float64}}(undef, ny, nx)
    for il = 1:ny
        lk_K[il, :] = K
    end

    # *** initialize l-coordinate ***
    L = [2pi/len_y * iy
         for iy = 0:div(ny,2)-1]
    if ny % 2 != 0
        push!(L, 0.0)
    end
    append!(L,
	    [2pi/len_y * (iy - div(ny, 2))
	     for iy = 0:div(ny,2)-1])

    lk_L = Matrix{Complex{Float64}}(undef, ny, nx)
    for ik = 1:nx
        lk_L[:, ik] = L
    end

    return FFT2D_config(nx, ny,
    	   		xmin, xmax, ymin, ymax,
			yx_X, yx_Y, lk_K, lk_L,
			dealiasing)

end


# ***** yx_Func *************************
mutable struct yx_Func

    vals::Matrix{Float64}
    config::FFT2D_config

end

Base.:+(f::yx_Func) = f
Base.:+(f::yx_Func,
	g::yx_Func) = yx_Func(f.vals + g.vals,
		      	      f.config)
Base.:+(f::yx_Func,
	g::Float64) = yx_Func(f.vals + g,
		      	      f.config)
Base.:+(f::Float64,
	g::yx_Func) = yx_Func(f + g.vals,
		      	      g.config)			      

Base.:-(f::yx_Func) = yx_Func(- f.vals, f.config)
Base.:-(f::yx_Func,
	g::yx_Func) = yx_Func(f.vals - g.vals,
		      	      f.config)
Base.:-(f::yx_Func,
	g::Float64) = yx_Func(f.vals - g,
		      	      f.config)
Base.:-(f::Float64,
	g::yx_Func) = yx_Func(f - g.vals,
		      	      g.config)

Base.:*(f::yx_Func,
	g::yx_Func) = yx_Func(f.vals .* g.vals,
		      	      f.config)
Base.:*(f::yx_Func,
	g::Float64) = yx_Func(f.vals * g,
		      	      f.config)
Base.:*(f::Float64,
	g::yx_Func) = yx_Func(f * g.vals,
		              g.config)

Base.:/(f::yx_Func,
	g::yx_Func) = yx_Func(f.vals ./ g.vals,
		      	      f.config)
Base.:/(f::yx_Func,
	g::Float64) = yx_Func(f.vals / g,
		      	      f.config)

Base.:sin(f::yx_Func) = yx_Func(sin.(f.vals), f.config)
Base.:cos(f::yx_Func) = yx_Func(cos.(f.vals), f.config)
Base.:exp(f::yx_Func) = yx_Func(exp.(f.vals), f.config)


# ***** lk_Func *************************
mutable struct lk_Func

    vals::Matrix{Complex{Float64}}
    config::FFT2D_config

end

Base.:+(f::lk_Func,
	g::lk_Func) = lk_Func(f.vals + g.vals,
		      	      f.config)
Base.:+(f::lk_Func,
	g::Complex{Float64}) = lk_Func(f.vals + g,
		      	               f.config)
Base.:+(f::Complex{Float64},
	g::lk_Func) = lk_Func(f + g.vals,
		      	      g.config)			      

Base.:*(f::lk_Func,
	g::lk_Func) = lk_Func(f.vals .* g.vals,
		      	      f.config)
Base.:*(f::lk_Func,
	g::Complex{Float64}) = lk_Func(f.vals * g,
		      	               f.config)
Base.:*(f::Complex{Float64},
	g::lk_Func) = lk_Func(f * g.vals,
		              g.config)

Base.:/(f::lk_Func,
	g::lk_Func) = lk_Func(f.vals ./ g.vals,
		      	      f.config)
Base.:/(f::lk_Func,
	g::Complex{Float64}) = lk_Func(f.vals / g,
		      	               f.config)


# ***** Fourier Transformation **********
function lk_yx(yx_func::yx_Func)::lk_Func

    return lk_Func(fft(complex(yx_func.vals)),
		   yx_func.config)

end

function yx_lk(lk_func::lk_Func)::yx_Func

    lk_vals = copy(lk_func.vals)
    if lk_func.config.dealiasing
        len_l, len_k = size(lk_vals)
	lk_vals[div(len_l, 3):div(2*len_l + 2, 3), :] .= 0.0
	lk_vals[:, div(len_k, 3):div(2*len_k + 2, 3)] .= 0.0
    end
    return yx_Func(real(ifft(lk_vals)), lk_func.config)

end

# **** differenciation by spectral method ****
function lk_Dx_lk(lk_func::lk_Func)::lk_Func

    c = lk_func.config
    lk_vals = im * c.lk_K .* lk_func.vals
    return lk_Func(lk_vals, c)

end


function lk_Dy_lk(lk_func::lk_Func)::lk_Func

    c = lk_func.config
    lk_vals = im * c.lk_L .* lk_func.vals
    return lk_Func(lk_vals, c)

end


# ***** test program ************************
function test_FFT2D_derivative(;dealiasing=true)

    # *** parameter setting ***
    nx = 128	    # # of grid points
    ny = 128	    # # of grid points
    xmin = 0.0	    # leftmost(lowest) value of x
    xmax = 4pi	    # rightmost(highest) value of x
    ymin = 0.0	    # leftmost(lowest) value of y
    ymax = 2pi	    # rightmost(highest) value of y
    c = configure_FFT2D(nx, ny,
		        xmin, xmax, ymin, ymax,
			dealiasing=dealiasing)

    yx_X = yx_Func(c.yx_X, c)   # x coordinate
    yx_Y = yx_Func(c.yx_Y, c)   # y coordinate

    # *** differentiate by spectral method ***
    yx_func = exp(sin(yx_X)) * cos(yx_Y)

    lk_func = lk_yx(yx_func)
    lk_dfunc = lk_Dy_lk(lk_Dx_lk(lk_func))
    yx_dfunc = yx_lk(lk_dfunc)

    # *** exact derivative ***
    yx_exact = cos(yx_X) * exp(sin(yx_X)) * (- sin(yx_Y))

    # *** output ***
    open("testcase.txt", "w") do f
        for ix = 1:nx
	    for iy = 1:ny
	        println(f, yx_X.vals[iy, ix], " ",
			   yx_Y.vals[iy, ix], " ",
			   yx_dfunc.vals[iy, ix], " ",
			   yx_exact.vals[iy, ix])
	    end
	    println(f)
	end
    end

end