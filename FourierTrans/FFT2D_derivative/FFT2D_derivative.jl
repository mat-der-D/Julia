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

end


function init_FFT2D(nx::Int64, ny::Int64,
                    xmin::Float64, xmax::Float64,
		    ymin::Float64, ymax::Float64)::FFT2D_config

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
			yx_X, yx_Y, lk_K, lk_L)

end


function lk_yx(yx_func::Matrix{Float64})::Matrix{Complex{Float64}}
    return fft(complex(yx_func))
end

function yx_lk(lk_func::Matrix{Complex{Float64}};
	       dealiasing=true)::Matrix{Float64}

    if dealiasing
        len_l, len_k = size(lk_func)
	lk_func[div(len_l, 3):div(2*len_l + 2, 3), :] .= 0.0
	lk_func[:, div(len_k, 3):div(2*len_k + 2, 3)] .= 0.0
    end
    return real(ifft(lk_func))

end


function lk_Dx_lk(lk_func::Matrix{Complex{Float64}},
		  config::FFT2D_config)::Matrix{Complex{Float64}}
		  
    return im * config.lk_K .* lk_func

end


function lk_Dy_lk(lk_func::Matrix{Complex{Float64}},
		  config::FFT2D_config)::Matrix{Complex{Float64}}

    return im * config.lk_L .* lk_func

end


function test_FFT2D_derivative(;dealiasing = true)

    # *** parameter setting ***
    nx = 128	    # # of grid points
    ny = 128	    # # of grid points
    xmin = 0.0	    # leftmost(lowest) value of x
    xmax = 4pi	    # rightmost(highest) value of x
    ymin = 0.0	    # leftmost(lowest) value of y
    ymax = 2pi	    # rightmost(highest) value of y
    c = init_FFT2D(nx, ny, xmin, xmax, ymin, ymax)
    yx_X = c.yx_X   # x coordinate
    yx_Y = c.yx_Y   # y coordinate

    # *** differentiate by spectral method ***
    yx_func = exp.(sin.(yx_X)) .* cos.(yx_Y)
    lk_func = lk_yx(yx_func)
    lk_dfunc = lk_Dy_lk(lk_Dx_lk(lk_func, c), c)
    yx_dfunc = yx_lk(lk_dfunc, dealiasing=dealiasing)

    # *** exact derivative ***
    yx_exact = cos.(yx_X) .* exp.(sin.(yx_X)) .* (- sin.(yx_Y))

    # *** output ***
    open("testcase.txt", "w") do f
        for ix = 1:nx
	    for iy = 1:ny
	        println(f, yx_X[iy, ix], " ",
			   yx_Y[iy, ix], " ",
			   yx_dfunc[iy, ix], " ",
			   yx_exact[iy, ix])
	    end
	    println(f)
	end
    end

end