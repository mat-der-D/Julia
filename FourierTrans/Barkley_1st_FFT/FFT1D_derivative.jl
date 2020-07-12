using FFTW

struct FFT1D_config

    nx::Int64			# # of grid points
    xl::Float64			# leftmost(lowest) value of x
    xr::Float64			# rightmost(highest) value of x
    x_X::Vector{Float64}	# x coordinates
    k_X::Vector{Float64}	# k corrdinates

end

function init_FFT1D(nx::Int64,
                    xl::Float64, xr::Float64)::FFT1D_config

    # *** initialize x-coordinate ***
    len_x = xr - xl
    x_X = [xl + len_x*n/nx for n = 0:nx-1]

    # *** initialize k-coordinate ***
    hx = 2pi/len_x
    k_X = [hx*n for n = 0:div(nx,2)-1]
    if nx % 2 != 0
        push!(k_X, 0.0)
    end    
    append!(k_X,
	  [hx*(n - div(nx,2)) for n = 0:div(nx,2)-1])
    return FFT1D_config(nx, xl, xr, x_X, k_X)

end

function k_x(x_func::Vector{Float64})::Vector{Complex{Float64}}
    return fft(complex(x_func))
end

function x_k(k_func::Vector{Complex{Float64}};
	     dealiasing = true)::Vector{Float64}
    if dealiasing
        len_k = length(k_func)
	k_func[div(len_k, 3):div(2*len_k + 2, 3)] .= 0.0
    end
    return real(ifft(k_func))
end

function k_DX_k(k_func::Vector{Complex{Float64}},
		config::FFT1D_config)::Vector{Complex{Float64}}
    return (im * config.k_X) .* k_func
end

function x_DX_x(x_func::Vector{Float64},
		config::FFT1D_config;
		dealiasing = true)::Vector{Float64}
    return x_k(k_DX_k(k_x(x_func), config), dealiasing = dealiasing)
end

function test_FFT1D_derivative(;dealiasing = true)

    # *** parameter setting ***
    nx = 50	    # # of grid points
    xl = 0.0	    # leftmost(lowest) value of x
    xr = 2pi	    # rightmost(highest) value of x
    c = init_FFT1D(nx, xl, xr)
    x_X = c.x_X	    # x coordinate

    # *** differentiate by spectral method ***
    x_func = exp.(sin.(x_X))
    k_func = k_x(x_func)
    k_dfunc = k_DX_k(k_func, c)
    x_dfunc = x_k(k_dfunc, dealiasing = dealiasing)

    # *** exact derivative ***
    x_exact = cos.(x_X) .* x_func

    # *** output ***
    open("testcase.txt", "w") do f
        for i = 1:nx
            println(f, x_X[i], " ", x_dfunc[i], " ", x_exact[i])
        end
    end

end