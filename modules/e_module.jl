using FFTW


# *******************************************
#  Configuration
# *******************************************
struct ConfigFFT1D

    nx::Int        # Number of Grids
    nk::Int        # Cut-off wavenumber, nk <= nx
    xmin::Float64  # Min value of x
    xmax::Float64  # Max value of x
    x_X::Array{Float64, 1} # x-corrdinate

end


function configure_FFT1D(; nx::Int, nk::Int,
	 		   xmin::Float64, xmax::Float64
			)::ConfigFFT1D
			
    dx = (xmax - xmin) / nx
    x_X = [xmin + ix*dx for ix = 1:nx-1]
    return ConfigFFT1D(nx, nk, xmin, xmax, x_X)

end


# *******************************************
#  x_Func, e_Func
# *******************************************
# +++++ x_Func +++++
mutable struct x_Func <: AbstractArray{Float64, 1}

    vals::Array{Float64, 1}
    config::ConfigFFT1D
    
    # constructor
    function x_Func(vals::Array{Float64, 1}, config::ConfigFFT1D)
        if size(vals) == (config.nx,)
	    return new(copy(vals), config)
	else
	    println("ERROR")
	end
    end

    function x_Func(vals::Array{Int, 1}, config::ConfigFFT1D)
        return x_Func(float(vals), config)
    end

end


Base.:size(f::x_Func) = size(f.vals)
Base.:getindex(f::x_Func, i::Int) = getindex(f.vals, i)
Base.:setindex!(f::x_Func, v, i::Int) = setindex!(f.vals, v, i)

Base.:copy(f::x_Func) = x_Func(f.vals, f.config)

# --- operators ---

# To Be Implemented!



# +++++ k_Func +++++
mutable struct k_Func <: AbstractArray{Float64, 1}

    vals::Array{Complex{Float64}, 1}
    config::ConfigFFT1D
    
    # constructor
    function k_Func(vals::Array{Complex{Float64}, 1},
    	     	    config::ConfigFFT1D)
        if size(vals) == (div(config.nx, 2) + 1,)
	    return new(copy(vals), config)
	else
	    println("ERROR")
	end
    end

    function k_Func(vals::Array{Real, 1}, config::ConfigFFT1D)
        return x_Func(float(complex(vals)), config)
    end

end


Base.:size(f::k_Func) = size(f.vals)
Base.:getindex(f::k_Func, i::Int) = getindex(f.vals, i)
Base.:setindex!(f::k_Func, v, i::Int) = setindex!(f.vals, v, i)

Base.:copy(f::k_Func) = k_Func(f.vals, f.config)