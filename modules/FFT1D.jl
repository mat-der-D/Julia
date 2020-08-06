import Base.size, Base.getindex, Base.setindex!
import Base.==, Base.copy
import Base.+, Base.-, Base.*, Base./
import Base.sin, Base.cos, Base.tan, Base.cot, Base.sec, Base.csc
import Base.asin, Base.acos, Base.atan, Base.acot, Base.asec, Base.acsc
import Base.sinh, Base.cosh, Base.tanh, Base.coth, Base.sech, Base.csch
import Base.asinh, Base.acosh, Base.atanh, Base.acoth, Base.asech, Base.acsch
import Base.exp, Base.log

using FFTW


# **************
# *** config ***
# **************
struct FFT1D_config

    nx::Int64	# num of grid Pts.
    xmin::Float64     # minimum val. of x
    xmax::Float64     # maximum val. of x
    x_X::Array{Float64, 1}	# x-coordinate
    k_K::Array{Complex{Float64}, 1}	# k-coordinate
    dealiasing::Bool

end


function ==(c1::FFT1D_config, c2::FFT1D_config)
    return ( c1.nx         == c2.nx   &&
    	     c1.xmin       == c2.xmin &&
	     c1.xmax       == c2.xmax &&
	     c1.dealiasing == c2.dealiasing )
end


function copy(c::FFT1D_config)
    return FFT1D_config(c.nx, c.xmin, c.xmax,
    	   		c.x_X, c.k_K,
			c.dealiasing)
end


function configure_FFT1D(; nx::Int64,
			   xmin::Float64, xmax::Float64,
			   dealiasing::Bool=true
			 )::FFT1D_config

    # *** initialize x-coordinate ***
    dx = (xmax - xmin) / nx
    x_X = [xmin + n*dx for n = 0:nx-1]

    # *** initialize k-coordinate ***
    hx = 2pi/(xmax - xmin)
    k_K = complex(zeros(nx))
    Nx = div(nx-1, 2)
    for n = 1 : Nx
    	k_K[begin + n] = n * hx
    end
    for n = Nx : -1 : 1
        k_K[end - n + 1] = - n * hx
    end

    return FFT1D_config(nx, xmin, xmax, x_X, k_K, dealiasing)

end


# **********************
# *** x_Func, k_Func ***
# **********************
# ----- x_Func -----
mutable struct x_Func <: AbstractArray{Float64, 1}

    vals::Array{Float64, 1}
    config::FFT1D_config

end

size(f::x_Func) = size(f.vals)
getindex(f::x_Func, i::Int) = getindex(f.vals, i)
setindex!(f::x_Func, v, i::Int) = setindex!(f.vals, v, i)

function x_Func_undef(c::FFT1D_config)
    return x_Func(Array{Float64, 1}(undef, c.nx), c)
end


function ==(f::x_Func, g::x_Func)
    return f.val == g.val && f.config == g.config
end


function copy(f::x_Func)
    return x_Func(f.vals, f.config)
end

# ----- OPERATORS -----
# + operator
+(f::x_Func) = f

function +(f::x_Func, g::x_Func)
    if f.config != g.config
        error("FFT1D_config of two x_Func are not equal.")
    else
        return x_Func(f.vals + g.vals, f.config)
    end
end

+(f::x_Func, a::Real) = x_Func(f.vals .+ a, f.config)
+(a::Real, f::x_Func) = f + a

# - operator
-(f::x_Func) = x_Func(- f.vals, f.config)

function -(f::x_Func, g::x_Func)
    if f.config != g.config
        error("FFT1D_config of two x_Func are not equal.")
    else
        return x_Func(f.vals - g.vals, f.config)
    end
end

-(f::x_Func, a::Real) = x_Func(f.vals .- a, f.config)
-(a::Real, f::x_Func) = x_Func(a .- f.vals, f.config)

# * operator
function *(f::x_Func, g::x_Func)
    if f.config != g.config
        error("FFT1D_config of two x_Func are not equal.")
    else
        return x_Func(f.vals .* g.vals, f.config)
    end
end

*(f::x_Func, a::Real) = x_Func(f.vals * a, f.config)
*(a::Real, f::x_Func) = f * a


# / operator
function /(f::x_Func, g::x_Func)
    if f.config != g.config
        error("FFT1D_config of two x_Func are not equal.")
    else
        return x_Func(f.vals ./ g.vals, f.config)
    end
end

/(f::x_Func, a::Real) = x_Func(f.vals / a, f.config)
/(a::Real, f::x_Func) = x_Func(a ./ f.vals, f.config)


# ----- k_Func -----
mutable struct k_Func <: AbstractArray{Complex{Float64}, 1}

    vals::Array{Complex{Float64}, 1}
    config::FFT1D_config

end

size(f::k_Func) = size(f.vals)
getindex(f::k_Func, i::Int) = getindex(f.vals, i)
setindex!(f::k_Func, v, i::Int) = setindex!(f.vals, v, i)


function k_Func_undef(c::FFT1D_config)
    return k_Func(Array{Complex{Float64}, 1}(undef, c.nx), c)
end


function ==(f::k_Func, g::k_Func)
    return f.val == g.val && f.config == g.config
end


function copy(f::k_Func)
    return k_Func(f.vals, f.config)
end

# ----- OPERATORS -----
# + operator
+(f::k_Func) = f

function +(f::k_Func, g::k_Func)
    if f.config != g.config
        error("FFT1D_config of two k_Func are not equal.")
    else
        return k_Func(f.vals + g.vals, f.config)
    end
end

+(f::k_Func, a::Number) = k_Func(k.vals .+ a, f.config)
+(a::Number, f::k_Func) = f + a

# - operator
-(f::k_Func) = k_Func(- f.vals, f.config)

function -(f::k_Func, g::k_Func)
    if f.config != g.config
        error("FFT1D_config of two k_Func are not equal.")
    else
        return k_Func(f.vals - g.vals, f.config)
    end
end

-(f::k_Func, a::Number) = k_Func(f.vals .- a, f.config)
-(a::Number, f::k_Func) = k_Func(a .- f.vals, f.config)

# * operator
function *(f::k_Func, g::k_Func)
    if f.config != g.config
        error("FFT1D_config of two k_Func are not equal.")
    else
        return k_Func(f.vals .* g.vals, f.config)
    end
end

*(f::k_Func, a::Number) = k_Func(f.vals * a, f.config)
*(a::Number, f::k_Func) = k_Func(a * f.vals, f.config)

# / operator
function /(f::k_Func, g::k_Func)
    if f.config != g.config
        error("FFT1D_config of two k_Func are not equal.")
    else
        return k_Func(f.vals ./ g.vals, f.config)
    end
end

/(f::k_Func, a::Number) = k_Func(f.vals / a, f.config)
/(a::Number, f::k_Func) = k_Func(a ./ f.vals, f.config)


# ----- FOURIER TRANSFORMATION -----
function k_x(f::x_Func)::k_Func
    return k_Func(fft(f.vals), f.config)
end

function x_k(f::k_Func)::x_Func

    if f.config.dealiasing
        f = copy(f)
	L = f.config.nx
	f.vals[div(L, 3) : div(2*L + 2, 3)] .= 0.0
    end
    return x_Func(real(ifft(f.vals)), f.config)
end


# ******************
# *** Derivative ***
# ******************
function k_Dx_k(k_func::k_Func)::k_Func
    c = k_func.config
    k_K = k_Func(c.k_K, c)
    return im * k_K * k_func
end


# ****************************
# *** Elementary Functions ***
# ****************************
# For x_Func
sin(f::x_Func) = x_Func(sin.(f.vals), f.config)
cos(f::x_Func) = x_Func(cos.(f.vals), f.config)
tan(f::x_Func) = x_Func(tan.(f.vals), f.config)
cot(f::x_Func) = x_Func(cot.(f.vals), f.config)
sec(f::x_Func) = x_Func(sec.(f.vals), f.config)
csc(f::x_Func) = x_Func(csc.(f.vals), f.config)

asin(f::x_Func) = x_Func(asin.(f.vals), f.config)
acos(f::x_Func) = x_Func(acos.(f.vals), f.config)
atan(f::x_Func) = x_Func(atan.(f.vals), f.config)
acot(f::x_Func) = x_Func(acot.(f.vals), f.config)
asec(f::x_Func) = x_Func(asec.(f.vals), f.config)
acsc(f::x_Func) = x_Func(acsc.(f.vals), f.config)

sinh(f::x_Func) = x_Func(sinh.(f.vals), f.config)
cosh(f::x_Func) = x_Func(cosh.(f.vals), f.config)
tanh(f::x_Func) = x_Func(tanh.(f.vals), f.config)
coth(f::x_Func) = x_Func(coth.(f.vals), f.config)
sech(f::x_Func) = x_Func(sech.(f.vals), f.config)
csch(f::x_Func) = x_Func(csch.(f.vals), f.config)

asinh(f::x_Func) = x_Func(asinh.(f.vals), f.config)
acosh(f::x_Func) = x_Func(acosh.(f.vals), f.config)
atanh(f::x_Func) = x_Func(atanh.(f.vals), f.config)
acoth(f::x_Func) = x_Func(acoth.(f.vals), f.config)
asech(f::x_Func) = x_Func(asech.(f.vals), f.config)
acsch(f::x_Func) = x_Func(acsch.(f.vals), f.config)

exp(f::x_Func) = x_Func(exp.(f.vals), f.config)
log(f::x_Func) = x_Func(log.(f.vals), f.config)

# For k_Func
sin(f::k_Func) = k_Func(sin.(f.vals), f.config)
cos(f::k_Func) = k_Func(cos.(f.vals), f.config)
tan(f::k_Func) = k_Func(tan.(f.vals), f.config)
cot(f::k_Func) = k_Func(cot.(f.vals), f.config)
sec(f::k_Func) = k_Func(sec.(f.vals), f.config)
csc(f::k_Func) = k_Func(csc.(f.vals), f.config)

asin(f::k_Func) = k_Func(asin.(f.vals), f.config)
acos(f::k_Func) = k_Func(acos.(f.vals), f.config)
atan(f::k_Func) = k_Func(atan.(f.vals), f.config)
acot(f::k_Func) = k_Func(acot.(f.vals), f.config)
asec(f::k_Func) = k_Func(asec.(f.vals), f.config)
acsc(f::k_Func) = k_Func(acsc.(f.vals), f.config)

sinh(f::k_Func) = k_Func(sinh.(f.vals), f.config)
cosh(f::k_Func) = k_Func(cosh.(f.vals), f.config)
tanh(f::k_Func) = k_Func(tanh.(f.vals), f.config)
coth(f::k_Func) = k_Func(coth.(f.vals), f.config)
sech(f::k_Func) = k_Func(sech.(f.vals), f.config)
csch(f::k_Func) = k_Func(csch.(f.vals), f.config)

asinh(f::k_Func) = k_Func(asinh.(f.vals), f.config)
acosh(f::k_Func) = k_Func(acosh.(f.vals), f.config)
atanh(f::k_Func) = k_Func(atanh.(f.vals), f.config)
acoth(f::k_Func) = k_Func(acoth.(f.vals), f.config)
asech(f::k_Func) = k_Func(asech.(f.vals), f.config)
acsch(f::k_Func) = k_Func(acsch.(f.vals), f.config)

exp(f::k_Func) = k_Func(exp.(f.vals), f.config)
log(f::k_Func) = k_Func(log.(f.vals), f.config)


# ********************
# *** test program ***
# ********************
function test_FFT1D(;dealiasing = true)

    # *** parameter setting ***
    nx = 50	    # # of grid points
    xmin = 0.0	    # leftmost(lowest) value of x
    xmax = 2pi	    # rightmost(highest) value of x

    c = configure_FFT1D(nx=nx, xmin=xmin, xmax=xmax,
      			dealiasing=dealiasing)
    x_X = x_Func(c.x_X, c) # x-coordinate

    # *** differentiate by spectral method ***
    x_func = exp(sin(x_X))
    k_func = k_x(x_func)
    k_dfunc = k_Dx_k(k_func)
    x_dfunc = x_k(k_dfunc)

    # *** exact derivative ***
    x_exact = cos(x_X) * x_func

    # *** output ***
    open("testcase.txt", "w") do f
        for ix = 1:nx
            println(f,
		    x_X.vals[ix], " ",
		    x_func.vals[ix], " ",
		    x_dfunc.vals[ix], " ",
		    x_exact.vals[ix])
	    # OUTPUT:
	    #   x, f(x), f'(x)(numerical), f'(x)(exact)
        end
    end

end