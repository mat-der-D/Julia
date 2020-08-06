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
struct FFT2D_config

    nx::Int64	# num of grid Pts. (x)
    ny::Int64	# num of grid Pts. (y)
    xmin::Float64     # minimum val. of x
    xmax::Float64     # maximum val. of x
    ymin::Float64     # minimum val. of y
    ymax::Float64     # maximum val. of y
    yx_X::Array{Float64, 2}	# x-coordinate
    yx_Y::Array{Float64, 2}	# y-coordinate
    lk_K::Array{Complex{Float64}, 2}	# k-coordinate
    lk_L::Array{Complex{Float64}, 2}	# l-coordinate
    dealiasing::Bool

end


function ==(c1::FFT2D_config, c2::FFT2D_config)
    return ( c1.nx         == c2.nx   &&
    	     c1.ny	   == c2.ny   &&
    	     c1.xmin       == c2.xmin &&
	     c1.xmax       == c2.xmax &&
	     c1.ymin	   == c2.ymin &&
	     c1.ymax	   == c2.ymax &&
	     c1.dealiasing == c2.dealiasing )
end


function copy(c::FFT2D_config)
    return FFT2D_config(c.nx, c.ny,
    	   		c.xmin, c.xmax, c.ymin, c.ymax,
			c.yx_X, c.yx_Y, c.lk_K, c.lk_K,
			c.config)
end


function configure_FFT2D(; nx::Int64, ny::Int64,
			   xmin::Float64, xmax::Float64,
			   ymin::Float64, ymax::Float64,
			   dealiasing::Bool=true
			 )::FFT2D_config

    # *** initialize x-coordinate ***
    dx = (xmax - xmin) / nx
    yx_X = zeros(ny, nx)
    yx_X .= [xmin + n*dx for n = 0:nx-1]

    # *** initialize x-coordinate ***
    dy = (ymax - ymin) / ny
    yx_Y = zeros(ny, nx)
    yx_Y .= [ymin + n*dy for n = 0:ny-1]'

    # *** initialize k-coordinate ***
    hx = 2pi/(xmax - xmin)
    lk_K = complex(zeros(ny, nx))
    Nx = div(nx-1, 2)
    for n = 1 : Nx
    	lk_K[begin + n, :] .= n * hx
    end
    for n = Nx : -1 : 1
        lk_K[end - n + 1, :] .= - n * hx
    end

    # *** initialize l-coordinate ***
    hy = 2pi/(ymax - ymin)
    lk_L = complex(zeros(ny, nx))
    Ny = div(ny-1, 2)
    for n = 1 : Ny
    	lk_L[:, begin + n] .= n * hx
    end
    for n = Ny : -1 : 1
        lk_L[:, end - n + 1] .= - n * hx
    end

    return FFT2D_config(nx, ny,
    	   		xmin, xmax, ymin, ymax,
			yx_X, yx_Y, lk_K, lk_L,
			dealiasing)

end


# ************************
# *** yx_Func, lk_Func ***
# ************************
# ----- yx_Func -----
mutable struct yx_Func

    vals::Array{Float64, 2}
    config::FFT2D_config

end


function yx_Func_undef(c::FFT2D_config)
    return yx_Func(Array{Float64, 2}(undef, c.ny, c.nx), c)
end


function ==(f::yx_Func, g::yx_Func)
    return f.val == g.val && f.config == g.config
end


function copy(f::yx_Func)
    return yx_Func(f.vals, f.config)
end

# ----- OPERATORS -----
# + operator
+(f::yx_Func) = f

function +(f::yx_Func, g::yx_Func)
    if f.config != g.config
        error("FFT2D_config of two yx_Func are not equal.")
    else
        return yx_Func(f.vals + g.vals, f.config)
    end
end

+(f::yx_Func, a::Real) = yx_Func(f.vals .+ a, f.config)
+(a::Real, f::yx_Func) = f + a

# - operator
-(f::yx_Func) = yx_Func(- f.vals, f.config)

function -(f::yx_Func, g::yx_Func)
    if f.config != g.config
        error("FFT2D_config of two yx_Func are not equal.")
    else
        return yx_Func(f.vals - g.vals, f.config)
    end
end

-(f::yx_Func, a::Real) = yx_Func(f.vals .- a, f.config)
-(a::Real, f::yx_Func) = yx_Func(a .- f.vals, f.config)

# * operator
function *(f::yx_Func, g::yx_Func)
    if f.config != g.config
        error("FFT2D_config of two_Func are not equal.")
    else
        return yx_Func(f.vals .* g.vals, f.config)
    end
end

*(f::yx_Func, a::Real) = yx_Func(f.vals * a, f.config)
*(a::Real, f::yx_Func) = f * a


# / operator
function /(f::yx_Func, g::yx_Func)
    if f.config != g.config
        error("FFT2D_config of two yx_Func are not equal.")
    else
        return yx_Func(f.vals ./ g.vals, f.config)
    end
end

/(f::yx_Func, a::Real) = yx_Func(f.vals / a, f.config)
/(a::Real, f::yx_Func) = yx_Func(a ./ f.vals, f.config)


# ----- lk_Func -----
mutable struct lk_Func

    vals::Array{Complex{Float64}, 2}
    config::FFT2D_config

end


function lk_Func_undef(c::FFT2D_config)
    return lk_Func(Array{Complex{Float64}, 1}(undef, c.nx), c)
end


function ==(f::lk_Func, g::lk_Func)
    return f.val == g.val && f.config == g.config
end


function copy(f::lk_Func)
    return lk_Func(f.vals, f.config)
end

# ----- OPERATORS -----
# + operator
+(f::lk_Func) = f

function +(f::lk_Func, g::lk_Func)
    if f.config != g.config
        error("FFT2D_config of two lk_Func are not equal.")
    else
        return lk_Func(f.vals + g.vals, f.config)
    end
end

+(f::lk_Func, a::Number) = lk_Func(k.vals .+ a, f.config)
+(a::Number, f::lk_Func) = f + a

# - operator
-(f::lk_Func) = lk_Func(- f.vals, f.config)

function -(f::lk_Func, g::lk_Func)
    if f.config != g.config
        error("FFT2D_config of two lk_Func are not equal.")
    else
        return lk_Func(f.vals - g.vals, f.config)
    end
end

-(f::lk_Func, a::Number) = lk_Func(f.vals .- a, f.config)
-(a::Number, f::lk_Func) = lk_Func(a .- f.vals, f.config)

# * operator
function *(f::lk_Func, g::lk_Func)
    if f.config != g.config
        error("FFT2D_config of two lk_Func are not equal.")
    else
        return lk_Func(f.vals .* g.vals, f.config)
    end
end

*(f::lk_Func, a::Number) = lk_Func(f.vals * a, f.config)
*(a::Number, f::lk_Func) = lk_Func(a * f.vals, f.config)

# / operator
function /(f::lk_Func, g::lk_Func)
    if f.config != g.config
        error("FFT2D_config of two lk_Func are not equal.")
    else
        return lk_Func(f.vals ./ g.vals, f.config)
    end
end

/(f::lk_Func, a::Number) = lk_Func(f.vals / a, f.config)
/(a::Number, f::lk_Func) = lk_Func(a ./ f.vals, f.config)


# ----- FOURIER TRANSFORMATION -----
function lk_yx(f::yx_Func)::lk_Func
    return lk_Func(fft(f.vals), f.config)
end

function yx_lk(f::lk_Func)::yx_Func

    if f.config.dealiasing
        f = copy(f)
	L = f.config.nx
	f.vals[div(L, 3) : div(2*L + 2, 3)] .= 0.0
    end
    return yx_Func(real(ifft(f.vals)), f.config)
end


# ******************
# *** Derivative ***
# ******************
function lk_Dx_lk(lk_func::lk_Func)::lk_Func
    c = lk_func.config
    lk_K = lk_Func(c.lk_K, c)
    return im * lk_K * lk_func
end


function lk_Dy_lk(lk_func::lk_Func)::lk_Func
    c = lk_func.config
    lk_L = lk_Func(c.lk_L, c)
    return im * lk_L * lk_func
end

# ****************************
# *** Elementary For ***
# ****************************
# Functions yx_Func
sin(f::yx_Func) = yx_Func(sin.(f.vals), f.config)
cos(f::yx_Func) = yx_Func(cos.(f.vals), f.config)
tan(f::yx_Func) = yx_Func(tan.(f.vals), f.config)
cot(f::yx_Func) = yx_Func(cot.(f.vals), f.config)
sec(f::yx_Func) = yx_Func(sec.(f.vals), f.config)
csc(f::yx_Func) = yx_Func(csc.(f.vals), f.config)

asin(f::yx_Func) = yx_Func(asin.(f.vals), f.config)
acos(f::yx_Func) = yx_Func(acos.(f.vals), f.config)
atan(f::yx_Func) = yx_Func(atan.(f.vals), f.config)
acot(f::yx_Func) = yx_Func(acot.(f.vals), f.config)
asec(f::yx_Func) = yx_Func(asec.(f.vals), f.config)
acsc(f::yx_Func) = yx_Func(acsc.(f.vals), f.config)

sinh(f::yx_Func) = yx_Func(sinh.(f.vals), f.config)
cosh(f::yx_Func) = yx_Func(cosh.(f.vals), f.config)
tanh(f::yx_Func) = yx_Func(tanh.(f.vals), f.config)
coth(f::yx_Func) = yx_Func(coth.(f.vals), f.config)
sech(f::yx_Func) = yx_Func(sech.(f.vals), f.config)
csch(f::yx_Func) = yx_Func(csch.(f.vals), f.config)

asinh(f::yx_Func) = yx_Func(asinh.(f.vals), f.config)
acosh(f::yx_Func) = yx_Func(acosh.(f.vals), f.config)
atanh(f::yx_Func) = yx_Func(atanh.(f.vals), f.config)
acoth(f::yx_Func) = yx_Func(acoth.(f.vals), f.config)
asech(f::yx_Func) = yx_Func(asech.(f.vals), f.config)
acsch(f::yx_Func) = yx_Func(acsch.(f.vals), f.config)

exp(f::yx_Func) = yx_Func(exp.(f.vals), f.config)
log(f::yx_Func) = yx_Func(log.(f.vals), f.config)

# For lk_Func
sin(f::lk_Func) = lk_Func(sin.(f.vals), f.config)
cos(f::lk_Func) = lk_Func(cos.(f.vals), f.config)
tan(f::lk_Func) = lk_Func(tan.(f.vals), f.config)
cot(f::lk_Func) = lk_Func(cot.(f.vals), f.config)
sec(f::lk_Func) = lk_Func(sec.(f.vals), f.config)
csc(f::lk_Func) = lk_Func(csc.(f.vals), f.config)

asin(f::lk_Func) = lk_Func(asin.(f.vals), f.config)
acos(f::lk_Func) = lk_Func(acos.(f.vals), f.config)
atan(f::lk_Func) = lk_Func(atan.(f.vals), f.config)
acot(f::lk_Func) = lk_Func(acot.(f.vals), f.config)
asec(f::lk_Func) = lk_Func(asec.(f.vals), f.config)
acsc(f::lk_Func) = lk_Func(acsc.(f.vals), f.config)

sinh(f::lk_Func) = lk_Func(sinh.(f.vals), f.config)
cosh(f::lk_Func) = lk_Func(cosh.(f.vals), f.config)
tanh(f::lk_Func) = lk_Func(tanh.(f.vals), f.config)
coth(f::lk_Func) = lk_Func(coth.(f.vals), f.config)
sech(f::lk_Func) = lk_Func(sech.(f.vals), f.config)
csch(f::lk_Func) = lk_Func(csch.(f.vals), f.config)

asinh(f::lk_Func) = lk_Func(asinh.(f.vals), f.config)
acosh(f::lk_Func) = lk_Func(acosh.(f.vals), f.config)
atanh(f::lk_Func) = lk_Func(atanh.(f.vals), f.config)
acoth(f::lk_Func) = lk_Func(acoth.(f.vals), f.config)
asech(f::lk_Func) = lk_Func(asech.(f.vals), f.config)
acsch(f::lk_Func) = lk_Func(acsch.(f.vals), f.config)

exp(f::lk_Func) = lk_Func(exp.(f.vals), f.config)
log(f::lk_Func) = lk_Func(log.(f.vals), f.config)


# ********************
# *** test program ***
# ********************
function test_FFT2D(;dealiasing=true)

    # *** parameter setting ***
    nx = 50	    # # of grid points
    ny = 50	    # # of grid points
    xmin = 0.0	    # minimum value of x
    xmax = 2pi	    # maximum value of x
    ymin = 0.0	    # minimum value of y
    ymax = 2pi	    # maximum value of y

    c = configure_FFT2D(nx=nx, ny=ny,
      			xmin=xmin, xmax=xmax,
			ymin=ymin, ymax=ymax,
      			dealiasing=dealiasing)
    yx_X = yx_Func(c.yx_X, c) # x-coordinate
    yx_Y = yx_Func(c.yx_Y, c) # y-coordinate

    # *** differentiate by spectral method ***
    yx_func = exp(sin(yx_X)) * cos(yx_Y)
    lk_func = lk_yx(yx_func)
    lk_dfunc = lk_Dy_lk(lk_Dx_lk(lk_func))
    yx_dfunc = yx_lk(lk_dfunc)

    # *** exact derivative ***
    yx_exact = - cos(yx_X) * exp(sin(yx_X)) * sin(yx_Y)

    # *** output ***
    open("testcase.txt", "w") do f
        for ix = 1:nx
	    for iy = 1:ny
                println(f,
		    yx_X.vals[iy, ix], " ",
		    yx_Y.vals[iy, ix], " ",
		    yx_func.vals[iy, ix], " ",
		    yx_dfunc.vals[iy, ix], " ",
		    yx_exact.vals[iy, ix])
	        # OUTPUT:
	        #   x, y, f(x, y),
		#   f'(x, y)(numerical), f'(x, y)(exact)
	    end
	    println(f)
	end
	println(f)
    end

end