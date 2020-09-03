# fft version
using FFTW

# *******************************************
#  Configuration
# *******************************************
struct ConfigFFT{N}

    ngrids::Tuple{Vararg{Int, N}} # Number Of Grids
    nwaves::Tuple{Vararg{Int, N}} # Cutoff wavenumber,
    		                  # ngrids[_] >= nwaves[_]
    xranges::Tuple{Vararg{Tuple{Float64, Float64}, N}}
    				  # tuple of (min, max) s
    Xcoords::Tuple{Vararg{Array{Float64, N}, N}}
    				  # X-coordinates
    Kcoords::Tuple{Vararg{Array{Complex{Float64}, N}, N}}
    				  # K-coordinates

    # CONSTRUCTOR
    function ConfigFFT{N}(
        ngrids::Tuple{Vararg{Int, N}},
        nwaves::Tuple{Vararg{Int, N}},
	xranges::Tuple{
	             Vararg{Tuple{Float64, Float64}, N}
		 },
	Xcoords::Tuple{
		     Vararg{Array{Float64, N}, N}
		 },
	Kcoords::Tuple{
		     Vararg{Array{Complex{Float64}, N}, N}
		 }) where N
		 
        if all(ngrids .>= nwaves)
	    return new(ngrids, nwaves, xranges, Xcoords, Kcoords)
	else
	    println("ERROR")
	end
    end

    function ConfigFFT(
        ngrids::Tuple{Vararg{Int, N}},
        nwaves::Tuple{Vararg{Int, N}},
	xranges::Tuple{
	             Vararg{Tuple{Float64, Float64}, N}
		 }) where N

	carts = CartesianIndices(ngrids)
	Xcoords = tuple(
	    (
	        (indices
		     -> xval(indices, dim, xranges[dim],
		             dxgen(ngrids[dim], xranges[dim]))
		).(carts)
		for dim = 1:N
	    )...
	)
	Kcoords = tuple(
	    (
	        (indices
		     -> kval(indices, dim, ngrids[dim])
		).(carts)
		for dim = 1:N
	    )...
	)

	return ConfigFFT{N}(
	           ngrids, nwaves, xranges, Xcoords, Kcoords)
    end

end

# --- helper function for constuctor ---
function xval(indices, dim, xrange, dx)
    return xrange[1] + (indices[dim] - 1)*dx
end

function dxgen(ngrid, xrange)
    return (xrange[2] - xrange[1]) / ngrid
end

function kval(indices, dim, ngrid)
    k = indices[dim]
    if k <= div(ngrid, 2)
        return Complex{Float64}(k - 1)
    else
        return Complex{Float64}(k - ngrid - 1)
    end
end


# *******************************************
#  XFunc, KFunc
# *******************************************
# +++++ XFunc +++++++++++++++++++++++++
mutable struct XFunc{N} <: AbstractArray{Float64, N}

    vals::Array{Float64, N}
    config::ConfigFFT{N}

    # constructor
    function XFunc{N}(vals::Array{Float64, N},
                      config::ConfigFFT{N}) where N
        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end
    end

    function XFunc(vals::Array{T, N},
    	           config::ConfigFFT{N}
		  ) where N where T <: Real
        return XFunc{N}(float(vals), config)
    end

    function XFunc(undef::UndefInitializer,
		   config::ConfigFFT{N}) where N
	f_undef = Array{Float64, N}(undef, config.ngrids)
        return XFunc{N}(f_undef, config)
    end

end

Base.:size(f::XFunc) = size(f.vals)
Base.:getindex(f::XFunc, i...) = getindex(f.vals, i...)
Base.:setindex!(f::XFunc, v, i...) = setindex!(f.vals, v, i...)
Base.:copy(f::XFunc) = XFunc(copy(f.vals), f.config)


# +++++ KFunc +++++++++++++++++++++++++
mutable struct KFunc{N} <: AbstractArray{Complex{Float64}, N}

    vals::Array{Complex{Float64}, N}
    config::ConfigFFT{N}

    # constructor
    function KFunc{N}(vals::Array{Complex{Float64}, N},
                      config::ConfigFFT{N}) where N
        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end
    end

    function KFunc(vals::Array{T, N},
    	           config::ConfigFFT{N}
		  ) where N where T <: Number
        return KFunc{N}(complex(float(vals)), config)
    end

    function KFunc(undef::UndefInitializer,
                   config::ConfigFFT{N}) where N
        f_undef = Array{Complex{Float64}, N}(undef, config.nwaves)
	return KFunc{N}(f_undef, config)
    end

end

Base.:size(f::KFunc) = size(f.vals)
Base.:getindex(f::KFunc, i...) = getindex(f.vals, i...)
Base.:setindex!(f::KFunc, v, i...) = setindex!(f.vals, v, i...)
Base.:copy(f::KFunc) = KFunc(copy(f.vals), f.config)


# *******************************************
#  Binomial Operators
# *******************************************

BINOP = ((:+, :.+), (:-, :.-), (:*, :.*), (:/, :./), (:\, :.\))

# +++++ XFunc +++++
Base.:+(f::XFunc) = f
Base.:-(f::XFunc) = XFunc(-f.vals, f.config)

for (op, opd) = BINOP
    @eval begin
        function Base.$op(f::XFunc, g::XFunc)
	    if f.config == g.config
	        return XFunc($opd(f.vals, g.vals), f.config)
	    else
	        println("ERROR")
	    end
	end

	Base.$op(f::XFunc, a::Real) = (
			   XFunc($opd(f.vals, a), f.config) )
	Base.$op(a::Real, f::XFunc) = (
			   XFunc($opd(a, f.vals), f.config) )
    end
end


# +++++ KFunc +++++
Base.:+(f::KFunc) = f
Base.:-(f::KFunc) = KFunc(-f.vals, f.config)

for (op, opd) = BINOP
    @eval begin
        function Base.$op(f::KFunc, g::KFunc)
	    if f.config == g.config
	        return KFunc($opd(f.vals, g.vals), f.config)
	    else
	        println("ERROR")
	    end
	end

	Base.$op(f::KFunc, a::Number) = (
			   KFunc($opd(f.vals, a), f.config) )
	Base.$op(a::Number, f::KFunc) = (
			   KFunc($opd(a, f.vals), f.config) )
    end
end


# *******************************************
#  Coordinate Tools
# *******************************************
# +++++ Coordinate in X-space ++++++++++
function x_Xgen(config::ConfigFFT{1})
    return XFunc(copy(config.Xcoords[1]), config)
end

function xy_Xgen(config::ConfigFFT{2})
    return XFunc(copy(config.Xcoords[1]), config)
end

function xy_Ygen(config::ConfigFFT{2})
    return XFunc(copy(config.Xcoords[2]), config)
end

function xyz_Xgen(config::ConfigFFT{3})
    return XFunc(copy(config.Xcoords[1]), config)
end

function xyz_Ygen(config::ConfigFFT{3})
    return XFunc(copy(config.Xcoords[2]), config)
end

function xyz_Zgen(config::ConfigFFT{3})
    return XFunc(copy(config.Xcoords[3]), config)
end


# +++++ Coordinate in K-space ++++++++++
function k_Kgen(config::ConfigFFT{1})
    return KFunc(copy(config.Kcoords[1]), config)
end

function kl_Kgen(config::ConfigFFT{2})
    return KFunc(copy(config.Kcoords[1]), config)
end

function kl_Lgen(config::ConfigFFT{2})
    return KFunc(copy(config.Kcoords[2]), config)
end

function klm_Kgen(config::ConfigFFT{3})
    return KFunc(copy(config.Kcoords[1]), config)
end

function klm_Lgen(config::ConfigFFT{3})
    return KFunc(copy(config.Kcoords[2]), config)
end

function klm_Mgen(config::ConfigFFT{3})
    return KFunc(copy(config.Kcoords[3]), config)
end


# *******************************************
#  Elementary Functions
# *******************************************
ELEMFUNC = (
    :sin, :cos, :tan, :cot, :sec, :csc,
    :sinh, :cosh, :tanh, :coth, :sech, :csch,
    :asin, :acos, :atan, :acot, :asec, :acsc,
    :asinh, :acosh, :atanh, :acoth, :asech, :acsch,
    :sinpi, :cospi, :sinc, :cosc,
    :exp, :log
)

for fn = ELEMFUNC
    @eval begin
        Base.$fn(X_F::XFunc) = XFunc($fn.(X_F.vals), X_F.config)
	Base.$fn(K_F::KFunc) = KFunc($fn.(K_F.vals), K_F.config)
    end
end

# *******************************************
#  Low/High-pass Filter
# *******************************************

# To Be Implemented!


# *******************************************
#  Fourier Transformation
# *******************************************
K_X(f::XFunc) = KFunc(fft(f.vals), f.config)
X_K(f::KFunc) = XFunc(real(ifft(f.vals)), f.config)

k_x(f::XFunc{1}) = K_X(f)
x_k(f::KFunc{1}) = X_K(f)

kl_xy(f::XFunc{2}) = K_X(f)
xy_kl(f::KFunc{2}) = X_K(f)

klm_xyz(f::XFunc{3}) = K_X(f)
xyz_klm(f::KFunc{3}) = X_K(f)


# *******************************************
#  Differentiation
# *******************************************
# +++++ destructive +++++
function DXd_K!(K_func::KFunc{N}, d::Int) where N

    vals = K_func.vals
    config = K_func.config
    ngrid = config.ngrids[d]
    xrange = config.xranges[d]
    xlen = xrange[2] - xrange[1]
    K_Kd = config.Kcoords[d]

    if ngrid % 2 == 0
        cut_highest_wave!(vals, d, ngrid)
    end
    vals .*= (2π*im/xlen)*K_Kd

end

function cut_highest_wave!(vals, d, ngrid)

    slice = [
        n == d ? div(ngrid, 2) : Colon()
	for n = 1:ndims(vals)
    ]
    if ndims(vals) == 1
        vals[slice...] = 0.
    else
        vals[slice...] .= 0.
    end

end

Dx_k!(k_func::KFunc{1}) = DXd_K!(k_func, 1)

Dx_kl!(kl_func::KFunc{2}) = DXd_K!(kl_func, 1)
Dy_kl!(kl_func::KFunc{2}) = DXd_K!(kl_func, 2)

Dx_klm!(klm_func::KFunc{3}) = DXd_K!(klm_func, 1)
Dy_klm!(klm_func::KFunc{3}) = DXd_K!(klm_func, 2)
Dz_klm!(klm_func::KFunc{3}) = DXd_K!(klm_func, 3)


# +++++ non-destructive +++++
function K_DXd_K(K_func::KFunc{N}, d::Int) where N

    K_func2 = copy(K_func)
    DXd_K!(K_func2, d)
    return K_func2

end

k_Dx_k(k_func::KFunc{1}) = K_DXd_K(k_func, 1)

kl_Dx_kl(kl_func::KFunc{2}) = K_DXd_K(kl_func, 1)
kl_Dy_kl(kl_func::KFunc{2}) = K_DXd_K(kl_func, 2)

klm_Dx_klm(klm_func::KFunc{3}) = K_DXd_K(klm_func, 1)
klm_Dy_klm(klm_func::KFunc{3}) = K_DXd_K(klm_func, 2)
klm_Dz_klm(klm_func::KFunc{3}) = K_DXd_K(klm_func, 3)