# fft version
using FFTW

# *******************************************
#  Configuration
# *******************************************
struct ConfigFFT{N}

    ngrids::NTuple{N, Int}        # Number Of Grids
    nwaves::NTuple{N, Int}   	  # Cutoff wavenumber
    
    xranges::NTuple{N, NTuple{2, Float64}}
    				  # tuple of (min, max) s
    Xcoords::NTuple{N, Array{Float64, N}}
    				  # X-coordinates
    Kcoords::NTuple{N, Array{Complex{Float64}, N}}
    				  # K-coordinates

    # CONSTRUCTOR
    function ConfigFFT{N}(
        ngrids, nwaves, xranges, Xcoords, Kcoords) where N
		 
        if all(@. ngrids ÷ 2 >= nwaves)
	    return new(ngrids, nwaves, xranges, Xcoords, Kcoords)
	else
	    println("ERROR")
	end
    end

    # EASY CONSTRUCTOR
    function ConfigFFT(
        ngrids::NTuple{N, Int},
	nwaves::NTuple{N, Int},
	xranges::NTuple{N, NTuple{2, Float64}}
	) where N

	carts = CartesianIndices(ngrids)
	Xcoords = Xcoordsgen(ngrids, xranges)
	Kcoords = Kcoordsgen(ngrids, xranges)

	return ConfigFFT{N}(
	           ngrids, nwaves, xranges, Xcoords, Kcoords)
    end

end

# --- helper function for constuctor ---
function Xcoordsgen(ngrids::NTuple{N, Int},
	 	    xranges::NTuple{N, NTuple{2, Float64}}
		   ) where N
    _Xcoordgen(axis) = Xcoordgen(axis, ngrids, xranges)
    return ntuple(_Xcoordgen, N)
end

function Xcoordgen(axis::Int,
		   ngrids::NTuple{N, Int},
		   xranges::NTuple{N, NTuple{2, Float64}}
		  ) where N
    ngrid = ngrids[axis]
    xrange = xranges[axis]
    _Xcoordgen(indices) = (
        (  (indices[axis] - 1)*xrange[2]
	 - (indices[axis] - 2)*xrange[1] ) / ngrid
    )
    return _Xcoordgen.(CartesianIndices(ngrids))
end

function Kcoordsgen(ngrids::NTuple{N, Int},
	 	    xranges::NTuple{N, NTuple{2, Float64}}
		   ) where N
    _Kcoordgen(axis) = Kcoordgen(axis, ngrids, xranges)
    return ntuple(_Kcoordgen, N)
end

function Kcoordgen(axis::Int,
		   ngrids::NTuple{N, Int},
		   xranges::NTuple{N, NTuple{2, Float64}}
		  ) where N
    ngrid = ngrids[axis]
    xrange = xranges[axis]
    _Kcoordgen(indices) = (
        kval(indices[axis], ngrid, xrange)
    )
    return _Kcoordgen.(CartesianIndices(ngrids))
end

function kval(index, ngrid, xrange)::Complex{Float64}
    index0 = index - 1  # start by 0
    if 2*index0 < ngrid
        return index0
    elseif 2*index0 == ngrid
        return 0
    elseif 2*index0 > ngrid
        return index0 - ngrid
    end
end


# *******************************************
#  XFunc, KFunc
# *******************************************
# +++++ XFunc +++++++++++++++++++++++++
mutable struct XFunc{N} <: AbstractArray{Float64, N}

    vals::Array{Float64, N}
    config::ConfigFFT{N}

    # CONSTRUCTOR
    function XFunc{N}(vals::Array{Float64, N},
                      config::ConfigFFT{N}) where N
        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end
    end

    # EASY CONSTRUCTOR
    function XFunc(vals::Array{T, N},
    	           config::ConfigFFT{N}
		  ) where N where T <: Real
        return XFunc{N}(float(vals), config)
    end

    # UNDEF CONSTRUCTOR
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

    # CONSTRUCTOR
    function KFunc{N}(vals::Array{Complex{Float64}, N},
                      config::ConfigFFT{N}) where N
        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end
    end

    # EASY CONSTRUCTOR
    function KFunc(vals::Array{T, N},
    	           config::ConfigFFT{N}
		  ) where N where T <: Number
        return KFunc{N}(complex(float(vals)), config)
    end

    # UNDEF CONSTRUCTOR
    function KFunc(undef::UndefInitializer,
                   config::ConfigFFT{N}) where N
        f_undef = Array{Complex{Float64}, N}(undef, config.grids)
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

BINOP = ((:+, :.+), (:-, :.-), (:*, :.*),
         (:/, :./), (:\, :.\), (:^, :.^))

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
#  Coordinate Tools for specific dimensions
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
    :exp, :log, :sqrt, :cbrt,
    :abs
)

for fn = ELEMFUNC
    @eval begin
        Base.$fn(X_F::XFunc) = XFunc($fn.(X_F.vals), X_F.config)
	Base.$fn(K_F::KFunc) = KFunc($fn.(K_F.vals), K_F.config)
    end
end

# *******************************************
#  Operators for Complex Numbers
# *******************************************
OPERATOR = (
    :real, :imag, :reim, :conj
)

for op = OPERATOR
    @eval Base.$op(K_F::KFunc) = KFunc($op(K_F.vals), K_F.config)
end


# *******************************************
#  Low/High-pass Filter
# *******************************************
# +++++ general pass filter +++++
function pass_K!(f::KFunc{N},
		 slices::NTuple{N, UnitRange{Int}}
		) where N
    vals = copy(f.vals)
    vals[slices...] .= 0.0 + 0.0im
    f.vals -= vals
end


# +++++ high-pass filter +++++
function highpass_K!(f::KFunc{N},
		     min_nwaves::NTuple{N, Int}
		    ) where N
    ngrids = f.config.ngrids
    if any(@. min_nwaves > ngrids ÷ 2)
        println("WARNING: all waves are suppressed")
	f.vals .= 0.
    else
        floors = tuple(ones(Int, N)...)
	ceils = ngrids
	min_indices = @. max(floors, floors + min_nwaves)
	max_indices = @. min(ceils, ceils - min_nwaves + 1)
	
	slices = ((x, y) -> x:y).(min_indices, max_indices)
	pass_K!(f, slices)
    end
end

function K_highpass_K(f::KFunc{N},
		      min_nwaves::NTuple{N, Int}
		     ) where N
    g = copy(f)
    highpass_K!(g, min_nwaves)
    return g
end


# +++++ low-pass filter +++++
function lowpass_K!(f::KFunc{N},
		    max_nwaves::NTuple{N, Int}
		   ) where N
		   
    ngrids = f.config.ngrids
    nshifts = @. (ngrids - 1) ÷ 2
    circshift!(f.vals, copy(f.vals), nshifts) # SHIFT!!!

    center_indices = @. (ngrids + 1) ÷ 2

    if any(max_nwaves .< 0)
        println("WARNING: all waves are suppressed")
	f.vals .= 0.
    else
        floors = tuple(ones(Int, N)...)
	ceils = ngrids
	min_indices = @. max(floors, center_indices - max_nwaves)
	max_indices = @. min(ceils, center_indices + max_nwaves)

	slices = ((x, y) -> x:y).(min_indices, max_indices)
	pass_K!(f, slices)
    end

    ndeshifts = .-(tuple(nshifts...))
    circshift!(f.vals, copy(f.vals), ndeshifts) # DESHIFT!!!

end

function K_lowpass_K(f::KFunc{N},
		     min_nwaves::NTuple{N, Int}
		    ) where N
    g = copy(f)
    lowpass_K!(g, min_nwaves)
    return g
end

# +++++ aliases for specific dimensins +++++
# @@@ high-pass filter @@@
function k_highpass_k(f::KFunc{1}, min_nwaves::Tuple{Int})
    return K_highpass_K(f, min_nwaves)
end

function kl_highpass_kl(f::KFunc{2}, min_nwaves::NTuple{2, Int})
    return K_highpass_K(f, min_nwaves)
end

function klm_highpass_klm(f::KFunc{3}, min_nwaves::NTuple{3, Int})
    return K_highpass_K(f, min_nwaves)
end

# @@@ low-pass filter @@@
function k_lowpass_k(f::KFunc{1}, max_nwaves::Tuple{Int})
    return K_lowpass_K(f, max_nwaves)
end

function kl_lowpass_kl(f::KFunc{2}, max_nwaves::NTuple{2, Int})
    return K_lowpass_K(f, max_nwaves)
end

function klm_lowpass_klm(f::KFunc{3}, max_nwaves::NTuple{3, Int})
    return K_lowpass_K(f, max_nwaves)
end


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

    vals .*= (2π*im/xlen)*K_Kd

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