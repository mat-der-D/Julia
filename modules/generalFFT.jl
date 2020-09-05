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
function Xcoordsgen(
			ngrids::NTuple{N, Int},
			xranges::NTuple{N, NTuple{2, Float64}}
		) where N

    _Xcoordgen(axis) = Xcoordgen(axis, ngrids, xranges)
    return ntuple(_Xcoordgen, N)

end

function Xcoordgen(
			axis::Int,
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

function Kcoordsgen(
			ngrids::NTuple{N, Int},
			xranges::NTuple{N, NTuple{2, Float64}}
		) where N

	_Kcoordgen(axis) = Kcoordgen(axis, ngrids, xranges)
    return ntuple(_Kcoordgen, N)

end

function Kcoordgen(
			axis::Int,
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
    function XFunc{N}(
				vals::Array{Float64, N},
    			config::ConfigFFT{N}
			) where N

        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end

    end

    # EASY CONSTRUCTOR
    function XFunc(
				vals::Array{T, N},
				config::ConfigFFT{N}
			) where N where T <: Real

        XFunc{N}(float(vals), config)

    end

    # UNDEF CONSTRUCTOR
    function XFunc(
				undef::UndefInitializer,
				config::ConfigFFT{N}
				) where N

		f_undef = Array{Float64, N}(undef, config.ngrids)
        return XFunc{N}(f_undef, config)

    end

end

Base.:size(f::XFunc) = size(f.vals)

Base.:getindex(f::XFunc, i::Int) = getindex(f.vals, i)
function Base.:getindex(f::XFunc{N}, I::Vararg{Int, N}) where N
	getindex(f.vals, I...)
end

Base.:setindex!(f::XFunc, v, i::Int) = setindex!(f.vals, v, i)
function Base.:setindex!(f::XFunc{N}, v, I::Vararg{Int, N}) where N
	setindex!(f, v, I...)
end

Base.:copy(f::XFunc) = XFunc(copy(f.vals), f.config)


# +++++ KFunc +++++++++++++++++++++++++
mutable struct KFunc{N} <: AbstractArray{Complex{Float64}, N}

    vals::Array{Complex{Float64}, N}
    config::ConfigFFT{N}

    # CONSTRUCTOR
    function KFunc{N}(
				vals::Array{Complex{Float64}, N},
            	config::ConfigFFT{N}
			) where N

        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end

    end

    # EASY CONSTRUCTOR
    function KFunc(
				vals::Array{T, N},
				config::ConfigFFT{N}
			) where N where T <: Number
        KFunc{N}(complex(float(vals)), config)
    end

    # UNDEF CONSTRUCTOR
    function KFunc(
				undef::UndefInitializer,
				config::ConfigFFT{N}
			) where N

		f_undef = Array{Complex{Float64}, N}(undef, config.ngrids)
		return KFunc{N}(f_undef, config)

    end

end

Base.:size(f::KFunc) = size(f.vals)

Base.:getindex(f::KFunc, i::Int) = getindex(f.vals, i)
function Base.:getindex(f::KFunc{N}, I::Vararg{Int, N}) where N
	getindex(f.vals, I...)
end

Base.:setindex!(f::KFunc, v, i::Int) = setindex!(f.vals, v, i)
function Base.:setindex!(f::KFunc{N}, v, I::Vararg{Int, N}) where N
	setindex!(f, v, I...)
end

Base.:copy(f::KFunc) = KFunc(copy(f.vals), f.config)


# *******************************************
#  Binomial Operators
# *******************************************
BINOP = (
	(:+, :.+), (:-, :.-), (:*, :.*),
	(:/, :./), (:\, :.\), (:^, :.^)
)

# +++++ XFunc +++++
Base.:+(f::XFunc) = f
Base.:-(f::XFunc) = XFunc(-f.vals, f.config)
Base.:^(f::XFunc, n::Integer) = XFunc(f.vals .^ n, f.config)
Base.:inv(f::XFunc) = f \ 1.0

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
			XFunc($opd(f.vals, a), f.config)
		)
		Base.$op(a::Real, f::XFunc) = (
			   XFunc($opd(a, f.vals), f.config)
		)
    end
end



# +++++ KFunc +++++
Base.:+(f::KFunc) = f
Base.:-(f::KFunc) = KFunc(-f.vals, f.config)
Base.:^(f::KFunc, n::Integer) = KFunc(f.vals .^ n, f.config)
Base.:inv(f::KFunc) = f \ 1.0

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
			KFunc($opd(f.vals, a), f.config)
		)
		Base.$op(a::Number, f::KFunc) = (
			KFunc($opd(a, f.vals), f.config)
		)
    end
end


# *******************************************
#  Coordinate Tools for specific dimensions
# *******************************************
# +++++ Coordinate in X-space ++++++++++
function x_Xgen(config::ConfigFFT{1})
    XFunc(copy(config.Xcoords[1]), config)
end

function xy_Xgen(config::ConfigFFT{2})
    XFunc(copy(config.Xcoords[1]), config)
end

function xy_Ygen(config::ConfigFFT{2})
    XFunc(copy(config.Xcoords[2]), config)
end

function xyz_Xgen(config::ConfigFFT{3})
    XFunc(copy(config.Xcoords[1]), config)
end

function xyz_Ygen(config::ConfigFFT{3})
    XFunc(copy(config.Xcoords[2]), config)
end

function xyz_Zgen(config::ConfigFFT{3})
    XFunc(copy(config.Xcoords[3]), config)
end


# +++++ Coordinate in K-space ++++++++++
function k_Kgen(config::ConfigFFT{1})
    KFunc(copy(config.Kcoords[1]), config)
end

function kl_Kgen(config::ConfigFFT{2})
    KFunc(copy(config.Kcoords[1]), config)
end

function kl_Lgen(config::ConfigFFT{2})
    KFunc(copy(config.Kcoords[2]), config)
end

function klm_Kgen(config::ConfigFFT{3})
    KFunc(copy(config.Kcoords[1]), config)
end

function klm_Lgen(config::ConfigFFT{3})
    KFunc(copy(config.Kcoords[2]), config)
end

function klm_Mgen(config::ConfigFFT{3})
    KFunc(copy(config.Kcoords[3]), config)
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
        Base.$fn(f::XFunc) = XFunc($fn.(f.vals), f.config)
		Base.$fn(f::KFunc) = KFunc($fn.(f.vals), f.config)
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
function pass_K!(
			f::KFunc{N},
			slices::NTuple{N, UnitRange{Int}}
		) where N

    vals = copy(f.vals)
    vals[slices...] .= 0.0 + 0.0im
    f.vals -= vals

end


# +++++ high-pass filter +++++
function highpass_K!(
			f::KFunc{N},
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
		max_indices = @. min(ceils, ceils + 1 - min_nwaves)

		slices = make_slice.(min_indices, max_indices)
		pass_K!(f, slices)
    end
end

function K_highpass_K(
			f::KFunc{N},
			min_nwaves::NTuple{N, Int}
		) where N

    g = copy(f)
    highpass_K!(g, min_nwaves)
    return g

end


# +++++ low-pass filter +++++
function lowpass_K!(
			f::KFunc{N},
			max_nwaves::NTuple{N, Int}
		) where N

    ngrids = f.config.ngrids
    nshifts = @. (ngrids - 1) ÷ 2
    circshift!(f.vals, copy(f.vals), nshifts) # SHIFT

    center_indices = @. (ngrids + 1) ÷ 2

    if any(max_nwaves .< 0)
        println("WARNING: all waves are suppressed")
		f.vals .= 0.
    else
        floors = tuple(ones(Int, N)...)
		ceils = ngrids
		min_indices = (
			@. max(floors, center_indices - max_nwaves)
		)
		max_indices = (
			@. min(ceils, center_indices + max_nwaves)
		)

		slices = make_slice.(min_indices, max_indices)
		pass_K!(f, slices)
    end

    ndeshifts = .-(tuple(nshifts...))
    circshift!(f.vals, copy(f.vals), ndeshifts) # DESHIFT

end

function K_lowpass_K(
			f::KFunc{N},
			min_nwaves::NTuple{N, Int}
		) where N

    g = copy(f)
    lowpass_K!(g, min_nwaves)
    return g

end

# helper funtcion
make_slice(x::Int, y::Int) = x:y


# +++++ aliases for specific dimensins +++++
# @@@ high-pass filter @@@
function k_highpass_k(f::KFunc{1}, min_nwaves::Tuple{Int})
    K_highpass_K(f, min_nwaves)
end

function kl_highpass_kl(f::KFunc{2}, min_nwaves::NTuple{2, Int})
    K_highpass_K(f, min_nwaves)
end

function klm_highpass_klm(f::KFunc{3}, min_nwaves::NTuple{3, Int})
    K_highpass_K(f, min_nwaves)
end

# @@@ low-pass filter @@@
function k_lowpass_k(f::KFunc{1}, max_nwaves::Tuple{Int})
    K_lowpass_K(f, max_nwaves)
end

function kl_lowpass_kl(f::KFunc{2}, max_nwaves::NTuple{2, Int})
    K_lowpass_K(f, max_nwaves)
end

function klm_lowpass_klm(f::KFunc{3}, max_nwaves::NTuple{3, Int})
    K_lowpass_K(f, max_nwaves)
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
function ∂Xaxis_K!(f::KFunc, axis::Int)

    config = f.config
    xrange = config.xranges[axis]
    xlen = xrange[2] - xrange[1]
    Kcoord = config.Kcoords[axis]

    f.vals .*= (2π*im/xlen)*Kcoord

end

∂x_k!(k_func::KFunc{1}) = ∂Xaxis_K!(k_func, 1)

∂x_kl!(kl_func::KFunc{2}) = ∂Xaxis_K!(kl_func, 1)
∂y_kl!(kl_func::KFunc{2}) = ∂Xaxis_K!(kl_func, 2)

∂x_klm!(klm_func::KFunc{3}) = ∂Xaxis_K!(klm_func, 1)
∂y_klm!(klm_func::KFunc{3}) = ∂Xaxis_K!(klm_func, 2)
∂z_klm!(klm_func::KFunc{3}) = ∂Xaxis_K!(klm_func, 3)

# +++++ non-destructive +++++
function K_∂Xaxis_K(f::KFunc, axis::Int)

    g = copy(f)
    ∂Xaxis_K!(g, axis)
    return g

end

X_∂Xaxis_X = X_K ∘ K_∂Xaxis_X ∘ K_X

k_∂x_k(k_func::KFunc{1}) = K_∂Xaxis_K(k_func, 1)
x_∂x_x(x_func::XFunc{1}) = X_∂Xaxis_X(x_func, 1)

kl_∂x_kl(kl_func::KFunc{2}) = K_∂Xaxis_K(kl_func, 1)
kl_∂y_kl(kl_func::KFunc{2}) = K_∂Xaxis_K(kl_func, 2)
xy_∂x_xy(xy_func::XFunc{2}) = X_∂Xaxis_X(xy_func, 1)
xy_∂y_xy(xy_func::XFunc{2}) = X_∂Xaxis_X(xy_func, 2)

klm_∂x_klm(klm_func::KFunc{3}) = K_∂Xaxis_K(klm_func, 1)
klm_∂y_klm(klm_func::KFunc{3}) = K_∂Xaxis_K(klm_func, 2)
klm_∂z_klm(klm_func::KFunc{3}) = K_∂Xaxis_K(klm_func, 3)
xyz_∂x_xyz(xyz_func::KFunc{3}) = X_∂Xaxis_X(xyz_func, 1)
xyz_∂y_xyz(xyz_func::KFunc{3}) = X_∂Xaxis_X(xyz_func, 2)
xyz_∂z_xyz(xyz_func::KFunc{3}) = X_∂Xaxis_X(xyz_func, 3)

# +++ tools for vector analysis +++
# 2-dimensional
function kl2_grad_kl(kl_func::KFunc{2})::Vector{KFunc{2}}
	return [
		kl_∂x_kl(kl_func)
		kl_∂y_kl(kl_func)
	]
end

function kl_rot_kl2(kl2_func::Vector{KFunc{2}})::KFunc{2}

	if length(kl2_func) != 2
		return println("ERROR")
	end
	return (
		kl_∂x_kl(kl2_func[2])
		- kl_∂y_kl(kl2_func[1])
	)
end

function kl_div_kl2(kl2_func::Vector{KFunc{2}})::KFunc{2}

	if length(kl2_func) != 2
		return println("ERROR")
	end
	return (
		kl_∂x_kl(kl2_func[1])
		+ kl_∂y_kl(kl2_func[2])
	)
end

# 3-dimensional
function klm3_grad_klm(
			klm_func::KFunc{3}
		)::Vector{KFunc{3}}

	return [
		klm_∂x_klm(klm_func)
		klm_∂y_klm(klm_func)
		klm_∂z_klm(klm_func)
	]

end

function klm3_rot_klm(
			klm3_func::Vector{KFunc{3}}
		)::Vector{KFunc{3}}

	if length(klm3_func) != 3
		return println("ERROR")
	end if
	return [
		klm_∂y_klm(klm3_func[3]) - klm_∂z_klm(klm3_func[2])
		klm_∂z_klm(klm3_func[1]) - klm_∂x_klm(klm3_func[3])
		klm_∂x_klm(klm3_func[2]) - klm_∂y_klm(klm3_func[1])
	]

end

function klm_div_klm3(
			klm3_func::Vector{KFunc{3}}
		)::KFunc{3}

	if length(klm3_func) != 3
		return println("ERROR")
	end
	return (
		klm_∂x_klm(klm3_func[1])
		+ klm_∂y_klm(klm3_func[2])
		+ klm_∂z_klm(klm3_func[3])
	)

end

# *******************************************
#  Analysis Tools
# *******************************************
function dVgen(config::ConfigFFT)
    xlens = (x -> -(.-(x...))).(config.xranges)
    return prod(xlens ./ config.ngrids)
end

function integ_X(f::XFunc)
    sum(f) * dVgen(f.config)
end

∫(f) = integ_X(f)

function norm_X(f::XFunc, p::Real=2)
	if p == Inf
		return max(abs(f)...)
	else
    	return ( ∫( abs(f)^p ) )^(1/p)
	end
end

function l2inpr_X_X(f::XFunc{N}, g::XFunc{N}) where N
    if f.config === g.config
		return ∫(f * g)
    else
        println("ERROR")
    end
end

# aliaces
integ_x(x_func::XFunc{1}) = integ_X(x_func)
integ_xy(xy_func::XFunc{2}) = integ_X(xy_func)
integ_xyz(xyz_func::XFunc{3}) = integ_X(xyz_func)

norm_x(x_func::XFunc{1}) = norm_X(x_func)
norm_xy(xy_func::XFunc{2}) = norm_X(xy_func)
norm_xyz(xyz_func::XFunc{3}) = norm_X(xyz_func)

l2inpr_x_x(x_f::XFunc{1}, x_g::XFunc{1}) = l2inpr_X_X(x_f, x_g)
l2inpr_xy_xy(xy_f::XFunc{2}, xy_g::XFunc{2}) = l2inpr_X_X(xy_f, xy_g)
l2inpr_xyz_xyz(xyz_f::XFunc{3}, xyz_g::XFunc{3}) = l2inpr_X_X(xyz_f, xyz_g)
