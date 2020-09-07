# fft version
using FFTW

# *******************************************
#  Configuration
# *******************************************
struct ConfigFFT{T<:Union{Float64,Complex{Float64}},N}

    ngrids::NTuple{N,Int}
                    # Number Of Grids

    xranges::NTuple{N,NTuple{2,Float64}}
                    # tuple of (min, max) s
    Xcoords::NTuple{N,Array{Float64,N}}
                    # X-coordinates
    Kcoords::NTuple{N,Array{Complex{Float64},N}}
                    # K-coordinates
    cut_zigzag_mode::Bool

    # CONSTRUCTOR
    function ConfigFFT{T,N}(
                ngrids, xranges, Xcoords, Kcoords;
                cut_zigzag_mode=true
            ) where N where T
        new(ngrids, xranges, Xcoords, Kcoords,
            cut_zigzag_mode)
    end

    # EASY CONSTRUCTOR
    function ConfigFFT(
            ngrids::NTuple{N,Int},
            xranges::NTuple{N,NTuple{2,Float64}};
            use_complex::Bool=false,
            cut_zigzag_mode::Bool=true
        ) where N

        T = use_complex ? Complex{Float64} : Float64

        carts = CartesianIndices(ngrids)
        Xcoords = Xcoordsgen(ngrids, xranges)
        Kcoords = Kcoordsgen(ngrids, xranges)

        return ConfigFFT{T,N}(
            ngrids, xranges, Xcoords, Kcoords,
            cut_zigzag_mode=cut_zigzag_mode)
    end

end

# --- helper function for constuctor ---
function Xcoordsgen(
            ngrids::NTuple{N,Int},
            xranges::NTuple{N,NTuple{2,Float64}}
        ) where N

    _Xcoordgen(axis) = Xcoordgen(axis, ngrids, xranges)
    return ntuple(_Xcoordgen, N)

end

function Xcoordgen(
            axis::Int,
            ngrids::NTuple{N,Int},
            xranges::NTuple{N,NTuple{2,Float64}}
        ) where N

    ngrid = ngrids[axis]
    xrange = xranges[axis]
    _Xcoordgen(indices) = (
        (  (indices[axis] - 1)*xrange[2]
         + (ngrid - indices[axis] + 1)*xrange[1] ) / ngrid
    )
    return _Xcoordgen.(CartesianIndices(ngrids))

end

function Kcoordsgen(
            ngrids::NTuple{N,Int},
            xranges::NTuple{N,NTuple{2,Float64}}
        ) where N

    _Kcoordgen(axis) = Kcoordgen(axis, ngrids, xranges)
    return ntuple(_Kcoordgen, N)

end

function Kcoordgen(
            axis::Int,
            ngrids::NTuple{N,Int},
            xranges::NTuple{N,NTuple{2,Float64}}
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
mutable struct XFunc{T,N} <: AbstractArray{T,N}

    vals::Array{T,N}
    config::ConfigFFT{T,N}

    # CONSTRUCTOR
    function XFunc{T,N}(
                vals::Array{T,N},
                config::ConfigFFT{T,N}
            ) where N where T

        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end

    end

    # EASY CONSTRUCTOR
    function XFunc(
                vals::Array{Tv,N},
                config::ConfigFFT{Tc,N}
            ) where N where Tv <: Number where Tc

        XFunc{Tc,N}(Tc.(vals), config)

    end

    # UNDEF CONSTRUCTOR
    function XFunc(
                undef::UndefInitializer,
                config::ConfigFFT{T,N}
                ) where N where T

        f_undef = Array{T,N}(undef, config.ngrids)
        return XFunc{T,N}(f_undef, config)

    end

end

Base.:size(f::XFunc) = size(f.vals)

Base.:getindex(f::XFunc, i::Int) = getindex(f.vals, i)
function Base.:getindex(
            f::XFunc{T,N}, I::Vararg{Int,N}
        ) where N where T
    getindex(f.vals, I...)
end

Base.:setindex!(f::XFunc, v, i::Int) = setindex!(f.vals, v, i)
function Base.:setindex!(
            f::XFunc{T,N}, v, I::Vararg{Int,N}
        ) where N where T
    setindex!(f, v, I...)
end

Base.:copy(f::XFunc) = XFunc(copy(f.vals), f.config)


# +++++ KFunc +++++++++++++++++++++++++
mutable struct KFunc{T,N} <: AbstractArray{Complex{Float64},N}

    vals::Array{Complex{Float64},N}
    config::ConfigFFT{T,N}

    # CONSTRUCTOR
    function KFunc{T,N}(
                vals::Array{Complex{Float64},N},
                config::ConfigFFT{T,N}
            ) where N where T

        if size(vals) == config.ngrids
            return new(vals, config)
        else
            println("ERROR")
        end

    end

    # EASY CONSTRUCTOR
    function KFunc(
                vals::Array{Tv,N},
                config::ConfigFFT{Tc,N}
            ) where N where Tv <: Number where Tc
        KFunc{Tc,N}(complex(float(vals)), config)
    end

    # UNDEF CONSTRUCTOR
    function KFunc(
                undef::UndefInitializer,
                config::ConfigFFT{T,N}
            ) where N where T

        f_undef = Array{Complex{Float64},N}(undef, config.ngrids)
        return KFunc{T,N}(f_undef, config)

    end

end

Base.:size(f::KFunc) = size(f.vals)

Base.:getindex(f::KFunc, i::Int) = getindex(f.vals, i)
function Base.:getindex(
            f::KFunc{T,N}, I::Vararg{Int,N}
        ) where N where T
    getindex(f.vals, I...)
end

Base.:setindex!(f::KFunc, v, i::Int) = setindex!(f.vals, v, i)
function Base.:setindex!(
            f::KFunc{T,N}, v, I::Vararg{Int,N}
        ) where N where T
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
            if f.config === g.config
                return XFunc($opd(f.vals, g.vals), f.config)
            else
                println("ERROR")
            end
        end

        Base.$op(f::XFunc, a::Number) = (
            XFunc($opd(f.vals, a), f.config)
        )
        Base.$op(a::Number, f::XFunc) = (
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
function x_Xgen(config::ConfigFFT{T,1}) where T
    XFunc(copy(config.Xcoords[1]), config)
end

function xy_Xgen(config::ConfigFFT{T,2}) where T
    XFunc(copy(config.Xcoords[1]), config)
end

function xy_Ygen(config::ConfigFFT{T,2}) where T
    XFunc(copy(config.Xcoords[2]), config)
end

function xyz_Xgen(config::ConfigFFT{T,3}) where T
    XFunc(copy(config.Xcoords[1]), config)
end

function xyz_Ygen(config::ConfigFFT{T,3}) where T
    XFunc(copy(config.Xcoords[2]), config)
end

function xyz_Zgen(config::ConfigFFT{T,3}) where T
    XFunc(copy(config.Xcoords[3]), config)
end


# +++++ Coordinate in K-space ++++++++++
function k_Kgen(config::ConfigFFT{T,1}) where T
    KFunc(copy(config.Kcoords[1]), config)
end

function kl_Kgen(config::ConfigFFT{T,2}) where T
    KFunc(copy(config.Kcoords[1]), config)
end

function kl_Lgen(config::ConfigFFT{T,2}) where T
    KFunc(copy(config.Kcoords[2]), config)
end

function klm_Kgen(config::ConfigFFT{T,3}) where T
    KFunc(copy(config.Kcoords[1]), config)
end

function klm_Lgen(config::ConfigFFT{T,3}) where T
    KFunc(copy(config.Kcoords[2]), config)
end

function klm_Mgen(config::ConfigFFT{T,3}) where T
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
    @eval begin
        Base.$op(f::XFunc) = XFunc($op(f.vals), f.config)
        Base.$op(f::KFunc) = KFunc($op(f.vals), f.config)
    end
end


# *******************************************
#  Low/High-pass Filter
# *******************************************
# +++++ general pass filter +++++
function pass_K!(
            f::KFunc{T,N},
            slices::NTuple{N,UnitRange{Int}}
        ) where N where T

    vals = copy(f.vals)
    vals[slices...] .= 0.0 + 0.0im
    f.vals -= vals

end


# +++++ high-pass filter +++++
function highpass_K!(
            f::KFunc{T,N},
            min_nwaves::NTuple{N,Int}
        ) where N where T
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
            f::KFunc{T,N},
            min_nwaves::NTuple{N,Int}
        ) where N where T

    g = copy(f)
    highpass_K!(g, min_nwaves)
    return g

end


# +++++ low-pass filter +++++
function lowpass_K!(
            f::KFunc{T,N},
            max_nwaves::NTuple{N,Int}
        ) where N where T

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
            f::KFunc{T,N},
            max_nwaves::NTuple{N,Int}
        ) where N where T

    g = copy(f)
    lowpass_K!(g, max_nwaves)
    return g

end

# helper funtcion
make_slice(x::Int, y::Int) = x:y


# +++++ aliases for specific dimensins +++++
# @@@ high-pass filter @@@
function k_highpass_k(
            f::KFunc{T,1}, min_nwaves::Tuple{Int}
        ) where T
    K_highpass_K(f, min_nwaves)
end

function kl_highpass_kl(
            f::KFunc{T,2}, min_nwaves::NTuple{2, Int}
        ) where T
    K_highpass_K(f, min_nwaves)
end

function klm_highpass_klm(
            f::KFunc{T,3}, min_nwaves::NTuple{3,Int}
        ) where T
    K_highpass_K(f, min_nwaves)
end

# @@@ low-pass filter @@@
function k_lowpass_k(
            f::KFunc{T,1}, max_nwaves::Tuple{Int}
        ) where T
    K_lowpass_K(f, max_nwaves)
end

function kl_lowpass_kl(
            f::KFunc{T,2}, max_nwaves::NTuple{2,Int}
        ) where T
    K_lowpass_K(f, max_nwaves)
end

function klm_lowpass_klm(
            f::KFunc{T,3}, max_nwaves::NTuple{3,Int}
        ) where T
    K_lowpass_K(f, max_nwaves)
end


# *******************************************
#  De-aliasing
# *******************************************
# 3/2-rule ... zero-padding
# 2/3-rule ... truncating

# de-aliased product by 3/2-rule (zero padding)
function K_dealiasedprod_32_K_K(
            f::KFunc{T,N}, g::KFunc{T,N}
        ) where N where T

    if !(f.config === g.config)
        return println("ERROR")
    end

    fvals_pad = padding(f)
    gvals_pad = padding(g)
    if T == Float64
        fgvals_pad = fft(
            real(ifft(fvals_pad)) .* real(ifft(gvals_pad))
        )
    else
        fgvals_pad = fft(
            ifft(fvals_pad) .* ifft(gvals_pad)
        )
    end
    fgvals = truncate(fgvals_pad, f.config)

    return KFunc(fgvals, f.config)

end

function padding(f::KFunc)
    config = f.config

    ngrids = config.ngrids
    pad_ngrids = @. ngrids ÷ 2 * 3

    nshifts = @. pad_ngrids÷2 - ngrids÷2
    min_nwaves = @. nshifts + 1
    max_nwaves = @. nshifts + ngrids
    slices = make_slice.(min_nwaves, max_nwaves)

    vals_shift = fftshift(f.vals)
    padded = zeros(Complex{Float64}, pad_ngrids)
    padded[slices...] = vals_shift
    padded .= ifftshift(padded)

    return padded
end

function truncate(
            padded::Array{Complex{Float64},N},
            config::ConfigFFT{T,N}
        ) where N where T

    ngrids = config.ngrids
    pad_ngrids = @. ngrids ÷ 2 * 3

    nshifts = @. pad_ngrids÷2 - ngrids÷2
    min_nwaves = @. nshifts + 1
    max_nwaves = @. nshifts + ngrids
    slices = make_slice.(min_nwaves, max_nwaves)

    vals = ifftshift(fftshift(padded)[slices...])

    N_origin = prod(ngrids)
    N_padded = prod(pad_ngrids)
    vals *= N_padded / N_origin

    return vals

end

# de-aliased product by 2/3-rule (truncation)
function K_dealiasedprod_23_K_K(f::KFunc, g::KFunc)

    if !(f.config === g.config)
        return println("ERROR")
    end
    config = f.config
    ngrids = config.ngrids
    max_nwaves = @. ngrids ÷ 3
    f_trunc = K_lowpass_K(f, max_nwaves)
    g_trunc = K_lowpass_K(g, max_nwaves)
    return K_X(X_K(f_trunc) * X_K(g_trunc))

end


# *******************************************
#  Fourier Transformation
# *******************************************
function K_X(f::XFunc{T,N}) where T where N

    g = KFunc(fft(f.vals), f.config)

    if f.config.cut_zigzag_mode
        ngrids = f.config.ngrids
        max_nwaves = @. ngrids ÷ 2 - 1
        lowpass_K!(g, max_nwaves)
    end

    return g
end


function X_K(f::KFunc{T,N}) where N where T
    if T <: Real
        XFunc(real(ifft(f.vals)), f.config)
    else
        XFunc(ifft(f.vals), f.config)
    end
end

k_x(f::XFunc{T,1} where T) = K_X(f)
x_k(f::KFunc{T,1} where T) = X_K(f)

kl_xy(f::XFunc{T,2} where T) = K_X(f)
xy_kl(f::KFunc{T,2} where T) = X_K(f)

klm_xyz(f::XFunc{T,3} where T) = K_X(f)
xyz_klm(f::KFunc{T,3} where T) = X_K(f)


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

∂x_k!(k_func::KFunc{T,1} where T) = ∂Xaxis_K!(k_func, 1)

∂x_kl!(kl_func::KFunc{T,2} where T) = ∂Xaxis_K!(kl_func, 1)
∂y_kl!(kl_func::KFunc{T,2} where T) = ∂Xaxis_K!(kl_func, 2)

∂x_klm!(klm_func::KFunc{T,3} where T) = ∂Xaxis_K!(klm_func, 1)
∂y_klm!(klm_func::KFunc{T,3} where T) = ∂Xaxis_K!(klm_func, 2)
∂z_klm!(klm_func::KFunc{T,3} where T) = ∂Xaxis_K!(klm_func, 3)

# +++++ non-destructive +++++
function K_∂Xaxis_K(f::KFunc, axis::Int)

    g = copy(f)
    ∂Xaxis_K!(g, axis)
    return g

end

function X_∂Xaxis_X(f::XFunc, axis::Int)
    return X_K(K_∂Xaxis_K(K_X(f), axis))
end

function K_laplacian_K(f::KFunc{T,N} where T where N)

    return sum(
        K_∂Xaxis_K(K_∂Xaxis_K(f, axis), axis)
        for axis = 1:N
    )

end

K_Δ_K = K_laplacian_K
X_laplacian_X = X_K ∘ K_laplacian_K ∘ K_X
X_Δ_X = X_laplacian_X

# aliases
# 1-dimensional
k_∂x_k(k_func::KFunc{T,1} where T) = K_∂Xaxis_K(k_func, 1)
x_∂x_x(x_func::XFunc{T,1} where T) = X_∂Xaxis_X(x_func, 1)

# 2-dimensional
kl_∂x_kl(kl_func::KFunc{T,2} where T) = K_∂Xaxis_K(kl_func, 1)
kl_∂y_kl(kl_func::KFunc{T,2} where T) = K_∂Xaxis_K(kl_func, 2)
kl_laplacian_kl(kl_func::KFunc{T,2} where T) = K_laplacian_K(kl_func)
kl_Δ_kl = kl_laplacian_kl

xy_∂x_xy(xy_func::XFunc{T,2} where T) = X_∂Xaxis_X(xy_func, 1)
xy_∂y_xy(xy_func::XFunc{T,2} where T) = X_∂Xaxis_X(xy_func, 2)
xy_laplacian_xy(xy_func::XFunc{T,2} where T) = X_laplacian_X(xy_func)
xy_Δ_xy = xy_laplacian_xy

# 3-dimensional
klm_∂x_klm(klm_func::KFunc{T,3} where T) = K_∂Xaxis_K(klm_func, 1)
klm_∂y_klm(klm_func::KFunc{T,3} where T) = K_∂Xaxis_K(klm_func, 2)
klm_∂z_klm(klm_func::KFunc{T,3} where T) = K_∂Xaxis_K(klm_func, 3)
klm_laplacian_klm(klm_func::KFunc{T,3} where T) = K_laplacian_K(klm_func)
klm_Δ_klm = klm_laplacian_klm

xyz_∂x_xyz(xyz_func::KFunc{T,3} where T) = X_∂Xaxis_X(xyz_func, 1)
xyz_∂y_xyz(xyz_func::KFunc{T,3} where T) = X_∂Xaxis_X(xyz_func, 2)
xyz_∂z_xyz(xyz_func::KFunc{T,3} where T) = X_∂Xaxis_X(xyz_func, 3)
xyz_laplacian_xyz(xyz_func::XFunc{T,3} where T) = X_laplacian_X(xyz_func)
xyz_Δ_xyz = xyz_laplacian_xyz

# +++ tools for vector analysis +++
# 2-dimensional
function kl2_grad_kl(kl_func::KFunc{T,2})::Vector{KFunc{T,2}} where T
    return [
        kl_∂x_kl(kl_func)
        kl_∂y_kl(kl_func)
    ]
end

function kl_rot_kl2(
            kl2_func::Vector{KFunc{T,2}}
        )::KFunc{T,2} where T

    if length(kl2_func) != 2
        return println("ERROR")
    end
    return (
        kl_∂x_kl(kl2_func[2])
        - kl_∂y_kl(kl2_func[1])
    )
end

function kl_div_kl2(
            kl2_func::Vector{KFunc{T,2}}
        )::KFunc{T,2} where T

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
            klm_func::KFunc{T,3}
        )::Vector{KFunc{T,3}} where T

    return [
        klm_∂x_klm(klm_func)
        klm_∂y_klm(klm_func)
        klm_∂z_klm(klm_func)
    ]

end

function klm3_rot_klm(
            klm3_func::Vector{KFunc{T,3}}
        )::Vector{KFunc{T,3}} where T

    if length(klm3_func) != 3
        return println("ERROR")
    end
    return [
        klm_∂y_klm(klm3_func[3]) - klm_∂z_klm(klm3_func[2])
        klm_∂z_klm(klm3_func[1]) - klm_∂x_klm(klm3_func[3])
        klm_∂x_klm(klm3_func[2]) - klm_∂y_klm(klm3_func[1])
    ]

end

function klm_div_klm3(
            klm3_func::Vector{KFunc{T,3}}
        )::KFunc{T,3} where T

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

function l2inpr_X_X(f::XFunc{T,N}, g::XFunc{T,N}) where N where T

    if f.config === g.config
        return ∫(f * g)
    else
        println("ERROR")
    end

end

# aliases
integ_x(x_func::XFunc{T,1} where T) = integ_X(x_func)
integ_xy(xy_func::XFunc{T,2} where T) = integ_X(xy_func)
integ_xyz(xyz_func::XFunc{T,3} where T) = integ_X(xyz_func)

norm_x(x_func::XFunc{T,1} where T, p::Real=2) = norm_X(x_func, p)
norm_xy(xy_func::XFunc{T,2} where T, p::Real=2) = norm_X(xy_func, p)
norm_xyz(xyz_func::XFunc{T,3} where T, p::Real=2) = norm_X(xyz_func, p)

l2inpr_x_x(x_f::XFunc{T,1}, x_g::XFunc{T,1}) where T = l2inpr_X_X(x_f, x_g)
l2inpr_xy_xy(xy_f::XFunc{T,2}, xy_g::XFunc{T,2}) where T = l2inpr_X_X(xy_f, xy_g)
l2inpr_xyz_xyz(xyz_f::XFunc{T,3}, xyz_g::XFunc{T,3}) where T = l2inpr_X_X(xyz_f, xyz_g)
