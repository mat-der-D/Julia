# fft version

using FFTW

# *******************************************
#  Configuration
# *******************************************
struct ConfigFFT{N}

    ngrids::Tuple{Vararg{Int, N}} # Number Of Grids
    nwaves::Tuple{Vararg{Int, N}} # Cutoff wavenumber,
    		                  # ngrids[_} >= nwaves[_]
    xmins::Tuple{Vararg{Float64, N}} # Min values of space
    xmaxs::Tuple{Vararg{Float64, N}} # Max values of space

    # constructor
    function ConfigFFT{N}(ngrids::Tuple{Vararg{Int, N}},
                          nwaves::Tuple{Vararg{Int, N}},
			  xmins::Tuple{Vararg{Float64, N}},
			  xmaxs::Tuple{Vararg{Float64, N}}
			  ) where N
        if all(ngrids .>= nwaves)
	    return new(ngrids, nwaves, xmins, xmaxs)
	else
	    println("ERROR")
	end
    end

end


function configure_FFT(; ngrids::Array{Int, 1},
	 		 nwaves::Array{Int, 1},
			 xmins::Array{Float64, 1},
			 xmaxs::Array{Float64, 1})

    N = length(ngrids)
    to_tup(x...) = x
    ngrids_tup = to_tup(ngrids...)
    nwaves_tup = to_tup(nwaves...)
    xmins_tup = to_tup(xmins...)
    xmaxs_tup = to_tup(xmaxs...)
    return ConfigFFT{N}(ngrids_tup, nwaves_tup,
                        xmins_tup, xmaxs_tup)

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
            return new(copy(vals), config)
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
Base.:copy(f::XFunc{N}) where N = XFunc{N}(f.vals, f.config)

# --- operators ---

# To Be Implemented!


# +++++ KFunc +++++++++++++++++++++++++
mutable struct KFunc{N} <: AbstractArray{Complex{Float64}, N}

    vals::Array{Complex{Float64}, N}
    config::ConfigFFT{N}

    # constructor
    function KFunc{N}(vals::Array{Complex{Float64}, N},
                      config::ConfigFFT{N}) where N
        if size(vals) == config.nwaves
            return new(copy(vals), config)
        else
            println("ERROR")
        end
    end

    function KFunc(vals::Array{T, N},
    	           config::ConfigFFT{N}) where N where T <: Number
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

# --- operators ---

# To Be Implemented!

# *******************************************
#  Basic Tools
# *******************************************
# +++++ Coordinate in X-space ++++++++++
function Xcoord(index::Int, config::ConfigFFT{N}) where N
    ngrid = config.ngrids[index]
    xmin = config.xmins[index]
    xmax = config.xmaxs[index]
    dx = (xmax - xmin) / ngrid
    xvalue(inds) = xmin + (inds[index] - 1)*dx

    UNDEF = Array{Any, N}(undef, config.ngrids)
    coord = xvalue.( CartesianIndices(UNDEF) )
    return XFunc(coord, config)
end

function x_Xgen(config::ConfigFFT{1})
    return Xcoord(1, config)
end

function xy_Xgen(config::ConfigFFT{2})
    return Xcoord(1, config)
end

function xy_Ygen(config::ConfigFFT{2})
    return Xcoord(2, config)
end

function xyz_Xgen(config::ConfigFFT{3})
    return Xcoord(1, config)
end

function xyz_Ygen(config::ConfigFFT{3})
    return Xcoord(2, config)
end

function xyz_Zgen(config::ConfigFFT{3})
    return Xcoord(3, config)
end

# +++++ Coordinate in K-space ++++++++++
function Kcoord(index::Int, config::ConfigFFT{N}) where N
    ngrid = config.ngrids[index]
    half_ngrid = div(ngrid, 2)
    
    function kvalue(inds)
        ind = inds[index]
	if ind <= half_ngrid
	    return ind - 1
	else
	    return ind - ngrids - 1
	end
    end

    UNDEF = Array{Any, N}(undef, config.ngrids)
    coord = kvalue.( CartesianIndices(UNDEF) )
    return KFunc(coord, config)
end

function k_Kgen(config::ConfigFFT{1})
    return Kcorrd(1, config)
end

function kl_Kgen(config::ConfigFFT{2})
    return Kcorrd(1, config)
end

function kl_Lgen(config::ConfigFFT{2})
    return Kcorrd(2, config)
end

function klm_Kgen(config::ConfigFFT{3})
    return Kcorrd(1, config)
end

function klm_Lgen(config::ConfigFFT{3})
    return Kcorrd(2, config)
end

function klm_Mgen(config::ConfigFFT{3})
    return Kcorrd(3, config)
end


# *******************************************
#  Fourier Transformation
# *******************************************
K_X(f::XFunc) = KFunc(fft(f.vals), f.config)
X_K(f::KFunc) = XFunc(ifft(f.vals), f.config)

k_x(f::XFunc{1}) = K_X(f)
x_k(f::KFunc{1}) = X_K(f)

kl_xy(f::XFunc{2}) = K_X(f)
xy_kl(f::KFunc{2}) = X_K(f)

klm_xyz(f::XFunc{3}) = K_X(f)
xyz_klm(f::KFunc{3}) = X_K(f)