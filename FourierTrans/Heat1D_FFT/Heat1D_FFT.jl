include("FFT1D_derivative.jl")

struct heat1D_setup

    D::Float64		# diffusion coefficient
    dt::Float64		# time step
    nt::Int64		# output each nt step

end


function k_time_evolve_heat1D_k(k_u0, 
	 		        setup::heat1D_setup,
	 		        config::FFT1D_config)
    return exp.(- setup.dt*setup.nt*
    	   	  setup.D*(config.k_X .^ 2)).* k_u0
end

function output_x(f::IOStream, x_func, config::FFT1D_config)
    for ix in 1:config.nx
        println(f, config.x_X[ix], " ", x_func[ix])
    end
end

function test_heat1D()
    c = init_FFT1D(128, 0., 1.)
    s = heat1D_setup(0.01, 1e-4, 100)

    x_u0 = zeros(size(c.x_X))
    x_u0[div(c.nx, 2)+1] = 1.

    open("testcase.txt", "w") do f
        output_x(f, x_u0, c)
    	x_u = x_u0
    	for _ = 1:100
            k_u = k_x(x_u)
	    k_u = k_time_evolve_heat1D_k(k_u, s, c)
	    x_u = x_k(k_u)
	    println(f, "\n")
	    output_x(f, x_u, c)
    	end
    end
end