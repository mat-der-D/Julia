include("FFT1D_derivative.jl")

struct setup_Barkley

    zeta::Float64
    D::Float64
    delta::Float64
    eps_1::Float64
    eps_2::Float64
    U_0::Float64
    U_bar::Float64
    r::Float64
    dt::Float64
    nt::Int64

end

function init_setup_Barkley(;zeta, D, delta, eps_1, eps_2,
	 		     U_0, U_bar, r, dt, nt)
    return setup_Barkley(zeta, D, delta, eps_1, eps_2,
    	   		 U_0, U_bar, r, dt, nt)
end

function x_x_time_evolve_Barkley_x_x(
		x_q, x_u,
		s::setup_Barkley,
		c::FFT1D_config)

    for _ = 1:s.nt
      x_q_mid, x_u_mid =
        (x_q, x_u) .+ x_x_nonlin_terms_x_x(x_q, x_u, s, c) .* s.dt
      k_q_mid, k_u_mid = k_x(x_q_mid), k_x(x_u_mid)
      k_q, k_u = k_k_linear_evolve_k_k(k_q_mid, k_u_mid, s, c)
      x_q, x_u = x_k(k_q), x_k(k_u)
    end

    return x_q, x_u

end

function k_k_linear_evolve_k_k(
		k_q, k_u,
		s::setup_Barkley,
		c::FFT1D_config)

    A = @. - (s.delta + s.U_0) - s.D * c.k_X ^ 2 +
      	     (im * s.zeta) * c.k_X
    B = ( s.eps_2 * s.U_bar ) .* ones(size(c.k_X))
    C = - s.eps_1 .* ones(size(c.k_X))
    k_qNew = @. exp(s.dt*A) * k_q
    k_uNew = @. B *
    	     	(exp(s.dt*A) - exp(s.dt*C)) / (A - C) *
		k_q + exp(s.dt*C) * k_u

    return k_qNew, k_uNew
end

function x_x_nonlin_terms_x_x(
		x_q, x_u,
		s::setup_Barkley,
		c::FFT1D_config)

    x_qFunc = x_q
    x_qFunc = x_q.*x_u .-
    	      (s.r + s.delta).*x_q.^2 .*(x_q .- 2.) .-
	      x_u .* x_DX_x(x_q, c)
    x_uFunc = s.eps_1 * s.U_0 .- x_u .* x_DX_x(x_u, c) .-
    	      s.eps_2 .* x_u .* x_q

    return x_qFunc, x_uFunc
end

function output_x(f::IOStream, x_func, c::FFT1D_config)
    for ix in 1:c.nx
        println(f, c.x_X[ix], " ", x_func[ix])
    end
end

function output_x_x(f::IOStream, x_func1, x_func2, c::FFT1D_config)
    for ix in 1:c.nx
        println(f, c.x_X[ix], " ", x_func1[ix], " ", x_func2[ix])
    end
end

function test_Barkley()

    c = init_FFT1D(512, 0., 30.)
    s = init_setup_Barkley(
      	zeta=0.8, D=0.5, delta=0.1, eps_1=0.1, eps_2=0.2,
	U_0=2., U_bar=1., r=10.0, dt=1e-4, nt=1000)

    x_q = @. exp(sin(16pi*c.x_X/c.xr))
    x_u = @. exp(cos(16pi*c.x_X/c.xr))

    open("testcase.txt", "w") do f
        output_x_x(f, x_q, x_u, c)
	for _ = 1:100
	    x_q, x_u =
	    	x_x_time_evolve_Barkley_x_x(x_q, x_u, s, c)
	    println(f, "\n")
	    output_x_x(f, x_q, x_u, c)
	end
    end

end