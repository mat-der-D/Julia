function RK4_solve(f,
	 t0::Float64, x0::Float64,
	 N::Int64, dt::Float64)::Tuple{Float64, Float64}
  t, x = t0, x0
  for _ in 1:N
    k1 = f(t       , x          )
    k2 = f(t + dt/2, x + k1/2*dt)
    k3 = f(t + dt/2, x + k2/2*dt)
    k4 = f(t + dt  , x + k3  *dt)
    t += dt
    x += dt * (k1 + 2k2 + 2k3 + k4)/6
  end
  return t, x
end


function RK4_solve(f,
	 t0::Float64, x0::Array{Float64, 1},
	 N::Int64, dt::Float64)::Tuple{Float64, Array{Float64, 1}}

  t, x = t0, x0
  for _ in 1:N
    k1 = f(t       , x          )
    k2 = f(t + dt/2, x + k1/2*dt)
    k3 = f(t + dt/2, x + k2/2*dt)
    k4 = f(t + dt  , x + k3  *dt)
    t += dt
    x += dt * (k1 + 2k2 + 2k3 + k4)/6
  end
  return t, x
end

function RK4_solve_dot(f,
	 t0::Float64, x0::Array{Float64, 1},
	 N::Int64, dt::Float64)::Tuple{Float64, Array{Float64, 1}}

  t, x = t0, x0
  for _ in 1:N
    k1 = f(t       , x             )
    k2 = f(t + dt/2, @. x + k1/2*dt)
    k3 = f(t + dt/2, @. x + k2/2*dt)
    k4 = f(t + dt  , @. x + k3  *dt)
    t += dt
    @. x += dt * (k1 + 2k2 + 2k3 + k4)/6
  end
  return t, x
end

function sin_solve()

  # Parameters
  F(t, x) = [x[2], -x[1]]
  t0 = 0.0
  x0 = [0.0, 1.0]
  dt = 0.001
  N  = Int(div(3pi, dt))
  N_out = div(N, 200)

  f = open("Sin.txt", "w")

  # initial condition
  t, x = t0, x0
  println(f, t, " ", x[1])
  for _ = 1:div(N, N_out)
    t, x = RK4_solve(F, t, x, N_out, dt)
    println(f, t, " ", x[1])
  end
  close(f)

end

function sin_solve_dot()

  # Parameters
  F(t, x) = [x[2], -x[1]]
  t0 = 0.0
  x0 = [0.0, 1.0]
  dt = 0.001
  N  = Int(div(3pi, dt))
  N_out = div(N, 200)

  f = open("Sin.txt", "w")

  # initial condition
  t, x = t0, x0
  println(f, t, " ", x[1])
  for _ = 1:div(N, N_out)
    t, x = RK4_solve_dot(F, t, x, N_out, dt)
    println(f, t, " ", x[1])
  end
  close(f)

end

function exp_solve()

  # Parameters
  F(t, x) = x
  t0 = 0.0
  x0 = 1.0
  dt = 0.001
  N  = 2000
  N_out = 10

  f = open("Exp.txt", "w")

  # initial condition
  t, x = t0, x0
  println(f, t, " ", x)
  for _ = 1:div(N, N_out)
    t, x = RK4_solve(F, t, x, N_out, dt)
    println(f, t, " ", x)
  end
  close(f)

end