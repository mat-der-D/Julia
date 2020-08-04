struct simple

  x::Float64
  s::String

end

Base.:+(f::simple, g::simple) = simple(f.x + g.x, f.s*g.s)
Base.:sin(f::simple) = simple(sin(f.x), "sin"*f.s)

a = simple(1., "one")
print(sin(a))