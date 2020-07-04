function sum_by_view(lst::Array{Float64,1})::Float64
  return sum(view(lst, :))
end

function sum_by_slice(lst::Array{Float64,1})::Float64
  return sum(lst[:])
end

function sum_by_for(lst::Array{Float64,1})::Float64
  sumval = Float64(0)
  for x in lst
    sumval += x
  end
  return sumval
end