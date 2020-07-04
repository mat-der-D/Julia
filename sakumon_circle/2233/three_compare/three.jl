function is_one_three_power(lst::Array{Int64,1}; type=1)::Bool

  if type == 1
    return is_one_three_power_1(lst)
  elseif type == 2
    return is_one_three_power_2(lst)
  end

end

function is_one_three_power_1(lst::Array{Int64,1})::Bool
  local sumval::Rational{Int64} = 0
  local N::Int64 = length(lst)
  for i = 1:N
    sumval += (N - i + 1)//(3^lst[i])
    if sumval > 1
      return false
    end
  end
  return sumval == 1
end

function is_one_three_power_2(lst::Array{Int64,1})::Bool
  return sum(map(x -> x//3^lst[x], 1:length(lst))) == 1
end

function test(type::Int64)
  an = [1,1,5,1,14,1,1,43,1,24,56]
  b = is_one_three_power(an, type=1)
end