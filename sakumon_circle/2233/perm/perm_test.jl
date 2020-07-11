using Combinatorics

function perm(lst::Array{Int64,1};key=1)

  if key == 1
    return perm_1(lst)
  elseif key == 2
    return perm_2(lst)
  elseif key == 3

  end

end

function perm_1(lst::Array{Int64,1})::Set{Array{Int64,1}}

  return Set(permutations(lst))

end

function perm_2(lst::Array{Int64,1})::Set{Array{Int64,1}}

  all_perm = permutations(lst)
  seen = Set(Array{Int64,1}[])
  for x in all_perm
    push!(seen, x)
  end
  return seen

end

function test(key::Int64)

  lst = [1, 2, 2, 5, 5, 5]
  perm(lst, key = key)

end