using Combinatorics

function split!(lsts::Set{Array{Int64,1}})
  new_lsts = Vector{Int64}[]
  for lst in lsts
    for elem in Set(lst)
      lst_dummy = copy(lst)
      x = findfirst(isequal(elem), lst_dummy)
      if x == 1
        popfirst!(lst_dummy)
  	new_lst = Int64[elem+1, elem+1]
  	append!(new_lst, lst_dummy)
  	push!(new_lsts, new_lst)
      elseif x == length(lst)
        pop!(lst_dummy)
  	append!(lst_dummy, Int64[elem+1, elem+1])
  	push!(new_lsts, lst_dummy)
      else
        new_lst = lst[1:x-1]
  	append!(new_lst, Int64[elem+1, elem+1])
  	append!(new_lst, lst[x+1:length(lst)])
  	push!(new_lsts, new_lst)
      end
    end
  end
  empty!(lsts)
  union!(lsts, Set(new_lsts))
end

function is_one_three_power(lst::Array{Int64,1})::Bool
  return sum(map(x -> x//3^lst[x], 1:length(lst))) == 1
end

function save_good_an_lst(lsts::Set{Array{Int64,1}}, f)
  for lst in lsts
    for an_lst in Set(permutations(lst))
      if is_one_three_power(an_lst)
        println(f, reverse(an_lst))
      end
    end
  end
end

function lucky(lst::Array{Bool,1})::Bool
  return !prod(map(!,lst))
end

function is_good_lst(lst_in::Array{Int64,1})::Bool
  return lucky([ is_one_three_power(lst)
		 for lst in Set(permutations(lst_in)) ])
end

function is_good_lsts(lsts::Set{Array{Int64,1}})::Bool
  return lucky([ is_good_lst(lst) for lst in lsts])
end

function main()
  an_lsts = Set([[0]])

  max_dim = 20
  start_dim = 11

  # f = open("an_"*string(max_dim)*".txt", "w")
  f = open("an_once2.txt", "w")

  # for dim = 1:max_dim
  #   println("dim=", dim)
  #   println(f, "---------- dim = ",dim," ----------")
  #   if dim > 1
  #     split!(an_lsts)
  #   end
  #   save_good_an_lst(an_lsts, f)
  # end
  # println(f, "========== END ==========")

  for dim = 1:max_dim
    if dim > 1
      split!(an_lsts)
    end
    if dim >= start_dim
      println("dim=", dim)
      println(f, "---------- dim = ",dim," ----------")
      println(f, is_good_lsts(an_lsts))
    end
  end
  println(f, "========== END ==========")

  close(f)
end

main()