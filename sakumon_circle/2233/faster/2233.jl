using Combinatorics

function split!(lsts::Set{Array{Int64,1}})
  N_lst = length( pop!(copy(lsts)) )
  new_lst = Array{Int64,1}(undef, N_lst+1)
  new_lsts = Array{Int64,1}[]

  for lst in lsts
    for elem in Set(lst)
      x = findfirst(isequal(elem), lst)
      for n in 1:N_lst+1
        if n <= x-1
	  new_lst[n] = lst[n]
	elseif n >= x+2
	  new_lst[n] = lst[n-1]
	else
	  new_lst[n] = elem+1
	end
      end
      push!(new_lsts, new_lst[:])
    end
  end
  
  empty!(lsts)
  union!(lsts, Set(new_lsts))
end

function f(x::Int64, lst::Array{Int64,1})::Rational{Int64}
  return x//3^lst[x]
end

@inline function is_one_three_power(lst::Array{Int64,1})::Bool
  # return sum(map(x::Int64 -> x//3^lst[x], 1:length(lst))) == 1//1
  return sum(map(x -> x//3^lst[x], 1:length(lst))) == 1//1
end

function alert_good_an_lst(an_lsts::Array{Array{Int64,1},1})
  for lst in an_lsts
    good_an_lsts::Array{Array{Int64,1},1} = []
    done_an_lsts::Array{Array{Int64,1},1} = []
    for x in permutations(lst)
      if !(x in done_an_lsts)
        push!(done_an_lsts, x)
    	if is_one_three_power(x)
    	  print(reverse(x), "\n")
    	  push!(good_an_lsts, x)
    	end
      end
    end
  end
end

function save_good_an_lst(lsts::Set{Array{Int64,1}}, f::IOStream)
  for lst in lsts
    for an_lst in Set(permutations(lst))
      if is_one_three_power(an_lst)
        println(f, an_lst)
      end
    end
  end
end

function main()
  start_dim = 1
  max_dim = 10
  an_lsts = Set([[0]])

  filename = "an_"*string(start_dim)*"-"*string(max_dim)*".txt"
  f = open(filename, "w")

  for dim = 1:max_dim

    if dim > 1
      split!(an_lsts)
    end
    if dim >= start_dim
      println("dim="*string(dim))
      println(f, "---------- dim = ",dim," ----------")
      save_good_an_lst(an_lsts, f)
    end
  end
  println(f, "========== END ==========")

  close(f)

end

main()