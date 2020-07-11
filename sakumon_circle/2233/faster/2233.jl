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

@inline function is_one_three_power(lst::Array{Int64,1})::Bool
  return sum(map(x -> x//3^lst[x], 1:length(lst))) == 1//1
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

function save_good_an_lst_perm(lsts::Set{Array{Int64,1}}, f::IOStream)
  for lst in lsts
    for an_lst in permutations(lst)
      if is_one_three_power(an_lst)
        println(f, an_lst)
      end
    end
  end
end

function save_good_an_lst_exit(lsts::Set{Array{Int64,1}}, f::IOStream)
  for lst in lsts
    for an_lst in Set(permutations(lst))
      if is_one_three_power(an_lst)
        println(f, an_lst)
	quit()
      end
    end
  end
end

function main()
  start_dim = 9
  max_dim = 9
  an_lsts = Set([[0]])

  filename = "an_"*string(start_dim)*"-"*string(max_dim)*"_exit.txt"
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

# main()