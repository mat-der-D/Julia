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

function save_good_reduced(lsts::Set{Array{Int64,1}}, f::IOStream)
  for lst in lsts
    if is_one_three_power(reverse(lst))
      println(f, reverse(lst))
    end
  end
end

function save_good_try(lsts::Set{Array{Int64,1}}, f::IOStream)
  for lst in lsts
    flag_YEAH = false
    cntr = 0
    for an_lst in permutations(lst)
      if is_one_three_power(reverse(an_lst))
        println("YEAH!!!!!!!!!!")
        println(f, reverse(an_lst))
	flag_YEAH = true
	break
      end
      cntr += 1
      if cntr > 100
        break
      end
    end
    if flag_YEAH
      break
    end
  end
end

function save_minmax(lsts::Set{Array{Int64,1}}, f::IOStream)
  min_v = 10.0
  max_v = 0.0
  for lst in lsts
    min_c = Float64(sum(map(x -> x//3^reverse(lst)[x],
    	    		      1:length(lst))))
    min_v = min(min_v, min_c)
    max_c = Float64(sum(map(x -> x//3^lst[x],
    	    		      1:length(lst))))
    max_v = max(max_v, max_c)
  end
  println(f, "min=", min_v, ", max=", max_v)
end

function save_count_chance(lsts::Set{Array{Int64,1}}, f::IOStream)

  all_num = length(lsts)
  cand_num = 0
  for lst in lsts
    min_v = Float64(sum(map(x -> x//3^reverse(lst)[x],
    	    		      1:length(lst))))
    max_v = Float64(sum(map(x -> x//3^lst[x],
    	    		      1:length(lst))))
    if min_v <= 1//1 && 1//1 <= max_v
      cand_num += 1
    end
  end
  println(f, "all=", all_num, ", candidate=", cand_num)
end

function main()
  start_dim = 11
  max_dim = 27
  an_lsts = Set([[0]])

  filename = "an_"*string(start_dim)*"-"*string(max_dim)*"_reduced.txt"
  f = open(filename, "w")

  for dim = 1:max_dim

    if dim > 1
      split!(an_lsts)
    end
    if dim >= start_dim
      println("dim="*string(dim))
      println(f, "---------- dim = ",dim," ----------")
      # save_minmax(an_lsts, f)
      save_good_reduced(an_lsts, f)
    end
  end
  println(f, "========== END ==========")

  close(f)

end

function main_try()
  start_dim = 29
  max_dim = 30
  an_lsts = Set([[0]])

  filename = "an_"*string(start_dim)*"-"*string(max_dim)*"_try.txt"
  f = open(filename, "w")

  for dim = 1:max_dim

    if dim > 1
      if dim >= start_dim && dim % 4 in (1, 2)
        println("dim="*string(dim))
	println(f, "---------- dim = ",dim," ----------")
      end
      split!(an_lsts)
    end
    if dim >= start_dim && dim % 4 in (1, 2)
      println("Ready.")
      save_good_try(an_lsts, f)
    end
  end
  println(f, "========== END ==========")

  close(f)

end

function main_trytry()
  start_dim = 30
  max_dim = 30
  an_lsts = Set([[0]])

  filename = "an_"*string(start_dim)*"-"*string(max_dim)*"_trytry.txt"
  f = open(filename, "w")

  for dim = 1:max_dim

    if dim > 1
      if dim >= start_dim && dim % 4 in (1, 2)
        println("dim="*string(dim))
	println(f, "---------- dim = ",dim," ----------")
      end
      split!(an_lsts)
      an_lsts_max = 1000000
      if length(an_lsts) > an_lsts_max
        dummy_lsts = Array{Int64, 1}[]
        for _ = 1:an_lsts_max
          lst = pop!(an_lsts)
	  push!(dummy_lsts, lst)
        end
        an_lsts = Set(dummy_lsts)
      end
    end
    if dim >= start_dim && dim % 4 in (1, 2)
      println("Ready.")
      save_good_try(an_lsts, f)
    end
  end
  println(f, "========== END ==========")

  close(f)

end