function split!(lsts::Set{Array{Int64,1}}; type::Int64=1)
  if type == 3
    split_fast!(lsts)
  elseif type == 4
    split_faster!(lsts)
  elseif type == 5
    split_view!(lsts)
  else
    new_lsts = [split_lst(lst, elem, type = type)
	        for lst in lsts
	        for elem in Set(lst)]
    empty!(lsts)
    union!(lsts, Set(new_lsts))
  end
end

function split_fast!(lsts::Set{Array{Int64,1}})
  N_lst = length( pop!(copy(lsts)) )
  new_lst = Array{Int64,1}(undef, N_lst+1)
  new_lsts = Array{Int64,1}[]

  for lst in lsts
    for elem in Set(lst)
      x = findfirst(isequal(elem), lst)
      if x != 1
        new_lst[1:x-1] = lst[1:x-1]
      end
      new_lst[x:x+1] = Int64[elem+1, elem+1]
      if x != N_lst
      	new_lst[x+2:N_lst+1] = lst[x+1:N_lst]
      end
      push!(new_lsts, new_lst[:])
    end
  end
  
  empty!(lsts)
  union!(lsts, Set(new_lsts))
end

function split_view!(lsts::Set{Array{Int64,1}})
  N_lst = length( pop!(copy(lsts)) )
  new_lst = Array{Int64,1}(undef, N_lst+1)
  new_lsts = Array{Int64,1}[]

  for lst in lsts
    for elem in Set(lst)
      x = findfirst(isequal(elem), lst)
      if x != 1
        new_lst[1:x-1] = view(lst, 1:x-1)
      end
      new_lst[x:x+1] = Int64[elem+1, elem+1]
      if x != N_lst
      	new_lst[x+2:N_lst+1] = view(lst, x+1:N_lst)
      end
      push!(new_lsts, new_lst[:])
    end
  end
  
  empty!(lsts)
  union!(lsts, Set(new_lsts))
end

function split_faster!(lsts::Set{Array{Int64,1}})
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



function split_lst(lst::Array{Int64,1}, val::Int64;
	 	type::Int64=1)::Array{Int64,1}

  if !(val in lst)
    throw(ErrorException("'lst' must include 'val'. [split_lst]"))
  end

  if type == 1
    return split_lst_1(lst, val)
  elseif type == 2
    return split_lst_2(lst, val)
  else
    throw(ErrorException("'type' is invalid [split_lst]"))
  end

end

function split_lst_1(lst::Array{Int64,1},
		val::Int64)::Array{Int64,1}
  x = findfirst(isequal(val), lst)
  if x == 1
    lst_dummy = copy(lst)
    popfirst!(lst_dummy)
    new_lst = Int64[val+1, val+1]
    append!(new_lst, lst_dummy)
    return new_lst
  elseif x == length(lst)
    lst_dummy = copy(lst)
    pop!(lst_dummy)
    append!(lst_dummy, Int64[val+1, val+1])
    return lst_dummy
  else
    new_lst = lst[1:x-1]
    append!(new_lst, Int64[val+1, val+1])
    append!(new_lst, lst[x+1:length(lst)])
    return new_lst
  end
end

function split_lst_2(lst::Array{Int64, 1},
	 	val::Int64)::Array{Int64,1}
  N_lst = length(lst)
  new_lst = Array{Int64, 1}(undef, N_lst+1)
  x = findfirst(isequal(val), lst)

  if x != 1
    new_lst[1:x-1] = lst[1:x-1]
  end
  new_lst[x:x+1] = Int64[val+1, val+1]
  if x != N_lst
    new_lst[x+2:N_lst+1] =lst[x+1:N_lst]
  end

  return new_lst

end

function test(type::Int64)
  an = Set([[0]])
  for i = 1:18
    split!(an, type=type)
  end
end


function main()


end
