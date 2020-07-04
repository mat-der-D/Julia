using Combinatorics

function split!(lsts::Set{Array{Int16,1}})
  new_lsts = Vector{Int16}[]
  for lst in lsts
    for elem in Set(lst)
      lst_dummy = copy(lst)
      x = findfirst(isequal(elem), lst_dummy)
      if x == 1
        popfirst!(lst_dummy)
  	new_lst = Int16[elem+1, elem+1]
  	append!(new_lst, lst_dummy)
  	push!(new_lsts, new_lst)
      elseif x == length(lst)
        pop!(lst_dummy)
  	append!(lst_dummy, Int16[elem+1, elem+1])
  	push!(new_lsts, lst_dummy)
      else
        new_lst = lst[1:x-1]
  	append!(new_lst, Int16[elem+1, elem+1])
  	append!(new_lst, lst[x+1:length(lst)])
  	push!(new_lsts, new_lst)
      end
    end
  end
  empty!(lsts)
  union!(lsts, Set(new_lsts))
end

function num_perm(lst::Array{Int16,1})::UInt16
  return length(Set(permutations(lst)))
end

function main()

  an_lsts = Set([[0]])
  println(typeof(an_lsts))
  split!(an_lsts)
  # f = open("perm_num.txt", "w")

  # for dim = 1:2
  #   if dim > 1
  #     split!(an_lsts)
  #   end
  #   # println(f, dim, " ", sum([num_perm.(lst) for lst in an_lsts]))
  # end

  # close(f)

end

main()