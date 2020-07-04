fuction split_lst!(lst::Array{Int64,1}, val::Int64, type::Int64=0)
  if !(val in lst)
    throw("[splt_lst!] 'lst' must include 'val'.")
  end
 
  if type == 0
    x = findfirwt(isequal(val), lst_dummy)
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
  elseif type == 1

  else
    throw("[split_lst!] Exception in 'type'."
  end
end