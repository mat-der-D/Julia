function subs_by_view(lst::Array{Float64,1})::Array{Float64,1}
  new_lst = Array{Float64, 1}(undef, length(lst))
  new_lst = view(lst, 2:(length(lst) - 1))
  return new_lst
end

function subs_by_slice(lst::Array{Float64,1})::Array{Float64,1}
  new_lst = Array{Float64, 1}(undef, length(lst))
  new_lst = lst[2:end-1]
  return new_lst
end

function subs_by_for(lst::Array{Float64,1})::Array{Float64,1}
  new_lst = Array{Float64, 1}(undef, length(lst))
  for n = 2:length(lst)-1
    new_lst[n] = lst[n]
  end
  return new_lst
end