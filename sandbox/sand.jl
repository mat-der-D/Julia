using IterTools

n = 4
I = collect(1:n)

for ss in subsets(I)
  length(ss) == 0 && continue
  @show ss
end