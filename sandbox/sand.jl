function my_func(x::NTuple{N,Int}) where N
    x[1]
end

OPNAMES = ((:my1, 1), (:my2, 2))

for (op, dim) = OPNAMES
    @eval begin
        $op(x::NTuple{$dim,Int}) = my_func(x)
    end
end
