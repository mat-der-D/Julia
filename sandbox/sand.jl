using LinearAlgebra

mutable struct MyArray <: AbstractArray{Int, 2}

        vals::Matrix{Int}

end
Base.size(A::MyArray) = size(A.vals)
Base.getindex(A::MyArray, i::Int) = A.vals[i]
Base.getindex(A::MyArray, I::Vararg{Int, 2}) = A.vals[I...]
Base.setindex!(A::MyArray, v, i::Int) = setindex!(A, v, i)
Base.setindex!(A::MyArray, v, I::Vararg{Int, 2}) = setindex!(A, v, I...)
a = MyArray([1 2; 3 4])
b = MyArray([5 6; 7 8])

BINOP = ((:*, :.*), (:+, :.+), (:^, :.^))

for (op, opd) = BINOP
    @eval begin
        Base.$op(A::MyArray, B::MyArray) = (
            MyArray($opd(A.vals, B.vals))
        )
        Base.$op(A::MyArray, a::Real) = (
            MyArray($opd(A.vals, a))
        )
        Base.$op(a::Real, B::MyArray) = (
            MyArray($opd(a, B.vals))
        )
    end
end
