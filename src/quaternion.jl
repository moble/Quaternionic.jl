# We'll need this awkward way of getting the `components` field when we set `getproperty`
components(q::AbstractQuaternion) = getfield(q, :components)

normalize(v::AbstractVector) = v ./ hypot(v...)


"""
    Quaternion{T<:Number} <: Number

Quaternionic number type with elements of type `T`.

`QuaternionF16`, `QuaternionF32` and `QuaternionF64` are aliases for `Quaternion{Float16}`,
`Quaternion{Float32}` and `Quaternion{Float64}` respectively.  See also [`Rotor`](@ref) and
[`QuatVec`](@ref).

The functions

    quaternion(w, x, y, z)
    quaternion(x, y, z)
    quaternion(w)

create a new quaternion with the given components.  The argument `w` is the scalar
component, and `x`, `y`, and `z` are the corresponding "vector" components.  If any of these
arguments is missing, it will be set to zero.  The type of the returned quaternion will be
inferred from the input arguments, or can be specified, by passing the type parameter `T` as
above.

Note that the constants [`imx`](@ref), [`imy`](@ref), and [`imz`](@ref) can also be used
like the complex `im` to create new `Quaternion` object.

# Examples

```jldoctest
julia> quaternion(1, 2, 3, 4)
1 + 2𝐢 + 3𝐣 + 4𝐤
julia> Quaternion{Float64}(1, 2, 3, 4)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> quaternion(1.0, 2.0, 3.0, 4.0)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> quaternion(2, 3, 4)
0 + 2𝐢 + 3𝐣 + 4𝐤
julia> quaternion(1)
1 + 0𝐢 + 0𝐣 + 0𝐤
```
"""
struct Quaternion{T<:Number} <: AbstractQuaternion{T}
    components::SVector{4,T}
    Quaternion{T}(a::SVector{4,T}) where {T<:Number} = new{T}(a)
    Quaternion{T}(a::AbstractVector) where {T<:Number} = new{T}(SVector{4,T}(a))
    Quaternion{T}(a::AbstractQuaternion) where {T<:Number} = new{T}(SVector{4,T}(components(a)))
    Quaternion{T}(w,x,y,z) where {T<:Number} = new{T}(SVector{4,T}(w,x,y,z))
    Quaternion{T}(x,y,z) where {T<:Number} = new{T}(SVector{4,T}(false,x,y,z))
    Quaternion{T}(w::Number) where {T<:Number} = new{T}(SVector{4,T}(w, false, false, false))
end

quaternion(a::SVector{4,T}) where {T<:Number} = Quaternion{T}(a)
#quaternion(a::AbstractVector) = Quaternion{T}(SVector{4,T}(a))  # See below
quaternion(a::AbstractQuaternion) = Quaternion{basetype(a)}(components(a))
quaternion(w,x,y,z) = (v=SVector{4}(w,x,y,z); Quaternion{eltype(v)}(v))
quaternion(x,y,z) = (v=SVector{4}(false,x,y,z); Quaternion{eltype(v)}(v))
quaternion(w::T) where {T<:Number} = Quaternion{T}(SVector{4,T}(w, false, false, false))

Quaternion(w::SVector{4,T}) where {T} = quaternion(w)
Quaternion(w::AbstractVector) = quaternion(w)
Quaternion(w::AbstractQuaternion) = quaternion(w)
Quaternion(w,x,y,z) = quaternion(w,x,y,z)
Quaternion(x,y,z) = quaternion(x,y,z)
Quaternion(w::T) where {T<:Number} = quaternion(w)


@doc raw"""
    Rotor{T<:Number} <: Number

Quaternion of unit magnitude with elements of type `T`.  These objects can be significantly
faster *and* more accurate in certain operations representing rotations.

A rotor is typically considered to be an element of the group ``\mathrm{Spin}(3) ≃
\mathrm{SU}(2)``, which can be thought of as the subgroup of quaternions with norm 1.  They
are particularly useful as representations of rotations because a rotor ``R`` acts on a
vector ``\vec{v}`` by "conjugation" as

```math
\vec{v}' = R\, \vec{v}\, R^{-1}.
```

(which can be represented in code as `R * v / R` or, more efficiently, as `R(v)`).  This
operation preserves the inner product between any two vectors conjugated in this way, and so
is a rotation.  Note that, because there are two factors of ``R`` here, the sign of ``R``
does not affect the result.  Therefore, ``\mathrm{Spin}(3)`` forms a *double* cover of the
rotation group ``\mathrm{SO}(3)``.  For this reason, it will occasionally be useful to
disregard or arbitrarily change the sign of a `Rotor` (as in [`distance`](@ref) functions) —
though this is not generally the default, and may cause problems if the input rotors change
sign when the corresponding rotations are not so different (cf. [`unflip`](@ref)).

`RotorF16`, `RotorF32` and `RotorF64` are aliases for `Rotor{Float16}`, `Rotor{Float32}` and
`Rotor{Float64}` respectively.  See also [`Quaternion`](@ref) and [`QuatVec`](@ref).

The functions

    rotor(w, x, y, z)
    rotor(w)

create a new rotor with the given components (where the components are as described in
[`Quaternion`](@ref)), automatically normalizing them on input.  Note that this
normalization step is the key difference between the `Rotor` and `rotor` functions; if you
would like to bypass normalization, you can call

    Rotor{T}(w, x, y, z)
    Rotor{T}(w)

in the same way as `rotor`, and `w, x, y, z` will be converted to type `T`.  Alternatively,
you can call

    Rotor{T}(v)

where `v<:AbstractArray` can be converted to an `SVector{4, T}`.  If you want to handle the
normalization step, you can use [`normalize`](@ref).

However, once a `Rotor` is created, its norm will often be *assumed* to be precisely 1.  So
if its true norm is significantly different, you will likely see weird results — including
vectors with very different lengths after "rotation" by a non-unit `Rotor`.

Note that simply creating a `Quaternion` that happens to have norm 1 does not make it a
`Rotor`.  However, you can pass such a `Quaternion` to the `rotor` function and get the
desired result.

# Examples

```jldoctest
julia> rotor(1, 2, 3, 4)
rotor(0.18257418583505536 + 0.3651483716701107𝐢 + 0.5477225575051661𝐣 + 0.7302967433402214𝐤)
julia> rotor(quaternion(1, 2, 3, 4))
rotor(0.18257418583505536 + 0.3651483716701107𝐢 + 0.5477225575051661𝐣 + 0.7302967433402214𝐤)
julia> Rotor{Float16}(1, 2, 3, 4)
rotor(1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤)
julia> normalize(Rotor{Float16}(1, 2, 3, 4))
rotor(0.1826 + 0.3652𝐢 + 0.548𝐣 + 0.7305𝐤)
julia> rotor(1.0)
rotor(1.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤)
```
"""
struct Rotor{T<:Number} <: AbstractQuaternion{T}
    components::SVector{4,T}
    Rotor{T}(a::SVector{4,T}) where {T<:Number} = new{T}(a)
    Rotor{T}(a::AbstractVector) where {T<:Number} = new{T}(SVector{4,T}(a))
    Rotor{T}(a::AbstractQuaternion) where {T<:Number} = new{T}(SVector{4,T}(components(a)))
    Rotor{T}(w,x,y,z) where {T<:Number} = new{T}(SVector{4,T}(w,x,y,z))
    Rotor{T}(x,y,z) where {T<:Number} = new{T}(SVector{4,T}(false,x,y,z))
    Rotor{T}(w::Number) where {T<:Number} = new{T}(SVector{4,T}(w, false, false, false))
end

function rotor(a::SVector{4,T}) where {T<:Number}
    â = normalize(a)
    Rotor{eltype(â)}(â)
end
#rotor(a::AbstractVector) = Rotor{T}(SVector{4,T}(a))  # See below
rotor(a::AbstractQuaternion) = rotor(components(a))
function rotor(w, x, y, z)
    v = normalize(SVector{4}(w, x, y, z))
    Rotor{eltype(v)}(v)
end
rotor(x,y,z) = rotor(false, x,y,z)
rotor(w::T) where {T<:Number} = Rotor{T}(SVector{4,T}(copysign(one(w), w), false, false, false))

Rotor(w::SVector{4,T}) where {T} = rotor(w)
Rotor(w::AbstractVector) = rotor(w)
Rotor(w::AbstractQuaternion) = rotor(w)
Rotor(w,x,y,z) = rotor(w,x,y,z)
Rotor(x,y,z) = rotor(x,y,z)
Rotor(w::T) where {T<:Number} = rotor(w)

"""
    QuatVec{T<:Number} <: Number

Pure-vector quaternion with elements of type `T`.  These objects can be significantly faster
*and* more accurate in certain operations than general `Quaternion`s.

`QuatVecF16`, `QuatVecF32` and `QuatVecF64` are aliases for `QuatVec{Float16}`,
`QuatVec{Float32}` and `QuatVec{Float64}` respectively.  See also [`Quaternion`](@ref) and
[`Rotor`](@ref).

The functions

    quatvec(w, x, y, z)
    quatvec(x, y, z)
    quatvec(w)

create a new `QuatVec` with the given components (where the components are as described in
[`Quaternion`](@ref)), except that the scalar argument `w` is always set to 0.

# Examples

```jldoctest
julia> quatvec(1, 2, 3, 4)
 + 2𝐢 + 3𝐣 + 4𝐤
julia> quatvec(quaternion(1, 2, 3, 4))
 + 2𝐢 + 3𝐣 + 4𝐤
julia> quatvec(2, 3, 4)
 + 2𝐢 + 3𝐣 + 4𝐤
julia> quatvec(1)
 + 0𝐢 + 0𝐣 + 0𝐤
```
"""
struct QuatVec{T<:Number} <: AbstractQuaternion{T}
    components::SVector{4,T}
    QuatVec{T}(a::SVector{4,T}) where {T<:Number} = new{T}(a)
    QuatVec{T}(a::AbstractVector) where {T<:Number} = new{T}(SVector{4,T}(a))
    QuatVec{T}(a::AbstractQuaternion) where {T<:Number} = new{T}(SVector{4,T}(components(a)))
    QuatVec{T}(w,x,y,z) where {T<:Number} = new{T}(SVector{4,T}(w,x,y,z))
    QuatVec{T}(x,y,z) where {T<:Number} = new{T}(SVector{4,T}(false,x,y,z))
    QuatVec{T}(w::Number) where {T<:Number} = new{T}(SVector{4,T}(w, false, false, false))
end

function quatvec(v::SVector{4,T}) where {T<:Number}
    v′ = SVector{4,T}(false, v[2], v[3], v[4])
    QuatVec{T}(v′)
end
#quatvec(a::AbstractVector) = QuatVec{T}(SVector{4,T}(a))  # See below
quatvec(a::AbstractQuaternion) = quatvec(components(a))
function quatvec(w, x, y, z)
    v = SVector{4}(false, x, y, z)
    QuatVec{eltype(v)}(v)
end
quatvec(x,y,z) = quatvec(false,x,y,z)
quatvec(w::T) where {T<:Number} = QuatVec{T}(SVector{4,T}(false, false, false, false))

QuatVec(w::SVector{4,T}) where {T} = quatvec(w)
QuatVec(w::AbstractVector) = quatvec(w)
QuatVec(w::AbstractQuaternion) = quatvec(w)
QuatVec(w,x,y,z) = quatvec(w,x,y,z)
QuatVec(x,y,z) = quatvec(x,y,z)
QuatVec(w::T) where {T<:Number} = quatvec(w)

# Constructor from AbstractVector
for q ∈ (:quaternion, :rotor, :quatvec)
    @eval begin
        function $q(v::AbstractVector)
            if length(v) == 4
                $q(v[begin], v[begin+1], v[begin+2], v[begin+3])
            elseif length(v) == 3
                $q(v[begin], v[begin+1], v[begin+2])
            elseif length(v) == 1
                $q(v[begin])
            else
                throw(DimensionMismatch("Quaternion must have 1, 3, or 4 inputs"))
            end
        end
    end
end

# Type constructors
(::Type{QT})(::Type{T}) where {T<:Number,QT<:AbstractQuaternion} = QT{T}
(::Type{QT})(::Type{<:AbstractQuaternion{T}}) where {T<:Number,QT<:AbstractQuaternion} = QT{T}

# Handy aliases like `ComplexF64`, etc.
const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}
const RotorF64 = Rotor{Float64}
const RotorF32 = Rotor{Float32}
const RotorF16 = Rotor{Float16}
const QuatVecF64 = QuatVec{Float64}
const QuatVecF32 = QuatVec{Float32}
const QuatVecF16 = QuatVec{Float16}

# Handy constants like `im`
"""
    imx
    𝐢

The quaternionic unit associated with rotation about the `x` axis.  Can also be entered as
Unicode bold `𝐢` (which can be input as `\\bfi<tab>`).

Note that — just as `im` is a `Complex{Bool}` — `imx` is a `QuatVec{Bool}`, and as soon as
you multiply by a scalar of any other number type (e.g., a `Float64`) it will be promoted to
a `QuatVec` of that number type, and once you *add* a scalar it will be promoted to a
`Quaternion`.

See also [`imy`](@ref) and [`imz`](@ref).

# Examples
```jldoctest
julia> imx * imx
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imx
 + 1.2𝐢 + 0.0𝐣 + 0.0𝐤
julia> 1.2 + 3.4imx
1.2 + 3.4𝐢 + 0.0𝐣 + 0.0𝐤
julia> 1.2 + 3.4𝐢
1.2 + 3.4𝐢 + 0.0𝐣 + 0.0𝐤
```
"""
const imx = QuatVec{Bool}(false, true, false, false)
const 𝐢 = imx

"""
    imy
    𝐣

The quaternionic unit associated with rotation about the `y` axis.  Can also be entered as
Unicode bold `𝐣` (which can be input as `\\bfj<tab>`).

Note that — just as `im` is a `Complex{Bool}` — `imy` is a `QuatVec{Bool}`, and as soon as
you multiply by a scalar of any other number type (e.g., a `Float64`) it will be promoted to
a `QuatVec` of that number type, and once you *add* a scalar it will be promoted to a
`Quaternion`.

See also [`imx`](@ref) and [`imz`](@ref).

# Examples
```jldoctest
julia> imy * imy
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imy
 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤
julia> 1.2 + 3.4imy
1.2 + 0.0𝐢 + 3.4𝐣 + 0.0𝐤
julia> 1.2 + 3.4𝐣
1.2 + 0.0𝐢 + 3.4𝐣 + 0.0𝐤
```
"""
const imy = QuatVec{Bool}(false, false, true, false)
const 𝐣 = imy

"""
    imz
    𝐤

The quaternionic unit associated with rotation about the `z` axis.  Can also be entered as
Unicode bold `𝐤` (which can be input as `\\bfk<tab>`).

Note that — just as `im` is a `Complex{Bool}` — `imz` is a `QuatVec{Bool}`, and as soon as
you multiply by a scalar of any other number type (e.g., a `Float64`) it will be promoted to
a `QuatVec` of that number type, and once you *add* a scalar it will be promoted to a
`Quaternion`.

See also [`imx`](@ref) and [`imy`](@ref).

# Examples
```jldoctest
julia> imz * imz
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imz
 + 0.0𝐢 + 0.0𝐣 + 1.2𝐤
julia> 1.2 + 3.4imz
1.2 + 0.0𝐢 + 0.0𝐣 + 3.4𝐤
julia> 1.2 + 3.4𝐤
1.2 + 0.0𝐢 + 0.0𝐣 + 3.4𝐤
```
"""
const imz = QuatVec{Bool}(false, false, false, true)
const 𝐤 = imz

# Essential constructors
Base.zero(::Type{QT}) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(false, false, false, false)
Base.zero(::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Base.zero(QT)
Base.zero(::Type{Rotor{T}}) where {T} = zero(Quaternion{T})

Base.one(::Type{QT}) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(true, false, false, false)
Base.one(::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Base.one(QT)
Base.one(::Type{QuatVec{T}}) where {T} = one(Quaternion{T})

# Getting pieces of quaternions
@inline function Base.getindex(q::AbstractQuaternion, i::Integer)
    @boundscheck checkbounds(components(q), i)
    components(q)[i]
end
@inline function Base.getproperty(q::AbstractQuaternion, sym::Symbol)
    @inbounds begin
        if sym === :w
            return q[1]
        elseif sym === :x
            return q[2]
        elseif sym === :y
            return q[3]
        elseif sym === :z
            return q[4]
        elseif sym === :re
            return q[1]
        elseif sym === :im
            return q[2:4]
        elseif sym === :vec
            return q[2:4]
        else # fallback to getfield
            return getfield(q, sym)
        end
    end
end
Base.@propagate_inbounds Base.getindex(q::AbstractQuaternion, I) = [q[i] for i in I]
Base.real(::Type{T}) where {T<:AbstractQuaternion} = basetype(T)
Base.real(q::AbstractQuaternion{T}) where {T<:Number} = q[1]
Base.imag(q::AbstractQuaternion{T}) where {T<:Number} = @view components(q)[2:4]
Base.vec(q::AbstractQuaternion{T}) where {T<:Number} = @view components(q)[2:4]

# Type games
wrapper(::T) where {T} = wrapper(T)
wrapper(T::UnionAll) = T
wrapper(T::Type{Q}) where {S<:Number,Q<:AbstractQuaternion{S}} = wrapper(T.name.wrapper)
wrapper(::Type{T}, ::Type{T}) where {T<:AbstractQuaternion} = wrapper(T)  # COV_EXCL_LINE

for QT1 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
    for QT2 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
        if QT1 === QT2
            @eval wrapper(::Type{<:$QT1}, ::Type{<:$QT1}) = $QT1
        else
            @eval wrapper(::Type{<:$QT1}, ::Type{<:$QT2}) = Quaternion
        end
        @eval wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:$QT2}) where {OP} = Quaternion
    end
    @eval begin
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:Number}) where {OP} = Quaternion
        wrapper(::Type{<:Number}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = Quaternion
    end
end

wrapper(::Type{<:Rotor}, ::Val{*}, ::Type{<:Rotor}) = Rotor
wrapper(::Type{<:Rotor}, ::Val{/}, ::Type{<:Rotor}) = Rotor
wrapper(::Type{<:Rotor}, ::Val{+}, ::Type{<:Rotor}) = Quaternion
wrapper(::Type{<:Rotor}, ::Val{-}, ::Type{<:Rotor}) = Quaternion
for QT ∈ (AbstractQuaternion, QuatVec)  # Quaternion is handled below
    @eval begin
        wrapper(::Type{<:Rotor}, ::Val{+}, ::Type{<:$QT}) = Quaternion
        wrapper(::Type{<:Rotor}, ::Val{-}, ::Type{<:$QT}) = Quaternion
        wrapper(::Type{<:$QT}, ::Val{+}, ::Type{<:Rotor}) = Quaternion
        wrapper(::Type{<:$QT}, ::Val{-}, ::Type{<:Rotor}) = Quaternion
    end
end

wrapper(::Type{<:Rotor}, ::Val{*}, ::Type{<:QuatVec}) = Quaternion
wrapper(::Type{<:Rotor}, ::Val{/}, ::Type{<:QuatVec}) = Quaternion
wrapper(::Type{<:QuatVec}, ::Val{*}, ::Type{<:Rotor}) = Quaternion
wrapper(::Type{<:QuatVec}, ::Val{/}, ::Type{<:Rotor}) = Quaternion

wrapper(::Type{<:QuatVec}, ::Val{+}, ::Type{<:QuatVec}) = QuatVec
wrapper(::Type{<:QuatVec}, ::Val{-}, ::Type{<:QuatVec}) = QuatVec
wrapper(::Type{<:QuatVec}, ::Val{*}, ::Type{<:QuatVec}) = Quaternion
wrapper(::Type{<:QuatVec}, ::Val{/}, ::Type{<:QuatVec}) = Quaternion

let NT = Number
    for QT ∈ (QuatVec,)
        for OP ∈ (Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = QuatVec
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = QuatVec
            end
        end
    end
    for QT ∈ (Rotor,)
        for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = Quaternion
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = Quaternion
            end
        end
    end
end
for T ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor, Number)
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = Quaternion
        if T !== Quaternion
            @eval wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = Quaternion
        end
    end
end


basetype(::AbstractQuaternion{T}) where {T} = T
basetype(::Type{<:AbstractQuaternion{T}}) where {T} = T
Base.widen(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){widen(basetype(Q))}
Base.float(::Type{Q}) where {Q<:AbstractQuaternion{<:AbstractFloat}} = Q
Base.float(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){float(basetype(Q))}
Base.float(q::AbstractQuaternion{T}) where {T<:AbstractFloat} = q
Base.float(q::AbstractQuaternion{T}) where {T} = wrapper(q){float(T)}(float(components(q)))

Base.big(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){big(basetype(Q))}
Base.big(q::AbstractQuaternion{T}) where {T<:Number} = wrapper(q){big(T)}(q)

Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Number} =
    wrapper(Q){promote_type(basetype(Q), S)}
Base.promote_rule(::Type{QuatVec{T}}, ::Type{S}) where {T<:Number,S<:Number} =
    Quaternion{promote_type(T, S)}
Base.promote_rule(::Type{QuatVec{T}}, ::Type{S}) where {T<:Number,S<:AbstractQuaternion} =
    wrapper(wrapper(QuatVec), wrapper(S)){promote_type(T, basetype(S))}
Base.promote_rule(::Type{Q1}, ::Type{Q2}) where {Q1<:AbstractQuaternion,Q2<:AbstractQuaternion} =
    wrapper(wrapper(Q1), wrapper(Q2)){promote_type(basetype(Q1), basetype(Q2))}
