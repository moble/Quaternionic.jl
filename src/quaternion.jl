"""
    Quaternion{T<:Number} <: Number

Quaternionic number type with elements of type `T`.

`QuaternionF16`, `QuaternionF32` and `QuaternionF64` are aliases for `Quaternion{Float16}`,
`Quaternion{Float32}` and `Quaternion{Float64}` respectively.  See also [`Rotor`](@ref) and
[`QuatVec`](@ref).

The functions

    Quaternion(w, x, y, z)
    Quaternion(x, y, z)
    Quaternion(w)
    Quaternion{T}(w, x, y, z)

create a new quaternion with the given components.  The argument `w` is the
scalar component, and `x`, `y`, and `z` are the corresponding "vector"
components.  If any of these arguments is missing, it will be set to zero.  The
type of the returned quaternion will be inferred from the input arguments, or
can be specified, by passing the type parameter `T` as above.

Note that the constants [`imx`](@ref), [`imy`](@ref), and [`imz`](@ref) can
also be used like the complex `im` to create new `Quaternion` object.

# Examples

```jldoctest
julia> Quaternion(1, 2, 3, 4)
1 + 2𝐢 + 3𝐣 + 4𝐤
julia> Quaternion{Float64}(1, 2, 3, 4)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> Quaternion(1.0, 2.0, 3.0, 4.0)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> Quaternion(2, 3, 4)
0 + 2𝐢 + 3𝐣 + 4𝐤
julia> Quaternion(1)
1 + 0𝐢 + 0𝐣 + 0𝐤
```
"""
struct Quaternion{T<:Number} <: AbstractQuaternion{T}
    components::SVector{4,T}
    Quaternion{T}(a::SVector{4,T}) where {T<:Number} = new{T}(a)
    Quaternion{T}(a::A) where {T<:Number,A<:AbstractVector} = new{T}(SVector{4,T}(a))
end

@doc raw"""
    Rotor{T<:Number} <: Number

Quaternion of unit magnitude with elements of type `T`.  These objects can be
significantly faster *and* more accurate in certain operations representing
rotations.

A rotor is typically considered to be an element of the group
``\mathrm{Spin}(3) ≃ \mathrm{SU}(2)``, which can be thought of as the subgroup
of quaternions with norm 1.  They are particularly useful as representations of
rotations because a rotor ``R`` acts on a vector ``\vec{v}`` by "conjugation" as

```math
\vec{v}' = R\, \vec{v}\, R^{-1}.
```

This preserves the inner product between any two vectors conjugated in this
way, and so is a rotation.  Note that, because there are two factors of ``R``
here, the sign of ``R`` does not affect the result.  Therefore,
``\mathrm{Spin}(3)`` forms a *double* cover of the rotation group
``\mathrm{SO}(3)``.  For this reason, it will occasionally be useful to
disregard or arbitrarily change the sign of a `Rotor` (as in [`distance`](@ref)
functions) — though this is not generally the default, and may cause problems
if the input rotors change sign when the corresponding rotations are not so
different (cf. [`unflip`](@ref)).

`RotorF16`, `RotorF32` and `RotorF64` are aliases for `Rotor{Float16}`,
`Rotor{Float32}` and `Rotor{Float64}` respectively.  See also
[`Quaternion`](@ref) and [`QuatVec`](@ref).

The functions

    Rotor(w, x, y, z)
    Rotor(x, y, z)
    Rotor(w)
    Rotor{T}(w, x, y, z)

create a new rotor with the given components (where the components are as
described in [`Quaternion`](@ref)), automatically normalizing them on input.
If you would like to bypass this normalization step, you can call

    Rotor{T}(v)

where `v<:AbstractArray`, and can be converted to an `SVector{4, T}`.

However, once a `Rotor` is created, its norm will always be assumed to be
precisely 1.  So if its true norm is significantly different, you will like see
weird results — including vectors with very different lengths after "rotation"
by a non-unit `Rotor`.

Note that simply creating a `Quaternion` that happens to have norm 1 does not
make it a `Rotor`.  However, you can pass such a `Quaternion` to the `Rotor`
function and get the desired result.

# Examples

```jldoctest
julia> Rotor(1, 2, 3, 4)
Rotor(0.18257418583505536 + 0.3651483716701107𝐢 + 0.5477225575051661𝐣 + 0.7302967433402214𝐤)
julia> Rotor(Quaternion(1, 2, 3, 4))
Rotor(0.18257418583505536 + 0.3651483716701107𝐢 + 0.5477225575051661𝐣 + 0.7302967433402214𝐤)
julia> Rotor{Float16}(1, 2, 3, 4)
Rotor(0.1826 + 0.3652𝐢 + 0.548𝐣 + 0.7305𝐤)
julia> Rotor(2, 3, 4)
Rotor(0.0 + 0.3713906763541037𝐢 + 0.5570860145311556𝐣 + 0.7427813527082074𝐤)
julia> Rotor(1)
Rotor(1 + 0𝐢 + 0𝐣 + 0𝐤)
```
"""
struct Rotor{T<:Number} <: AbstractQuaternion{T}
    components::SVector{4,T}
    Rotor{T}(a::SVector{4,T}) where {T<:Number} = new{T}(a)
    Rotor{T}(a::A) where {T<:Number,A<:AbstractVector} = new{T}(SVector{4,T}(a))
end

"""
    QuatVec{T<:Number} <: Number

Pure-vector quaternion with elements of type `T`.  These objects can be
significantly faster *and* more accurate in certain operations than general
`Quaternion`s.

`QuatVecF16`, `QuatVecF32` and `QuatVecF64` are aliases for `QuatVec{Float16}`,
`QuatVec{Float32}` and `QuatVec{Float64}` respectively.  See also
[`Quaternion`](@ref) and [`Rotor`](@ref).

The functions

    QuatVec(w, x, y, z)
    QuatVec(x, y, z)
    QuatVec(w)
    QuatVec{T}(w, x, y, z)

create a new rotor with the given components (where the components are as
described in [`Quaternion`](@ref)), except that the scalar argument `w` is
always set to 0.

# Examples

```jldoctest
julia> QuatVec(1, 2, 3, 4)
 + 2𝐢 + 3𝐣 + 4𝐤
julia> QuatVec(Quaternion(1, 2, 3, 4))
 + 2𝐢 + 3𝐣 + 4𝐤
julia> QuatVec(2, 3, 4)
 + 2𝐢 + 3𝐣 + 4𝐤
julia> QuatVec(1)
 + 0𝐢 + 0𝐣 + 0𝐤
```
"""
struct QuatVec{T<:Number} <: AbstractQuaternion{T}
    components::SVector{4,T}
    QuatVec{T}(a::SVector{4,T}) where {T<:Number} = new{T}(a)
    QuatVec{T}(a::A) where {T<:Number,A<:AbstractVector} = new{T}(SVector{4,T}(a))
end

# We'll need this awkward way of getting the `components` field when we set `getproperty`
components(q::AbstractQuaternion) = getfield(q, :components)

# Untyped constructor from SVector
# (::Type{QT})(v::SVector{4, T}) where {T<:Number, QT<:AbstractQuaternion} = QT{eltype(v)}(v)
Quaternion(v::SVector{4,T}) where {T<:Number} = Quaternion{eltype(v)}(v)
Rotor(v::SVector{4,T}) where {T<:Number} = Rotor{eltype(v)}(v)
QuatVec(v::SVector{4,T}) where {T<:Number} = QuatVec{eltype(v)}(v)

# Constructor from AbstractVector
(::Type{QT})(v::AbstractVector{T}) where {QT<:AbstractQuaternion, T} = QT{T}(v...)

# Constructor from all 4 components
(::Type{QT})(w, x, y, z) where {QT<:AbstractQuaternion} = (v = SVector{4}(w, x, y, z); QT{eltype(v)}(v))
(::Type{QT})(w, x, y, z) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(SVector{4,T}(w, x, y, z))
Rotor(w, x, y, z) = (n = √(w^2 + x^2 + y^2 + z^2); Rotor(SVector{4}(w / n, x / n, y / n, z / n)))
Rotor{T}(w, x, y, z) where {T<:Number} = (n = √T(w^2 + x^2 + y^2 + z^2); Rotor{T}(SVector{4,T}(w / n, x / n, y / n, z / n)))
QuatVec(w, x, y, z) = QuatVec(SVector{4}(oftype(w, false), x, y, z))
QuatVec{T}(w, x, y, z) where {T<:Number} = QuatVec{T}(SVector{4,T}(oftype(w, false), x, y, z))

# Constructor from vector components
(::Type{QT})(x, y, z) where {QT<:AbstractQuaternion} = (v = SVector{4}(false, x, y, z); QT{eltype(v)}(v))
(::Type{QT})(x, y, z) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(SVector{4,T}(false, x, y, z))
Rotor(x, y, z) = (n = √(x^2 + y^2 + z^2); Rotor(SVector{4}(false, x / n, y / n, z / n)))
Rotor{T}(x, y, z) where {T<:Number} = (n = √T(x^2 + y^2 + z^2); Rotor{T}(SVector{4,T}(false, x / n, y / n, z / n)))
QuatVec(x, y, z) = QuatVec(SVector{4}(false, x, y, z))
QuatVec{T}(x, y, z) where {T<:Number} = QuatVec{T}(SVector{4,T}(false, x, y, z))

# Constructor from scalar component
# (::Type{QT})(w::Number) where {QT<:AbstractQuaternion} = (v=SVector{4}(w, false, false, false); QT{eltype(v)}(v))
# (::Type{QT})(w::Number) where {T<:Number, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(w, false, false, false))
Quaternion(w::Number) = Quaternion(SVector{4}(w, false, false, false))
Quaternion(w::Symbolics.Num) = Quaternion(SVector{4}(w, false, false, false))
Quaternion{T}(w::Number) where {T<:Number} = Quaternion{T}(SVector{4,T}(w, false, false, false))
Rotor(w::Number) = Rotor(SVector{4}(one(w), false, false, false))
Rotor(w::Symbolics.Num) = Rotor(SVector{4}(one(w), false, false, false))
Rotor{T}(w::Number) where {T<:Number} = Rotor{T}(SVector{4,T}(one(T), false, false, false))
QuatVec(w::Number) = QuatVec(SVector{4,typeof(w)}(false, false, false, false))
QuatVec(w::Symbolics.Num) = QuatVec(SVector{4,typeof(w)}(false, false, false, false))
QuatVec{T}(w::Number) where {T<:Number} = QuatVec{T}(SVector{4,T}(false, false, false, false))

# Copy constructor
# (::Type{QT})(q::AbstractQuaternion) where {QT<:AbstractQuaternion} = QT{eltype(components(q))}(components(q))
# (::Type{QT})(q::AbstractQuaternion{S}) where {T<:Number, S<:Number, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(components(q)))
# (::Type{QT})(q::AbstractQuaternion{T}) where {T<:Number, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(components(q)))
Quaternion(q::AbstractQuaternion{T}) where {T<:Number} = Quaternion(components(q)...)
Quaternion{T}(q::AbstractQuaternion{S}) where {T<:Number,S<:Number} = Quaternion{T}(components(q)...)
Rotor(q::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Rotor(components(q)...)
Rotor{T}(q::AbstractQuaternion{S}) where {T<:Number,S<:Number} = Rotor{T}(components(q)...)
QuatVec(q::AbstractQuaternion{T}) where {T<:Number} = QuatVec(components(q)...)
QuatVec{T}(q::AbstractQuaternion{S}) where {T<:Number,S<:Number} = QuatVec{T}(components(q)...)

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

The quaternionic unit associated with rotation about the `x` axis.  Can also be entered as Unicode
bold: `𝐢`.

# Examples
```jldoctest
julia> imx * imx
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imx
 + 1.2𝐢 + 0.0𝐣 + 0.0𝐤
```
"""
const imx = QuatVec(false, true, false, false)
const 𝐢 = imx

"""
    imy

The quaternionic unit associated with rotation about the `y` axis.  Can also be entered as Unicode
bold: `𝐣`.

# Examples
```jldoctest
julia> imy * imy
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imy
 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤
```
"""
const imy = QuatVec(false, false, true, false)
const 𝐣 = imy

"""
    imz

The quaternionic unit associated with rotation about the `z` axis.  Can also be entered as Unicode
bold: `𝐤`.

# Examples
```jldoctest
julia> imz * imz
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imz
 + 0.0𝐢 + 0.0𝐣 + 1.2𝐤
```
"""
const imz = QuatVec(false, false, false, true)
const 𝐤 = imz

# Essential constructors
Base.zero(::Type{QT}) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(false, false, false, false)
Base.zero(q::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Base.zero(QT)
Base.zero(::Type{Rotor}) = throw(DomainError("Rotor", "Zero is not a possible rotor."))
Base.zero(::Type{Rotor{T}}) where {T} = throw(DomainError("Rotor", "Zero is not a possible rotor."))

Base.one(::Type{QT}) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(true, false, false, false)
Base.one(q::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Base.one(QT)
Base.one(::Type{QuatVec}) = throw(DomainError("QuatVec", "One is not a possible 3-vector."))
Base.one(::Type{QuatVec{T}}) where {T} = throw(DomainError("QuatVec", "One is not a possible 3-vector."))

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
Base.real(::Type{T}) where {T<:AbstractQuaternion} = eltype(T)
Base.real(q::AbstractQuaternion{T}) where {T<:Number} = q[1]
Base.imag(q::AbstractQuaternion{T}) where {T<:Number} = @view components(q)[2:4]
Base.vec(q::AbstractQuaternion{T}) where {T<:Number} = @view components(q)[2:4]

# Type games
wrapper(::T) where {T} = wrapper(T)
wrapper(T::UnionAll) = T
wrapper(T::Type{Q}) where {S<:Number,Q<:AbstractQuaternion{S}} = T.name.wrapper
wrapper(::Type{T}, ::Type{T}) where {T<:AbstractQuaternion} = wrapper(T)

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
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:Symbolics.Num}) where {OP} = Quaternion
        wrapper(::Type{<:Symbolics.Num}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = Quaternion
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

for QT ∈ (AbstractQuaternion, QuatVec)
    for NT ∈ (Number, Symbolics.Num)
        for OP ∈ (Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = $QT
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = $QT
            end
        end
    end
end
for QT ∈ (Rotor,)
    for NT ∈ (Number, Symbolics.Num)
        for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = Quaternion
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = Quaternion
            end
        end
    end
end
for T ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor, Number, Symbolics.Num)
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = Quaternion
        if T !== Quaternion
            @eval wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = Quaternion
        end
    end
end


Base.eltype(::Type{<:AbstractQuaternion{T}}) where {T} = T
Base.widen(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){widen(eltype(Q))}
Base.float(::Type{Q}) where {Q<:AbstractQuaternion{<:AbstractFloat}} = Q
Base.float(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){float(eltype(Q))}
Base.float(q::AbstractQuaternion{T}) where {T<:AbstractFloat} = q
Base.float(q::AbstractQuaternion{T}) where {T} = wrapper(q){float(T)}(float(components(q)))

Base.big(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){big(eltype(Q))}
Base.big(q::AbstractQuaternion{T}) where {T<:Number} = wrapper(q){big(T)}(q)

Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Number} =
    wrapper(Q){promote_type(eltype(Q), S)}
Base.promote_rule(::Type{QuatVec{T}}, ::Type{S}) where {T<:Number,S<:Number} =
    Quaternion{promote_type(T, S)}
Base.promote_rule(::Type{QuatVec{T}}, ::Type{S}) where {T<:Number,S<:AbstractQuaternion} =
    wrapper(wrapper(QuatVec), wrapper(S)){promote_type(T, eltype(S))}
Base.promote_rule(::Type{Q1}, ::Type{Q2}) where {Q1<:AbstractQuaternion,Q2<:AbstractQuaternion} =
    wrapper(wrapper(Q1), wrapper(Q2)){promote_type(eltype(Q1), eltype(Q2))}
