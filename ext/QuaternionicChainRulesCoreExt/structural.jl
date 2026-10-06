# Rules for structural access (`components`, `getindex`, `real`, `imag`, and `vec`), for the
# constructors, and for the normalizing functions `rotor` and `normalize`.
#
# Each constructor rule is the exact derivative of what the constructor stores.  The type
# constructors with an explicit parameter, such as `Rotor{T}(w, x, y, z)`, store their
# arguments without normalizing them; `QuatVec` constructors discard any scalar part, so a
# scalar argument receives a zero cotangent; and `rotor`, `Rotor` without a parameter, and
# `normalize` divide by the norm computed with `hypot`.  The structural `q.components` (through
# `getfield`) needs no rule, because its cotangent arrives as a `Tangent` or a `NamedTuple`,
# which `cot` and `ProjectTo` accept.  Properties such as `q.w` and `q.im` reach `getindex`.

########################################################################################
## Structural access
########################################################################################

function rrule(::typeof(components), q::AQ)
    proj = ProjectTo(q)
    function components_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(Δ)
    end
    return components(q), components_pullback
end
function frule((_, q̇)::Tuple, ::typeof(components), q::AQ)
    Ω = components(q)
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    Δp = ProjectTo(q)(Δ)
    Δp isa AbstractZero && return Ω, Δp
    return Ω, components(Δp)
end

function rrule(::typeof(getindex), q::AQ, i::Integer)
    proj = ProjectTo(q)
    function getindex_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, NoTangent())
        return NoTangent(), proj(componentcot(i, Δ)), NoTangent()
    end
    return q[i], getindex_pullback
end
function frule((_, q̇, _)::Tuple, ::typeof(getindex), q::AQ, i::Integer)
    Ω = q[i]
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    return Ω, Δ[i]
end

# The cotangent of `real(q)` is placed in the scalar component as it is.  Taking `real` of
# it would be wrong for complex components, whose scalar part is itself complex.
function rrule(::typeof(real), q::AQ)
    proj = ProjectTo(q)
    function real_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(componentcot(1, Δ))
    end
    return real(q), real_pullback
end
function frule((_, q̇)::Tuple, ::typeof(real), q::AQ)
    Ω = real(q)
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    return Ω, Δ[1]
end

"""
    vectorpartcot(Δ)

Normalize a cotangent that a backend may pass for the three-component view returned by
`imag(q)` or `vec(q)` into the `Quaternion` with zero scalar part and that vector part, or
into `ZeroTangent()`.  It accepts a `Tuple` or an `AbstractVector` of length 3, and the
structural `Tangent` (or `NamedTuple`) of a `SubArray`, whose `parent` field holds the
cotangent of all four components.
"""
vectorpartcot(Δ::AbstractThunk) = vectorpartcot(unthunk(Δ))
vectorpartcot(::Union{Nothing,AbstractZero}) = ZeroTangent()
vectorpartcot(Δ::Tangent) = vectorpartcot(ChainRulesCore.backing(Δ))
function vectorpartcot(Δ::NamedTuple)
    haskey(Δ, :parent) ||
        throw(ArgumentError("Cannot interpret $(typeof(Δ)) as the cotangent of a vector part"))
    P = cot(Δ.parent)
    P isa AbstractZero && return P
    return quaternion(false, P[2], P[3], P[4])
end
function vectorpartcot(Δ::Tuple)
    length(Δ) == 3 ||
        throw(DimensionMismatch("A vector-part cotangent must have 3 components, not $(length(Δ))"))
    return cotcomponents(nothing, Δ[1], Δ[2], Δ[3])
end
function vectorpartcot(Δ::AbstractVector)
    length(Δ) == 3 ||
        throw(DimensionMismatch("A vector-part cotangent must have 3 components, not $(length(Δ))"))
    return cotcomponents(nothing, Δ[begin], Δ[begin+1], Δ[begin+2])
end

for f ∈ (:imag, :vec)
    @eval begin
        function rrule(::typeof($f), q::AQ)
            proj = ProjectTo(q)
            function vectorpart_pullback(ΔΩ)
                Δ = vectorpartcot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(Δ)
            end
            return $f(q), vectorpart_pullback
        end
        function frule((_, q̇)::Tuple, ::typeof($f), q::AQ)
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, vec(Δ)
        end
    end
end

########################################################################################
## Splitting and assembling the components of constructor arguments
########################################################################################

"""
    splitcot(QT, Δ, Val(N))

Split the `Quaternion` cotangent `Δ` of a quaternion of type `QT`, constructed from `N`
numbers, into a tuple of the cotangents of those numbers.  Four numbers are the components
`(w, x, y, z)`, three are `(x, y, z)`, and one is `w`.  A `QuatVec` discards its scalar
argument, which therefore receives `ZeroTangent()`.
"""
splitcot(::Type, Δ, ::Val{4}) = (Δ[1], Δ[2], Δ[3], Δ[4])
splitcot(::Type, Δ, ::Val{3}) = (Δ[2], Δ[3], Δ[4])
splitcot(::Type, Δ, ::Val{1}) = (Δ[1],)
splitcot(::Type{<:QuatVec}, Δ, ::Val{4}) = (ZeroTangent(), Δ[2], Δ[3], Δ[4])
splitcot(::Type{<:QuatVec}, Δ, ::Val{1}) = (ZeroTangent(),)

"""
    vectorcot(QT, v, Δ)

Return the cotangent of the vector `v` from which a quaternion of type `QT` was constructed,
given the `Quaternion` cotangent `Δ` of that quaternion, as an array projected onto `v`.  A
vector of length 4 holds the components `(w, x, y, z)`, one of length 3 holds `(x, y, z)`,
and one of length 1 holds `w`.  The scalar entry of a vector from which a `QuatVec` was
constructed receives a zero.
"""
function vectorcot(::Type{QT}, v::AbstractVector, Δ) where {QT}
    Δ isa AbstractZero && return Δ
    n = length(v)
    z = zero(Δ[1])
    dropscalar = QT <: QuatVec
    d = if n == 4
        [dropscalar ? z : Δ[1], Δ[2], Δ[3], Δ[4]]
    elseif n == 3
        [Δ[2], Δ[3], Δ[4]]
    else
        [dropscalar ? z : Δ[1]]
    end
    # `ProjectTo` of a `StaticArray` does not convert the element type, so the elements are
    # projected first.
    T = eltype(v)
    d′ = isconcretetype(T) && T <: Number ? map(ProjectTo(zero(T)), d) : d
    return ProjectTo(v)(d′)
end

"""
    assembletangent(ȧs)

Return the `Quaternion` tangent of a quaternion constructed from numbers whose tangents are
the entries of the tuple `ȧs`, interpreted as in `splitcot`, or `ZeroTangent()` if every
entry is zero.  A tangent may also be an `AbstractVector` of length 1, 3, or 4, for a
quaternion constructed from a vector.
"""
function assembletangent(ȧs::Tuple)
    all(x -> unthunk(x) isa Union{Nothing,AbstractZero}, ȧs) && return ZeroTangent()
    return assemblecomponents(map(zerofill, ȧs)...)
end
function assembletangent(v̇::AbstractVector)
    n = length(v̇)
    n ∈ (1, 3, 4) ||
        throw(DimensionMismatch("A quaternion tangent must have 1, 3, or 4 components, not $n"))
    if n == 4
        return assembletangent((v̇[begin], v̇[begin+1], v̇[begin+2], v̇[begin+3]))
    elseif n == 3
        return assembletangent((v̇[begin], v̇[begin+1], v̇[begin+2]))
    else
        return assembletangent((v̇[begin],))
    end
end
assembletangent(v̇::AbstractThunk) = assembletangent(unthunk(v̇))
assembletangent(::Union{Nothing,AbstractZero}) = ZeroTangent()

"""
    assemblecomponents(w, x, y, z)
    assemblecomponents(x, y, z)
    assemblecomponents(w)

Return the `Quaternion` with the given components, missing components being zero, as
`quaternion` does for numbers.
"""
assemblecomponents(w, x, y, z) = quaternion(w, x, y, z)
assemblecomponents(x, y, z) = quaternion(false, x, y, z)
assemblecomponents(w) = quaternion(w, false, false, false)

########################################################################################
## Constructors that store their arguments
########################################################################################

"""
    storing_rule(QT, Ω, args)

Return the primal `Ω` and the pullback of a constructor of the quaternion type `QT` that
stores the numbers `args` as components (subject only to the `QuatVec` discarding its scalar
part), as the type constructors and `quaternion` and `quatvec` do.
"""
function storing_rule(::Type{QT}, Ω, args::Tuple{Vararg{Number,N}}) where {QT,N}
    proj = ProjectTo(Ω)
    function constructor_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), map(a -> ZeroTangent(), args)...)
        Δp = proj(Δ)
        Δp isa AbstractZero && return (NoTangent(), map(a -> Δp, args)...)
        return (NoTangent(), projectargs(args, splitcot(QT, Δp, Val(N)))...)
    end
    return Ω, constructor_pullback
end

"""
    storing_vector_rule(QT, Ω, v)

Return the primal `Ω` and the pullback of a constructor of the quaternion type `QT` from the
vector `v`, of length 1, 3, or 4, which stores the entries of `v` as components (subject only
to the `QuatVec` discarding its scalar part).
"""
function storing_vector_rule(::Type{QT}, Ω, v::AbstractVector) where {QT}
    proj = ProjectTo(Ω)
    function constructor_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        Δp = proj(Δ)
        Δp isa AbstractZero && return (NoTangent(), Δp)
        return NoTangent(), vectorcot(QT, v, Δp)
    end
    return Ω, constructor_pullback
end

"""
    converting_rule(Ω, q)

Return the primal `Ω` and the pullback of a constructor that converts the quaternion `q` to
the type of `Ω`, without normalizing it.
"""
function converting_rule(Ω, q::AQ)
    projΩ = ProjectTo(Ω)
    projq = ProjectTo(q)
    function constructor_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        Δp = projΩ(Δ)
        Δp isa AbstractZero && return (NoTangent(), Δp)
        return NoTangent(), projq(Δp)
    end
    return Ω, constructor_pullback
end

for QT ∈ (:Quaternion, :Rotor, :QuatVec)
    @eval begin
        rrule(::Type{$QT{T}}, args::Vararg{Number,N}) where {T,N} =
            storing_rule($QT, $QT{T}(args...), args)
        rrule(::Type{$QT{T}}, v::AbstractVector) where {T} =
            storing_vector_rule($QT, $QT{T}(v), v)
        rrule(::Type{$QT{T}}, q::AQ) where {T} = converting_rule($QT{T}(q), q)

        function frule((_, ȧs...)::Tuple, ::Type{$QT{T}}, args::Vararg{Number,N}) where {T,N}
            Ω = $QT{T}(args...)
            Δ = assembletangent(ȧs)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(Δ)
        end
        function frule((_, v̇)::Tuple, ::Type{$QT{T}}, v::AbstractVector) where {T}
            Ω = $QT{T}(v)
            Δ = assembletangent(v̇)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(Δ)
        end
        function frule((_, q̇)::Tuple, ::Type{$QT{T}}, q::AQ) where {T}
            Ω = $QT{T}(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(ProjectTo(q)(Δ))
        end
    end
end

# The lower-case constructors, and the type constructors without a parameter, promote their
# arguments to a common element type and then store them.
for (f, QT) ∈ ((:quaternion, :Quaternion), (:quatvec, :QuatVec))
    @eval begin
        rrule(::Union{typeof($f),Type{$QT}}, args::Vararg{Number,N}) where {N} =
            storing_rule($QT, $f(args...), args)
        rrule(::Union{typeof($f),Type{$QT}}, v::AbstractVector) = storing_vector_rule($QT, $f(v), v)
        rrule(::Union{typeof($f),Type{$QT}}, q::AQ) = converting_rule($f(q), q)

        function frule((_, ȧs...)::Tuple, ::Union{typeof($f),Type{$QT}}, args::Vararg{Number,N}) where {N}
            Ω = $f(args...)
            Δ = assembletangent(ȧs)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(Δ)
        end
        function frule((_, v̇)::Tuple, ::Union{typeof($f),Type{$QT}}, v::AbstractVector)
            Ω = $f(v)
            Δ = assembletangent(v̇)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(Δ)
        end
        function frule((_, q̇)::Tuple, ::Union{typeof($f),Type{$QT}}, q::AQ)
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(ProjectTo(q)(Δ))
        end
    end
end

########################################################################################
## Normalizing constructors and `normalize`
########################################################################################

# `rotor(w, x, y, z)`, `rotor(x, y, z)`, `rotor(v)`, and `rotor(q)` divide the components by
# `abs` of the quaternion they form; `Rotor(...)` without a type parameter is `rotor(...)`.
# The result of `rotor(w)` with a single number is ±1 (for real `w`), so its derivative is
# zero.

"""
    normalizing_rule(Ω, a, args)

Return the primal `Ω` and the pullback of a function that normalizes the quaternion `a`,
whose components are formed from the numbers `args` as by `quaternion(args...)`.
"""
function normalizing_rule(Ω, a::AQ, args::Tuple{Vararg{Number,N}}) where {N}
    function rotor_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), map(x -> ZeroTangent(), args)...)
        return (NoTangent(), projectargs(args, splitcot(Quaternion, normalize_pullback(a, Δ), Val(N)))...)
    end
    return Ω, rotor_pullback
end

const RotorConstructor = Union{typeof(rotor),Type{Rotor}}

function rrule(::RotorConstructor, w::Number, x::Number, y::Number, z::Number)
    return normalizing_rule(rotor(w, x, y, z), quaternion(w, x, y, z), (w, x, y, z))
end
function rrule(::RotorConstructor, x::Number, y::Number, z::Number)
    return normalizing_rule(rotor(x, y, z), quaternion(false, x, y, z), (x, y, z))
end
function rrule(::RotorConstructor, w::Number)
    rotor_constant_pullback(ΔΩ) = (NoTangent(), ZeroTangent())
    return rotor(w), rotor_constant_pullback
end
function rrule(::RotorConstructor, v::AbstractVector)
    Ω = rotor(v)
    a = quaternion(v)
    n = length(v)
    function rotor_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        # A vector of length 1 gives ±1, as `rotor(w)` does, so its cotangent is zero.
        ∂a = n == 1 ? zero(normalize_pullback(a, Δ)) : normalize_pullback(a, Δ)
        return NoTangent(), vectorcot(Quaternion, v, ∂a)
    end
    return Ω, rotor_pullback
end

# `rotor(q)`, `Rotor(q)`, and `normalize(q)` (which is `q / abs(q)`, or `rotor(q)` for a
# `Rotor`) share the normalization pullback.
for F ∈ (:RotorConstructor, :(typeof(normalize)))
    @eval function rrule(f::$F, q::AQ)
        Ω = f(q)
        proj = ProjectTo(q)
        function normalize_q_pullback(ΔΩ)
            Δ = cot(ΔΩ)
            Δ isa AbstractZero && return (NoTangent(), Δ)
            return NoTangent(), proj(normalize_pullback(q, Δ))
        end
        return Ω, normalize_q_pullback
    end
    @eval function frule((_, q̇)::Tuple, f::$F, q::AQ)
        Ω = f(q)
        Δ = cot(q̇)
        Δ isa AbstractZero && return Ω, Δ
        Δq = ProjectTo(q)(Δ)
        Δq isa AbstractZero && return Ω, Δq
        return Ω, ProjectTo(Ω)(normalize_pushforward(q, Δq))
    end
end

function frule((_, ẇ, ẋ, ẏ, ż)::Tuple, ::RotorConstructor, w::Number, x::Number, y::Number, z::Number)
    Ω = rotor(w, x, y, z)
    ȧ = assembletangent(map(zerofill, projectargs((w, x, y, z), (ẇ, ẋ, ẏ, ż))))
    ȧ isa AbstractZero && return Ω, ȧ
    return Ω, ProjectTo(Ω)(normalize_pushforward(quaternion(w, x, y, z), ȧ))
end
function frule((_, ẋ, ẏ, ż)::Tuple, ::RotorConstructor, x::Number, y::Number, z::Number)
    Ω = rotor(x, y, z)
    ȧ = assembletangent(map(zerofill, projectargs((x, y, z), (ẋ, ẏ, ż))))
    ȧ isa AbstractZero && return Ω, ȧ
    return Ω, ProjectTo(Ω)(normalize_pushforward(quaternion(false, x, y, z), ȧ))
end
frule(::Tuple, ::RotorConstructor, w::Number) = rotor(w), ZeroTangent()
function frule((_, v̇)::Tuple, ::RotorConstructor, v::AbstractVector)
    Ω = rotor(v)
    ȧ = assembletangent(v̇)
    ȧ isa AbstractZero && return Ω, ȧ
    Ω̇ = ProjectTo(Ω)(normalize_pushforward(quaternion(v), ȧ))
    # A vector of length 1 gives ±1, as `rotor(w)` does, so its tangent is zero.
    return Ω, length(v) == 1 ? zero(Ω̇) : Ω̇
end
