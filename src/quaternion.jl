# We'll need this awkward way of getting the `components` field when we set `getproperty`
"""
    components(q::AbstractQuaternion{T})

Return the components of `q` as stored in the struct — as an `SVector{4, T}`.  The
components are ordered as `(w, x, y, z)`, where `w` is the scalar part and `x`, `y`, and `z`
are the vector parts.
"""
components(q::AbstractQuaternion) = getfield(q, :components)

# The square root of a complex number `z`, which is just `sqrt(z)`.  It is a separate
# function owned by this package, so that a derivative rule can be attached to it without
# type piracy.  Base's `sqrt(::Complex)` uses bit-level operations that Mooncake cannot
# differentiate, so the Mooncake extension imports the ChainRules rule for this function.
# Every complex square root in the package goes through it.
csqrt(z) = sqrt(z)

# The exponential, cosine, and sine of a complex number, computed from real functions.
# Base's methods for complex arguments take special branches when the real or imaginary part
# is exactly zero, and tracing AD backends (Mooncake, Enzyme, and ReverseDiff) differentiate
# those branches incorrectly.  For finite arguments, these formulas are the ones Base uses
# away from those branches.  The sine and cosine of each real number are computed separately
# rather than with `sincos`, because Enzyme cannot take second derivatives through `sincos`.
# For a real argument, each function is Base's own.
cexp(z) = exp(z)
function cexp(z::Complex)
    e = exp(real(z))
    Complex(e * cos(imag(z)), e * sin(imag(z)))
end
ccos(z) = cos(z)
ccos(z::Complex) = Complex(cos(real(z)) * cosh(imag(z)), -sin(real(z)) * sinh(imag(z)))
csin(z) = sin(z)
csin(z::Complex) = Complex(sin(real(z)) * cosh(imag(z)), cos(real(z)) * sinh(imag(z)))

# This helper function is mostly copied from Base.math, except that we restrict to Complex
# elements, and we omit `abs2` in the final sum.  This is crucial because it is the
# appropriate norm for complex quaternions, which are supposed to represent rotors in the
# spacetime algebra.
function _hypot(x)
    maxabs = maximum(y -> abs(value(y)), x)
    result_type = promote_type(eltype(x), typeof(maxabs))
    if isnan(maxabs) && any(isinf, x)
        return convert(result_type, Inf)
    elseif (iszero(maxabs) || isinf(maxabs))
        return convert(result_type, maxabs)
    else
        # The square is written as a product, because `^` of a `Complex` goes through
        # Base's `_cpow`, which reverse-mode AD packages such as ReverseDiff cannot handle.
        return maxabs * csqrt(sum(y -> (y / maxabs) * (y / maxabs), x))
    end
end


"""
    AbstractQuaternion{T<:Number} <: Number

Abstract supertype of the quaternion types [`Quaternion`](@ref), [`Rotor`](@ref), and
[`QuatVec`](@ref), with elements of type `T`, which may be retrieved with
[`basetype`](@ref).  Every subtype stores its four components `(w, x, y, z)` as an
`SVector{4, T}`, which may be retrieved with [`components`](@ref).
"""
AbstractQuaternion

"""
    Quaternion{T<:Number} <: AbstractQuaternion{T}

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

"""
    quaternion(w, x, y, z)
    quaternion(x, y, z)
    quaternion(w)
    quaternion(v)
    quaternion(q)

Create a [`Quaternion`](@ref) with the given components, whose element type is inferred from
the arguments.  Any missing components are set to zero.  Here, `v` is a vector of length 1,
3, or 4, which is interpreted as the corresponding list of arguments, and `q` is any
quaternion, whose components are copied.
"""
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
    Rotor{T<:Number} <: AbstractQuaternion{T}

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
    rotor(x, y, z)
    rotor(w)
    rotor(v)
    rotor(q)

create a new rotor with the given components (where the components are as described in
[`Quaternion`](@ref)), automatically normalizing them on input.  Here, `v` is a vector of
length 1, 3, or 4, and `q` is any quaternion.  The element type of the result is a
floating-point type, even for integer input.  The type constructors `Rotor(...)`, without a
type parameter, are identical to `rotor(...)` and also normalize.  So does
`convert(Rotor{T}, x)` when `x` is a number or a quaternion other than a `Rotor`, and
therefore so do the implicit conversions made by `push!`, `setindex!`, and typed array
literals.  Converting a `Rotor` to a `Rotor{T}` with a different element type `T` only
converts the components, without normalizing them again.

If you would like to bypass normalization, you can call the type constructor with an
explicit type parameter:

    Rotor{T}(w, x, y, z)
    Rotor{T}(x, y, z)
    Rotor{T}(w)
    Rotor{T}(v)
    Rotor{T}(q)

These simply convert the given components to type `T` and store them; here, `v` is any
`AbstractVector` that can be converted to an `SVector{4, T}`.  If you want to handle the
normalization step yourself, you can use `LinearAlgebra.normalize`.

However, once a `Rotor` is created, its norm will often be *assumed* to be precisely 1.  So
if its true norm is significantly different, you will likely see weird results — including
vectors with very different lengths after "rotation" by a non-unit `Rotor`.  For the same
reason, `abs` and `abs2` of a `Rotor` return exactly 1, and mixing a `Rotor` with any number
or quaternion that is not a `Rotor`, in arithmetic or in `promote`, gives a `Quaternion`.
The functions `zero(Rotor{T})` also returns a `Quaternion`, and `round` of a `Rotor` returns
a `Quaternion`, because neither result is a unit quaternion.  However, `zeros(Rotor{T},
dims...)` returns an array with element type `Rotor{T}`, filled with zero rotors that are
not normalized, so that it can be used to preallocate an array of rotors.

The product or quotient of two `Rotor`s is a `Rotor`.  For real element types, that result
is renormalized to remove the accumulated rounding error.  For complex element types (as for
[`Lorentz`](@ref) transformations), it is renormalized only when the sum of the squared
magnitudes of its components is less than 16.  The complex spinor norm of a boost with
rapidity ``η`` is computed as ``\cosh^2(η/2) - \sinh^2(η/2)``, which cancels
catastrophically for large ``η``, so renormalizing such a product would add error rather
than remove it.

For complex element types, normalization divides by the principal square root of the spinor
norm ``w^2 + x^2 + y^2 + z^2``, which is itself complex.  The overall sign of the normalized
result therefore changes discontinuously as the spinor norm crosses the negative real axis.

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
julia> rotor(Rotor{Float16}(1, 2, 3, 4))
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

# The division is written out component by component, rather than as a broadcast `a ./ n`,
# because ReverseDiff's tapes cannot replay a broadcast into an immutable `SVector`.
"""
    rotor(w, x, y, z)
    rotor(x, y, z)
    rotor(w)
    rotor(v)
    rotor(q)

Create a [`Rotor`](@ref) from the given components, normalizing them so that the result has
norm 1.  Missing components are set to zero, `v` is a vector of length 1, 3, or 4, which is
interpreted as the corresponding list of arguments, and `q` is any quaternion.  The element
type of the result is a floating-point type, even for integer input.  To store components
without normalizing them, use the type constructor `Rotor{T}(...)` with an explicit type
parameter.
"""
function rotor(a::SVector{4,T}) where {T<:Number}
    n = abs(Quaternion{T}(a))
    normalized = SVector{4}(a[1] / n, a[2] / n, a[3] / n, a[4] / n)
    Rotor{eltype(normalized)}(normalized)
end
#rotor(a::AbstractVector) = Rotor{T}(SVector{4,T}(a))  # See below
rotor(a::AbstractQuaternion) = rotor(components(a))
rotor(w, x, y, z) = rotor(SVector{4}(w, x, y, z))
rotor(x,y,z) = rotor(false, x,y,z)
# Every arity shares the same normalization, so that `rotor(0)` is NaN like `rotor(0, 0, 0,
# 0)`, and the element type does not depend on the number of arguments.
rotor(w::Number) = rotor(w, false, false, false)

# This builds the `Rotor` that results from a product or quotient of rotors, at least one of
# which is complex, as for `Lorentz` transformations; see `wrapper` below.  Renormalizing
# divides by the complex spinor norm, which is computed with an absolute error of roughly
# eps·S, where S is the sum of the squared magnitudes of the components (cosh η for a boost
# of rapidity η).  Below the bound S < 16, that error is at most several ulps, and
# renormalizing removes the drift that accumulates over long chains of products.  Above it,
# the error of renormalizing grows in proportion to S and exceeds the error of the product
# itself, so the components are stored as they are.  The test is made on `value(S)`, so that
# dual numbers take the same branch as their values.  A comparison that does not return a
# `Bool`, as for symbolic types, skips renormalization.
function complex_rotor_product(a::SVector{4})
    if (value(sum(abs2, a)) < 16) === true
        rotor(a)
    else
        Rotor{eltype(a)}(a)
    end
end
complex_rotor_product(w, x, y, z) = complex_rotor_product(SVector{4}(w, x, y, z))

Rotor(w::SVector{4,T}) where {T} = rotor(w)
Rotor(w::AbstractVector) = rotor(w)
Rotor(w::AbstractQuaternion) = rotor(w)
Rotor(w,x,y,z) = rotor(w,x,y,z)
Rotor(x,y,z) = rotor(x,y,z)
Rotor(w::T) where {T<:Number} = rotor(w)

"""
    QuatVec{T<:Number} <: AbstractQuaternion{T}

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
[`Quaternion`](@ref)), except that the scalar argument `w` is always set to 0.  The same is
true of the type constructors `QuatVec(...)` and `QuatVec{T}(...)`, and therefore also of
`convert`, `push!`, and other implicit conversions: any scalar part given to them is
discarded, so that a `QuatVec` never stores a nonzero scalar part.  For example,
`QuatVec{T}(w, x, y, z)` and `convert(QuatVec{T}, q)` ignore `w` and the scalar part of `q`,
respectively, and `QuatVec{T}(w)` is the zero vector.

A `QuatVec` is equal (by `==` and `isequal`) to another quaternion only if the other
quaternion has zero scalar part and the same vector part, and it is equal to a number only
if both are zero.  Because the multiplicative identity is not a pure vector,
`one(QuatVec{T})` and `oneunit(QuatVec{T})` return `Quaternion` values.  On the other hand,
`ones(QuatVec{T}, dims...)` returns an array with element type `QuatVec{T}`, as requested,
so the scalar part of each `one(QuatVec{T})` is discarded when it is stored, and every
element is the zero vector.

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
    # Every constructor discards the scalar part, so that a `QuatVec` never stores a
    # nonzero scalar component.  This includes the implicit conversions made by
    # `convert`, `push!`, `setindex!`, and typed array literals.
    QuatVec{T}(a::SVector{4,T}) where {T<:Number} = new{T}(SVector{4,T}(false, a[2], a[3], a[4]))
    QuatVec{T}(a::AbstractVector) where {T<:Number} = QuatVec{T}(SVector{4,T}(a))
    QuatVec{T}(a::AbstractQuaternion) where {T<:Number} = QuatVec{T}(SVector{4,T}(components(a)))
    QuatVec{T}(_,x,y,z) where {T<:Number} = new{T}(SVector{4,T}(false,x,y,z))
    QuatVec{T}(x,y,z) where {T<:Number} = new{T}(SVector{4,T}(false,x,y,z))
    QuatVec{T}(::Number) where {T<:Number} = new{T}(SVector{4,T}(false, false, false, false))
end

"""
    quatvec(w, x, y, z)
    quatvec(x, y, z)
    quatvec(w)
    quatvec(v)
    quatvec(q)

Create a [`QuatVec`](@ref) with the given vector components, whose element type is inferred
from the arguments.  The scalar argument `w`, or the scalar part of the quaternion `q`, is
discarded, so `quatvec(w)` is the zero vector.  Here, `v` is a vector of length 1, 3, or 4,
which is interpreted as the corresponding list of arguments.
"""
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
"""
    𝐢

Alias for [`imx`](@ref), the quaternionic unit `QuatVec{Bool}` along the `x` axis.
"""
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
"""
    𝐣

Alias for [`imy`](@ref), the quaternionic unit `QuatVec{Bool}` along the `y` axis.
"""
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
"""
    𝐤

Alias for [`imz`](@ref), the quaternionic unit `QuatVec{Bool}` along the `z` axis.
"""
const 𝐤 = imz

# Essential constructors
Base.zero(::Type{QT}) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(false, false, false, false)
Base.zero(::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Base.zero(QT)
Base.zero(::Type{Rotor{T}}) where {T} = zero(Quaternion{T})

Base.one(::Type{QT}) where {T<:Number,QT<:AbstractQuaternion{T}} = QT(true, false, false, false)
Base.one(::QT) where {T<:Number,QT<:AbstractQuaternion{T}} = Base.one(QT)
Base.one(::Type{QuatVec{T}}) where {T} = one(Quaternion{T})

# Base's fallbacks for the `UnionAll` types `Rotor` and `QuatVec` go through `convert`,
# which would give the wrong identity element.  These methods follow the convention of the
# concrete types above: a quantity that cannot be stored in the original type is returned as
# a `Quaternion`.
Base.zero(::Type{Rotor}) = zero(Quaternion)
Base.one(::Type{QuatVec}) = one(Quaternion)
Base.oneunit(::Type{QuatVec{T}}) where {T} = one(Quaternion{T})
Base.oneunit(::Type{QuatVec}) = one(Quaternion)
Base.oneunit(q::QuatVec) = oneunit(typeof(q))

# Base's `zeros` fills an array of the requested element type with `zero(T)`.  For `Rotor`,
# that value is a `Quaternion`, which would then be converted back to a `Rotor`, and that
# conversion normalizes, which would give NaN.  Instead, as in earlier versions, the array
# keeps the element type `Rotor{T}` and is filled with the zero rotor, constructed without
# normalization, so that it can be used for preallocation.  Every other method of `zeros`
# with a type argument reaches one of these two.
Base.zeros(::Type{Rotor{T}}, dims::NTuple{N,Integer}) where {T,N} =
    fill(Rotor{T}(false, false, false, false), dims)
Base.zeros(::Type{Rotor{T}}, dims::Tuple{}) where {T} =
    fill(Rotor{T}(false, false, false, false), dims)

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
Base.getindex(q::AbstractQuaternion, ::Colon) = [q[i] for i in 1:4]
Base.getindex(q::AbstractQuaternion, ::CartesianIndex{0}) = q
# Indexing with `begin` and `end` refers to the four components, as for `q[i]`, rather than
# treating the quaternion as a single scalar, as Base's methods for `Number` do.
Base.firstindex(::AbstractQuaternion) = 1
Base.lastindex(::AbstractQuaternion) = 4
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
# The product or quotient of two unit rotors is a unit rotor up to rounding error, so for
# real elements it is merely renormalized (through `Rotor(w, x, y, z)`) to remove drift.
# For complex elements, as for `Lorentz` transformations, renormalization divides by the
# complex spinor norm.  For a boost of rapidity η, the spinor norm is computed as cosh²(η/2)
# - sinh²(η/2), which cancels catastrophically: the relative error is about eps·cosh(η), and
# the result is NaN for η ≳ 38 in `Float64`.  So complex products are renormalized only when
# that error is small; see `complex_rotor_product` above.
for OP ∈ (Val{*}, Val{/})
    @eval begin
        wrapper(::Type{<:Rotor{<:Complex}}, ::$OP, ::Type{<:Rotor}) = complex_rotor_product
        wrapper(::Type{<:Rotor}, ::$OP, ::Type{<:Rotor{<:Complex}}) = complex_rotor_product
        wrapper(::Type{<:Rotor{<:Complex}}, ::$OP, ::Type{<:Rotor{<:Complex}}) = complex_rotor_product
    end
end
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
# These resolve ambiguities between the generic methods for `AbstractQuaternion` and the
# methods for `Number` above, which would otherwise make `*` and `/` between a `Rotor` or a
# `QuatVec` and a user-defined subtype of `AbstractQuaternion` throw.
for QT ∈ (Rotor, QuatVec)
    for OP ∈ (Val{*}, Val{/})
        @eval begin
            wrapper(::Type{<:AbstractQuaternion}, ::$OP, ::Type{<:$QT}) = Quaternion
            wrapper(::Type{<:$QT}, ::$OP, ::Type{<:AbstractQuaternion}) = Quaternion
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


"""
    basetype(q)
    basetype(QT)

Return the element type `T` of the quaternion `q`, or of the quaternion type `QT`, where `q
isa AbstractQuaternion{T}` or `QT <: AbstractQuaternion{T}`.  This is the type of each
component, as `eltype` is for an array.  For example, `basetype(QuaternionF64)` is
`Float64`.
"""
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
# A `Rotor` mixed with any other number (including `Complex`) is promoted to `Quaternion`,
# because the other number need not have unit norm.
Base.promote_rule(::Type{Rotor{T}}, ::Type{S}) where {T<:Number,S<:Number} =
    Quaternion{promote_type(T, S)}
Base.promote_rule(::Type{Rotor{T}}, ::Type{S}) where {T<:Number,S<:AbstractQuaternion} =
    wrapper(wrapper(Rotor), wrapper(S)){promote_type(T, basetype(S))}
Base.promote_rule(::Type{Q1}, ::Type{Q2}) where {Q1<:AbstractQuaternion,Q2<:AbstractQuaternion} =
    wrapper(wrapper(Q1), wrapper(Q2)){promote_type(basetype(Q1), basetype(Q2))}

# Conversion to `Rotor{T}` normalizes, as `convert(Rotor, x)` and `rotor(x)` do, so that the
# implicit conversions made by `push!`, `setindex!`, `fill!`, and typed array literals never
# create a non-unit `Rotor`.  The components are converted to `T` before normalizing, so
# that the norm is 1 to the precision of `T`.  A `Rotor` is already normalized, so only its
# element type is converted.  The raw constructor `Rotor{T}(...)` remains the way to store
# components without normalizing them.
Base.convert(::Type{Rotor{T}}, x::Number) where {T<:Number} = Rotor{T}(rotor(Quaternion{T}(x)))
Base.convert(::Type{Rotor{T}}, q::AbstractQuaternion) where {T<:Number} =
    Rotor{T}(rotor(Quaternion{T}(q)))
Base.convert(::Type{Rotor{T}}, q::Rotor) where {T<:Number} = Rotor{T}(q)
Base.convert(::Type{Rotor{T}}, q::Rotor{T}) where {T<:Number} = q

# Like a `Complex`, a quaternion can be converted to a real type only when it is real.  This
# is restricted to the standard real types, because packages such as ForwardDiff and
# Symbolics define constructors of their own real types from any argument, which would
# otherwise be ambiguous.
(::Type{T})(q::AbstractQuaternion) where {T<:Union{AbstractFloat,Integer,Rational}} =
    isreal(q) ? T(real(q))::T : throw(InexactError(nameof(T), T, q))
