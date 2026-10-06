# Support for Mooncake.  Mooncake differentiates the package's source natively, so this
# extension provides rules only for the two internal functions whose implementations
# Mooncake cannot trace:
#
#   * `Quaternionic.csqrt(z)`, the complex square root, because Base's `sqrt(::Complex)`
#     reinterprets the bits of floating-point numbers, which Mooncake refuses to
#     differentiate.  Every complex square root in the package goes through `csqrt`, which
#     covers complex quaternions and `Lorentz` rotors.
#   * `Quaternionic.dominant_eigenvector_lapack(A, uplo)`, because it calls LAPACK.  It is
#     used by `from_rotation_matrix` and `align`.
#
# It also replaces, for Mooncake only, the implementations of `abs` and `absvec` of
# quaternions with real floating-point components (see the overlays below).
#
# The rules themselves are the ChainRules `frule`s and `rrule`s of the ChainRulesCore
# extension, imported with `@from_chainrules`, which provides both forward and reverse
# mode.  That is why this extension is triggered by ChainRulesCore as well as by Mooncake;
# Mooncake depends on ChainRulesCore, so the second trigger is always present.  Real
# arguments of `csqrt` are left to Mooncake's own rule for `sqrt`.
#
# This extension is loaded only as a package extension (Julia 1.9 and later), never by
# Requires.
module QuaternionicMooncakeExt

import Mooncake
import ChainRulesCore
using Quaternionic: Quaternionic, AbstractQuaternion, Quaternion, QuatVec, quaternion, quatvec,
    components, absvec

for T ∈ (Float32, Float64)
    @eval Mooncake.@from_chainrules Mooncake.DefaultCtx Tuple{
        typeof(Quaternionic.csqrt), Complex{$T}
    }
    @eval Mooncake.@from_chainrules Mooncake.DefaultCtx Tuple{
        typeof(Quaternionic.dominant_eigenvector_lapack), Matrix{$T}, Char
    }
end

# For real components, `abs` and `absvec` call `hypot` with three or four arguments.
# Mooncake's reverse rule for `hypot` with more than two arguments calls Base's `_hypot`,
# which reinterprets the bits of floating-point numbers, so that forward mode over reverse
# mode (as in `DifferentiationInterface.hessian` with `SecondOrder(AutoMooncakeForward(),
# AutoMooncake())`) throws an `ArgumentError` for every function that calls them, which
# includes `log`, `exp`, and powers.  These overlays compute the same norms with the
# package's `_hypot`, which scales the components by the largest of their absolute values
# and takes a square root, so that every mode of Mooncake differentiates them natively.
# The overlays are used only while Mooncake differentiates, and their values agree with
# those of `hypot` to within roundoff.  A zero vector gets the norm zero as a constant, with
# a zero derivative.
const IEEE = Base.IEEEFloat
scalednorm(c) = iszero(maximum(abs, c)) ? zero(eltype(c)) : Quaternionic._hypot(c)
Mooncake.@mooncake_overlay Base.abs(q::Quaternion{T}) where {T<:IEEE} = scalednorm(components(q))
Mooncake.@mooncake_overlay Base.abs(q::QuatVec{T}) where {T<:IEEE} = scalednorm(vec(q))
Mooncake.@mooncake_overlay Quaternionic.absvec(q::AbstractQuaternion{T}) where {T<:IEEE} = scalednorm(vec(q))

# Friendly tangents.  By default, Mooncake reports the gradient with respect to a
# quaternion argument in its structural form, `Tangent{(components=Tangent{(data=…)},)}`.
# With `Mooncake.Config(friendly_tangents=true)`, the methods below report it as a
# `Quaternion` instead (or as a `QuatVec` for a `QuatVec` argument), which is also the
# tangent type that the ChainRules rules use.  The gradient with respect to a `Rotor` is
# therefore a `Quaternion`, not a `Rotor`.  These hooks were introduced in Mooncake 0.5, so
# they are defined only when Mooncake provides them.
@static if isdefined(Mooncake, :AsCustomised) && isdefined(Mooncake, :FriendlyTangentCache)
    const MooncakeScalar = Union{Base.IEEEFloat,Complex{<:Base.IEEEFloat}}

    function Mooncake.friendly_tangent_cache(::AbstractQuaternion{<:MooncakeScalar})
        return Mooncake.FriendlyTangentCache{Mooncake.AsCustomised}(nothing)
    end

    # `componenttangents(t)` returns the four components of the structural tangent `t` of a
    # quaternion.
    componenttangents(t::Mooncake.Tangent) = t.fields.components.fields.data

    function Mooncake.tangent_to_friendly_internal!!(
        ::Nothing, ::AbstractQuaternion{<:MooncakeScalar}, t::Mooncake.Tangent
    )
        return quaternion(componenttangents(t)...)
    end
    function Mooncake.tangent_to_friendly_internal!!(
        ::Nothing, ::QuatVec{<:MooncakeScalar}, t::Mooncake.Tangent
    )
        c = componenttangents(t)
        return quatvec(c[2], c[3], c[4])
    end
end

end  # module
