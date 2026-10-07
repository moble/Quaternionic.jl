# DifferentiationInterface treats every `Number` as a real or complex scalar.  For
# quaternions, that would silently give wrong results in two places, which this extension
# corrects.  Both concern only the input and output of the function being differentiated;
# quaternions inside the function are unaffected.
#
# 1. A reverse-mode backend computes `derivative` (or `pushforward`) of a function with a
#    scalar output by a single pullback, seeded with `oneunit(y)`.  For a quaternion output,
#    that gives only the derivative of its scalar part.  DifferentiationInterface handles
#    `Complex` outputs with a second pullback, seeded with `im`; the method below does the
#    same for quaternions, with one pullback for each of the four units (and four more,
#    seeded with `im` times the units, for complex components).
#
# 2. DifferentiationInterface builds seeds for arrays (for reverse-mode Jacobians of array
#    outputs, and for forward-mode derivatives with respect to array inputs) from `basis(a,
#    i)`, which puts `oneunit` of the element type at index `i`.  For an array of
#    quaternions, those seeds cover only the scalar parts, and the results would silently
#    omit the derivatives of (or with respect to) the vector parts.  The method below raises
#    an error instead.
#
# Both methods extend internal functions of DifferentiationInterface, and the tests check
# that they are still reached.  `Test.detect_ambiguities` reports the first method as
# ambiguous with the in-place methods of `_value_and_pushforward_via_pullback`, exactly as
# it does for DifferentiationInterface's own `Number` and `Complex` methods; no call can
# reach those ambiguities, because their first argument would have to be a quaternion and a
# function at once.
module QuaternionicDifferentiationInterfaceExt

import Quaternionic: AbstractQuaternion, Quaternion, QuatVec, basetype
import DifferentiationInterface as DI
import LinearAlgebra: dot

# The four units `c`, `c𝐢`, `c𝐣`, and `c𝐤`, as quaternions of type `S`.  The seeds of the
# pullbacks must have the type of `oneunit(y)`, the seed with which DifferentiationInterface
# prepared them; for a `Rotor` output, that is a `Rotor` (built here without normalization).
function units(::Type{S}, c) where {S<:AbstractQuaternion}
    z = zero(c)
    (S(c, z, z, z), S(z, c, z, z), S(z, z, c, z), S(z, z, z, c))
end

function DI._value_and_pushforward_via_pullback(
    y_ex::AbstractQuaternion,
    f::F,
    pullback_prep::DI.PullbackPrep,
    backend::DI.ADTypes.AbstractADType,
    x,
    tx::NTuple{B},
    contexts::Vararg{DI.Context,C},
) where {F,B,C}
    S = typeof(oneunit(y_ex))
    T = basetype(S)
    pullbackof(dy) = only(DI.pullback(f, pullback_prep, backend, x, (dy,), contexts...))
    y = f(x, map(DI.unwrap, contexts)...)
    # As with ForwardDiff, the derivative of a `QuatVec` is a `QuatVec`, and the derivative of
    # any other quaternion is a `Quaternion`.
    D = y_ex isa QuatVec ? QuatVec : Quaternion
    a = map(pullbackof, units(S, one(T)))
    ty = if T <: Complex
        b = map(pullbackof, units(S, im * one(T)))
        map(tx) do dx
            D(ntuple(k -> complex(real(dot(a[k], dx)), real(dot(b[k], dx))), 4)...)
        end
    else
        map(dx -> D(ntuple(k -> real(dot(a[k], dx)), 4)...), tx)
    end
    return y, ty
end

function DI.basis(::AbstractArray{<:AbstractQuaternion}, _)
    throw(ArgumentError(
        "DifferentiationInterface treats each element of an array as a real or complex " *
        "number, so it cannot compute derivatives of (or with respect to) the vector parts " *
        "of an array of quaternions.  Return (or pass) their components as a real array " *
        "instead, for example with `to_float_array` (or `from_float_array`)."
    ))
end

end # module
