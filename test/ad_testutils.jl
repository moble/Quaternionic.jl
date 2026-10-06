# Shared utilities for the automatic-differentiation tests (the `ad_*.jl` files)
#
# Every test item that needs these utilities declares `setup=[ADTestUtils]` and then calls
# `using .ADTestUtils`, which brings every name listed below into scope.  A test item that
# assigns one of these names itself (`R = ...`, for example) should instead import only what
# it uses, as in `using .ADTestUtils: bigfd, sc`.  The interface is as follows.
#
# References, computed by central differences in 256-bit `BigFloat` arithmetic.  None of
# them uses an AD backend, and their errors are far below Float64 roundoff (about 1e-47
# relative to the function values for first derivatives and 1e-29 for the Hessian).
#
#   bigfd(f, x::AbstractVector{<:Real})   Gradient of `f` (a `Vector`) if `f(x)` is real;
#                                         otherwise the Jacobian (a `Matrix`) of
#                                         `realcoords(f(x))` with respect to `x`.
#   bigfd(f, t::Real)                     Derivative of `f` at `t`, in the tangent type of
#                                         `f(t)` (a `Quaternion` for a `Rotor` output).
#   bigfdH(f, x::AbstractVector{<:Real})  Hessian (a `Matrix`) of the real function `f`.
#   bigjvp(f, args::Tuple, ẋs::Tuple)     Directional derivative of `f(args...)` along the
#                                         tangents `ẋs`, in the tangent type of the output:
#                                         the reference for an `frule`.
#   bigvjp(f, args::Tuple, ȳ)             Cotangents (a `Tuple`, one per argument) pulled
#                                         back from the output cotangent `ȳ`: the
#                                         reference for an `rrule`.  Integer, `Bool`, and
#                                         non-numeric arguments get `NoTangent()`.
#   bigjvp(f, x, ẋ), bigvjp(f, x, ȳ)      The same for a single argument `x` that is not a
#                                         `Tuple`; `bigvjp` then returns one cotangent.
#
#   Each result is converted to the floating-point type of the corresponding output or
#   argument at its original precision (`Float64` for `Float64` or integer input,
#   `BigFloat` for `BigFloat` input).  Arguments are perturbed in their own types, so a
#   `Rotor{T}` argument stays an unnormalized `Rotor` and a `QuatVec` stays a `QuatVec`.
#
# Losses and comparisons
#
#   sc(y)              A fixed weighted sum of all the real coordinates of `y`: a real
#                      number, a complex number, a quaternion with real or complex
#                      components, or an array or tuple of these.  It preserves the
#                      element type (`Float32`, `Dual`, `TrackedReal`, …).
#   realcoords(y)      The real coordinates of `y` as a `Vector`: 1 for a real number, 2
#                      for a complex number, 4 (or 8 if complex) for a `Quaternion` or
#                      `Rotor`, 3 (or 6) for a `QuatVec`, concatenated for arrays and tuples.
#                      Complex components are listed as all real parts, then all imaginary
#                      parts.
#   relerr(a, b, floor=1)
#                      `norm(a - b) / max(norm(b), floor)`, comparing numbers, quaternions,
#                      and arrays or tuples of them by their real coordinates.  Thunks are
#                      unthunked; `NoTangent()`, `ZeroTangent()`, and `nothing` count as
#                      zeros of the other side's shape; a `Tangent` or `NamedTuple` is
#                      compared by its field values.  With the default floor, the result is
#                      an absolute error wherever the reference is smaller than 1, which is
#                      true of many derivatives at :nearidentity and :nearminus1 (of order
#                      1e-9 to 1e-7).  To detect a relative error there, pass a small floor
#                      or 0.
#
# Builders, from consecutive entries of a real vector `x`, starting at index `i` (default 1)
#
#   Q(x, i)    quaternion(x[i], x[i+1], x[i+2], x[i+3])
#   R(x, i)    rotor(x[i], …, x[i+3]), which normalizes
#   Ru(x, i)   Rotor{eltype(x)}(x[i], …, x[i+3]), which stores without normalizing
#   V(x, i)    quatvec(x[i+1], x[i+2], x[i+3]); `x[i]` is ignored, so that `V` takes the
#              same four-entry points as the others (its derivative with respect to `x[i]`
#              is zero)
#   QC(x, i)   Quaternion{Complex}: real parts x[i:i+3], imaginary parts x[i+4:i+7]
#   L(x, i)    Lorentz rotor Boost(x[i], NHAT) * Lorentz(rotor(1, x[i+1], x[i+2], x[i+3])),
#              with a rapidity-form boost
#
# Constants
#
#   P0, R0, V0, NHAT   A fixed generic `Quaternion`, `Rotor`, `QuatVec`, and unit 3-vector
#                      (a `Vector{Float64}`), for the other argument of binary functions.
#
# Points: four-entry vectors, used as `x` for the builders above
#
#   POINTS               All points, as `name => x` pairs (a `Vector{Pair{Symbol,Vector{Float64}}}`):
#                          :generic        [1.2, -0.7, 0.5, 0.3]  (norm ≈ 1.5, so also off-sphere)
#                          :identity       [1, 0, 0, 0]
#                          :nearidentity   [1, 1e-9, -2e-9, 3e-9]
#                          :purevector     [0, 0.6, -0.48, 0.64]  (unit norm)
#                          :nearminus1     [-1, 1e-7, 2e-7, -1e-7]
#                          :minus1         [-1, 0, 0, 0]
#                          :offsphere      [0.5, 0.3, -0.2, 0.4]  (norm ≈ 0.73, w > 0)
#                          :offspherenegw  [-0.8, 0.6, 0.5, -0.7] (norm ≈ 1.32, w < 0)
#   points(kind=:all; exclude=())         The pairs of one kind, in the order above, minus
#   points(T, kind=:all; exclude=())      the names in `exclude`, optionally converted to
#                                         element type `T`.  The kinds are
#                          :all        every point
#                          :nocut      not on the negative real axis (excludes :minus1),
#                                      for `log`, `sqrt`, and non-integer powers
#                          :nonreal    nonzero vector part (excludes :identity and :minus1),
#                                      for `absvec`, `angle`, and the like
#                          :unit       unit norm (:identity, :nearidentity, :purevector,
#                                      :nearminus1, and :minus1, within 3e-14), where the
#                                      values (not the derivatives) of `R(x)` and `Ru(x)`
#                                      agree
#                          :offsphere  norm different from 1 (:generic, :offsphere, and
#                                      :offspherenegw), for `Ru`
#   point(name)                           The vector of one point.
#
#   Both `points` and `point` return fresh copies of the vectors, so that a test item may
#   modify its point in place without changing the point set for later items.
#
# Random tangents for ChainRulesTestUtils
#
#   This module adds one method, `ChainRulesTestUtils.rand_tangent(rng, ::AbstractQuaternion)`,
#   which returns a `Quaternion` for a `Quaternion` or a `Rotor` (never a `Rotor`), a
#   `QuatVec` for a `QuatVec`, and `NoTangent()` for `Bool` components (`𝐢`, `𝐣`, `𝐤`).  The
#   components are drawn with ChainRulesTestUtils' own methods (so complex components get
#   complex tangents), after conversion to a floating-point type, so that integer
#   components get `Float64` tangents, as `ProjectTo` gives them.  The method is defined
#   only here, once, because redefining `rand_tangent` in a running session has crashed it
#   with SIGILL; no other test file may define it.  Note that ChainRulesTestUtils' finite
#   differences cannot perturb an integer-valued primal (FiniteDifferences rebuilds it with
#   integer components and throws an `InexactError`), so rules with integer quaternion
#   arguments are checked against `bigvjp` and `bigjvp` instead.

@testmodule ADTestUtils begin
    using Random: AbstractRNG
    using LinearAlgebra: norm, normalize, dot
    using Quaternionic
    using ChainRulesCore: NoTangent, AbstractZero, AbstractThunk, Tangent, backing, unthunk
    import ChainRulesTestUtils

    export bigfd, bigfdH, bigjvp, bigvjp, sc, realcoords, relerr
    export Q, R, Ru, V, QC, L, P0, R0, V0, NHAT
    export POINTS, points, point

    # ── Real coordinates ─────────────────────────────────────────────────────────────────

    # The finite-difference references treat every argument and output as a vector of
    # real coordinates.  The layout is fixed by a template value `t` of the primal type,
    # and the coordinates are read from a value `x` with the same layout: the primal itself,
    # a perturbed primal, or a tangent of it.  A `Rotor` and a `Quaternion` both have four
    # real coordinates (eight if complex), and a `QuatVec` has three (six if complex),
    # because its scalar part is always zero.

    """
        ncoords(t)

    Return the number of real coordinates of a value with the type and shape of `t`.
    """
    ncoords(::Real) = 1
    ncoords(::Complex) = 2
    ncoords(::AbstractQuaternion{<:Real}) = 4
    ncoords(::QuatVec{<:Real}) = 3
    ncoords(::AbstractQuaternion{<:Complex}) = 8
    ncoords(::QuatVec{<:Complex}) = 6
    ncoords(t::Union{AbstractArray,Tuple}) = sum(ncoords, t; init=0)

    """
        coords(t, x)

    Return the real coordinates of `x` as a `Vector`, in the layout of the template `t`.
    Here, `x` is a value of the same type as `t`, or a tangent of it; a `Quaternion` tangent
    of a `QuatVec` contributes only its vector part.  A zero tangent (`NoTangent()`,
    `ZeroTangent()`, or `nothing`) gives zero coordinates.
    """
    coords(::Real, x::Real) = [x]
    coords(::Complex, x::Number) = [real(x), imag(x)]
    coords(::AbstractQuaternion{<:Real}, x::AbstractQuaternion) = [x[1], x[2], x[3], x[4]]
    coords(::QuatVec{<:Real}, x::AbstractQuaternion) = [x[2], x[3], x[4]]
    function coords(::AbstractQuaternion{<:Complex}, x::AbstractQuaternion)
        c = components(x)
        [real(c[1]), real(c[2]), real(c[3]), real(c[4]),
         imag(c[1]), imag(c[2]), imag(c[3]), imag(c[4])]
    end
    function coords(::QuatVec{<:Complex}, x::AbstractQuaternion)
        c = components(x)
        [real(c[2]), real(c[3]), real(c[4]), imag(c[2]), imag(c[3]), imag(c[4])]
    end
    coords(t::Union{AbstractArray,Tuple}, x) =
        reduce(vcat, [coords(tk, xk) for (tk, xk) in zip(t, x)]; init=Float64[])
    coords(t, ::Union{AbstractZero,Nothing}) = zeros(ncoords(t))
    coords(t::Union{AbstractArray,Tuple}, ::Union{AbstractZero,Nothing}) = zeros(ncoords(t))

    """
        realcoords(y)

    Return the real coordinates of `y` as a `Vector`; see `coords`.
    """
    realcoords(y) = coords(y, y)

    # The type that stores a perturbed primal with components of type `S`, and the type of
    # a tangent with components of type `S`.  A perturbed `Rotor` is stored without
    # normalization, and the tangent of a `Rotor` is a `Quaternion`.
    primaltype(::Quaternion, ::Type{S}) where {S} = Quaternion{S}
    primaltype(::Rotor, ::Type{S}) where {S} = Rotor{S}
    primaltype(::QuatVec, ::Type{S}) where {S} = QuatVec{S}
    tangenttype(::AbstractQuaternion, ::Type{S}) where {S} = Quaternion{S}
    tangenttype(::QuatVec, ::Type{S}) where {S} = QuatVec{S}

    """
        build(kind, t, c)

    Build a quaternion in the layout of the template `t` from its real coordinates `c`,
    with the type `kind(t, S)`, where `S` is the component type and `kind` is `primaltype`
    or `tangenttype`.
    """
    build(kind, t::AbstractQuaternion{<:Real}, c) = kind(t, eltype(c))(c[1], c[2], c[3], c[4])
    build(kind, t::QuatVec{<:Real}, c) = kind(t, eltype(c))(c[1], c[2], c[3])
    build(kind, t::AbstractQuaternion{<:Complex}, c) = kind(t, Complex{eltype(c)})(
        complex(c[1], c[5]), complex(c[2], c[6]), complex(c[3], c[7]), complex(c[4], c[8])
    )
    build(kind, t::QuatVec{<:Complex}, c) = kind(t, Complex{eltype(c)})(
        complex(c[1], c[4]), complex(c[2], c[5]), complex(c[3], c[6])
    )

    # Split the coordinates `c` of an array or tuple into one vector per entry of the
    # template `t`.
    function splitcoords(t, c)
        ts = vec(collect(t))
        ends = cumsum([ncoords(tk) for tk in ts])
        [c[e-ncoords(tk)+1:e] for (tk, e) in zip(ts, ends)]
    end

    """
        fromcoords(t, c)

    Return the value of the same type as the template `t` (but with components of the type
    of `c`) whose real coordinates are `c`.  An array is rebuilt as an `Array` of the same
    shape.
    """
    fromcoords(::Real, c) = c[1]
    fromcoords(::Complex, c) = complex(c[1], c[2])
    fromcoords(t::AbstractQuaternion, c) = build(primaltype, t, c)
    fromcoords(t::AbstractArray, c) = reshape(map(fromcoords, vec(collect(t)), splitcoords(t, c)), size(t))
    fromcoords(t::Tuple, c) = Tuple(map(fromcoords, collect(Any, t), splitcoords(t, c)))

    """
        tangentfromcoords(t, c)

    Return the tangent of the template `t` whose real coordinates are `c`, in the tangent
    type of the ChainRules conventions: a `Quaternion` for a `Quaternion` or `Rotor`, and a
    `QuatVec` for a `QuatVec`.
    """
    tangentfromcoords(::Real, c) = c[1]
    tangentfromcoords(::Complex, c) = complex(c[1], c[2])
    tangentfromcoords(t::AbstractQuaternion, c) = build(tangenttype, t, c)
    tangentfromcoords(t::AbstractArray, c) =
        reshape(map(tangentfromcoords, vec(collect(t)), splitcoords(t, c)), size(t))
    tangentfromcoords(t::Tuple, c) = Tuple(map(tangentfromcoords, collect(Any, t), splitcoords(t, c)))

    """
        isdifferentiable(x)

    Return `true` if the argument `x` is perturbed by the finite-difference references.
    Integers, `Bool`s, quaternions with `Bool` components (the constants `𝐢`, `𝐣`, and
    `𝐤`), and non-numeric values are not, in keeping with ChainRulesTestUtils, which gives
    them `NoTangent()`.  Quaternions with integer components are, because their tangents
    are floating-point quaternions.
    """
    isdifferentiable(::Any) = false
    isdifferentiable(::Real) = true
    isdifferentiable(::Integer) = false
    isdifferentiable(::Complex) = true
    isdifferentiable(::Complex{Bool}) = false
    isdifferentiable(::AbstractQuaternion) = true
    isdifferentiable(::AbstractQuaternion{Bool}) = false
    isdifferentiable(a::AbstractArray) = !isempty(a) && all(isdifferentiable, a)

    """
        realtype(x)

    Return the floating-point type of the real coordinates of `x` at its own precision.
    """
    realtype(x::Real) = float(typeof(x))
    realtype(x::Complex) = float(real(typeof(x)))
    realtype(::AbstractQuaternion{T}) where {T} = realtype(zero(T))
    realtype(a::Union{AbstractArray,Tuple}) =
        isempty(a) ? Float64 : mapreduce(realtype, promote_type, a)

    """
        tobig(x)

    Return `x` with its real coordinates converted to `BigFloat`, at the current precision,
    if `x` is differentiable; otherwise return `x` unchanged.
    """
    tobig(x) = isdifferentiable(x) ? fromcoords(x, BigFloat.(coords(x, x))) : x

    # ── Finite-difference references ─────────────────────────────────────────────────────

    const PRECISION = 256

    # The first-derivative stencil is the fourth-order central difference with the step
    # h = 2⁻¹⁰⁰ ≈ 8e-31, an exact power of two.  Its truncation error is of order
    # h⁴ f⁽⁵⁾ ≈ 4e-121 f⁽⁵⁾, and its rounding error is of order 2⁻²⁵⁶ |f| / h ≈ 1e-47 |f|.
    # The Hessian stencil is the second-order central difference with h = 2⁻⁸⁰ ≈ 8e-25.  Its
    # truncation error is of order h² f⁗ ≈ 7e-49 f⁗, and its rounding error is of order
    # 2⁻²⁵⁶ |f| / h² ≈ 1e-29 |f|.
    const JVPSTEPEXPONENT = -100
    const HESSIANSTEPEXPONENT = -80

    """
        coordjacobian(f, y, args, i)

    Return the Jacobian of the real coordinates of `f(args...)` with respect to the real
    coordinates of `args[i]`, as a `Matrix{BigFloat}`, where `y = f(args...)` has already
    been computed at the original precision.  It must be called at 256-bit precision.
    """
    function coordjacobian(f, y, args, i)
        bigargs = map(tobig, args)
        c = coords(args[i], bigargs[i])
        h = BigFloat(2)^JVPSTEPEXPONENT
        columns = map(eachindex(c)) do k
            function F(s)
                ck = copy(c)
                ck[k] += s
                perturbed = ntuple(j -> j == i ? fromcoords(args[i], ck) : bigargs[j], length(args))
                coords(y, f(perturbed...))
            end
            (8 .* (F(h) .- F(-h)) .- (F(2h) .- F(-2h))) ./ (12h)
        end
        isempty(columns) ? zeros(BigFloat, ncoords(y), 0) : reduce(hcat, columns)
    end

    """
        bigjvp(f, args::Tuple, ẋs::Tuple)
        bigjvp(f, x, ẋ)

    Return the directional derivative of `f(args...)` along the tangents `ẋs`, computed by
    central differences in 256-bit arithmetic, as a tangent of the output: a `Quaternion`
    for a `Quaternion` or `Rotor` output, a `QuatVec` for a `QuatVec` output, a real or
    complex number for a scalar output, and an array or tuple of these for an array or
    tuple output.  The tangent of a `Rotor` argument is a `Quaternion`, and the tangents of
    integer, `Bool`, and non-numeric arguments are ignored.  The second form is for a single
    argument `x` that is not a `Tuple`.
    """
    function bigjvp(f, args::Tuple, ẋs::Tuple)
        y = f(args...)
        T = realtype(y)
        D = setprecision(BigFloat, PRECISION) do
            D = zeros(BigFloat, ncoords(y))
            for i in eachindex(args)
                isdifferentiable(args[i]) || continue
                ċ = BigFloat.(coords(args[i], ẋs[i]))
                iszero(ċ) && continue
                D .+= coordjacobian(f, y, args, i) * ċ
            end
            T.(D)
        end
        tangentfromcoords(y, D)
    end
    bigjvp(f, x, ẋ) = bigjvp(f, (x,), (ẋ,))

    """
        bigvjp(f, args::Tuple, ȳ)
        bigvjp(f, x, ȳ)

    Return the cotangents of the arguments of `f(args...)` pulled back from the cotangent
    `ȳ` of its output, computed by central differences in 256-bit arithmetic, as a `Tuple`
    with one entry per argument.  This is the conjugate transpose of the Jacobian applied to
    `ȳ`, as in ChainRules: each real coordinate of a cotangent is the real inner product of
    `ȳ` with the derivative of the output along that coordinate.  The cotangent of a
    `Rotor` argument is a `Quaternion`, that of a `QuatVec` is a `QuatVec`, and that of an
    integer, `Bool`, or non-numeric argument is `NoTangent()`.  The second form is for a
    single argument `x` that is not a `Tuple`, and returns its cotangent alone.
    """
    function bigvjp(f, args::Tuple, ȳ)
        y = f(args...)
        cȳ = coords(y, ȳ)
        ntuple(length(args)) do i
            a = args[i]
            isdifferentiable(a) || return NoTangent()
            T = realtype(a)
            g = setprecision(BigFloat, PRECISION) do
                T.(transpose(coordjacobian(f, y, args, i)) * BigFloat.(cȳ))
            end
            tangentfromcoords(a, g)
        end
    end
    bigvjp(f, x, ȳ) = bigvjp(f, (x,), ȳ)[1]

    """
        bigfd(f, x::AbstractVector{<:Real})
        bigfd(f, t::Real)

    Return the derivative of `f`, computed by central differences in 256-bit arithmetic.
    For a vector `x`, this is the gradient (a `Vector`) if `f(x)` is a real number, and
    otherwise the Jacobian (a `Matrix` with one row per entry of `realcoords(f(x))`).  For a
    real number `t`, this is the derivative as a tangent of `f(t)`, as for `bigjvp`.  The
    result has the floating-point type of `f(x)` at its original precision.
    """
    function bigfd(f, x::AbstractVector{<:Real})
        xf = float.(collect(x))
        y = f(xf)
        T = realtype(y)
        J = setprecision(BigFloat, PRECISION) do
            T.(coordjacobian(f, y, (xf,), 1))
        end
        y isa Real ? vec(J) : J
    end
    bigfd(f, t::Real) = bigjvp(f, (float(t),), (one(float(t)),))

    """
        bigfdH(f, x::AbstractVector{<:Real})

    Return the Hessian of the real function `f` at `x`, computed by central second
    differences in 256-bit arithmetic, as a symmetric `Matrix` with the floating-point type
    of `f(x)` at its original precision.
    """
    function bigfdH(f, x::AbstractVector{<:Real})
        xf = float.(collect(x))
        y = f(xf)
        y isa Real || throw(ArgumentError("bigfdH needs a real-valued function; got $(typeof(y))"))
        T = realtype(y)
        n = length(xf)
        setprecision(BigFloat, PRECISION) do
            xb = BigFloat.(xf)
            h = BigFloat(2)^HESSIANSTEPEXPONENT
            e(i) = (v = zeros(BigFloat, n); v[i] = h; v)
            H = zeros(BigFloat, n, n)
            for i in 1:n, j in i:n
                H[i, j] = H[j, i] = (
                    f(xb .+ e(i) .+ e(j)) - f(xb .+ e(i) .- e(j))
                    - f(xb .- e(i) .+ e(j)) + f(xb .- e(i) .- e(j))
                ) / (4h^2)
            end
            T.(H)
        end
    end

    # ── Losses and comparisons ───────────────────────────────────────────────────────────

    # The weights are ratios of integers, so that the loss has the element type of its
    # argument: a `Float32` stays a `Float32`, and a `Dual` stays a `Dual`.  The weights of
    # the imaginary parts differ from those of the real parts, so that a derivative that
    # confuses the two is detected.

    """
        sc(y)

    Return a fixed weighted sum of the real coordinates of `y`, which may be a real or
    complex number, a quaternion with real or complex components, or an array or tuple of
    these.  This is the scalar loss of the AD tests.
    """
    sc(x::Real) = x
    sc(z::Complex) = real(z) + 37 * imag(z) / 100
    function sc(q::AbstractQuaternion{<:Real})
        c = components(q)
        (3c[1] - 11c[2] + 7c[3] + 5c[4]) / 10
    end
    function sc(q::AbstractQuaternion{<:Complex})
        c = components(q)
        (3real(c[1]) - 11real(c[2]) + 7real(c[3]) + 5real(c[4])
         + 5imag(c[1]) + 7imag(c[2]) - 11imag(c[3]) + 3imag(c[4])) / 10
    end
    # The entries of an array are weighted by 1, 9/8, 10/8, …, in linear-index order.  The
    # sum is an explicit loop, which every backend differentiates without any broadcast
    # machinery.
    function sc(a::AbstractArray)
        s = sc(a[1])
        for k in 2:length(a)
            s += (7 + k) * sc(a[k]) / 8
        end
        s
    end
    # The entries of a tuple are weighted in the same way.
    sc(t::Tuple) = sc(t, 1)
    sc(::Tuple{}, ::Integer) = false
    sc(t::Tuple, k) = (7 + k) * sc(first(t)) / 8 + sc(Base.tail(t), k + 1)

    """
        relerr(a, b, floor=1)

    Return the error of `a` relative to the reference `b`, as `norm(a - b) / max(norm(b),
    floor)`, so that with the default floor a vanishing reference does not demand an exact
    zero.  Values are compared by their real coordinates, as given by `cmpvec`, and arrays
    and tuples are compared entry by entry.  Thunks are unthunked first, and a zero tangent
    (`NoTangent()`, `ZeroTangent()`, or `nothing`) is compared as zeros of the shape of the
    other side.  A difference in shape throws a `DimensionMismatch`.  If the denominator is
    zero, the result is 0 when `a` equals `b` and `Inf` otherwise.
    """
    function relerr(a, b, floor=1)
        va, vb = pairedcoords(a, b)
        e, d = norm(va .- vb), max(norm(vb), floor)
        iszero(d) && return iszero(e) ? zero(e) : oftype(e, Inf)
        e / d
    end

    """
        pairedcoords(a, b)

    Return the real coordinates of `a` and of `b` as two `Vector`s of equal length, for
    `relerr`.  Arrays and tuples are paired entry by entry, and a zero tangent on one side
    is replaced by zeros of the length of the other side.
    """
    function pairedcoords(a, b)
        a, b = unthunk(a), unthunk(b)
        if a isa Union{AbstractArray,Tuple} && b isa Union{AbstractArray,Tuple}
            length(a) == length(b) || throw(DimensionMismatch(
                "cannot compare containers of lengths $(length(a)) and $(length(b))"
            ))
            pairs = [pairedcoords(ak, bk) for (ak, bk) in zip(a, b)]
            return (reduce(vcat, first.(pairs); init=Float64[]),
                    reduce(vcat, last.(pairs); init=Float64[]))
        end
        va, vb = cmpvec(a), cmpvec(b)
        isempty(va) && a isa Union{AbstractZero,Nothing} && (va = zero(vb))
        isempty(vb) && b isa Union{AbstractZero,Nothing} && (vb = zero(va))
        length(va) == length(vb) || throw(DimensionMismatch(
            "cannot compare $(typeof(a)) with $(length(va)) coordinates " *
            "and $(typeof(b)) with $(length(vb)) coordinates"
        ))
        va, vb
    end

    # The real coordinates of a value, for comparisons: a `Rotor` or `QuatVec` gradient is
    # compared by all four of its components, like a `Quaternion`.  A structural tangent
    # (a ChainRulesCore `Tangent` or a `NamedTuple`, such as Zygote may return for a
    # quaternion) is compared by the coordinates of its field values, so that the tangent
    # of a quaternion's `components` field is compared like the quaternion itself.
    cmpvec(x::Real) = [x]
    cmpvec(x::Complex) = [real(x), imag(x)]
    cmpvec(q::AbstractQuaternion) = coords(quaternion(q), q)
    cmpvec(a::Union{AbstractArray,Tuple}) = reduce(vcat, [cmpvec(x) for x in a]; init=Float64[])
    cmpvec(::Union{AbstractZero,Nothing}) = Float64[]
    cmpvec(t::Tangent) = cmpvec(Tuple(values(backing(t))))
    cmpvec(t::NamedTuple) = cmpvec(Tuple(t))
    cmpvec(x::AbstractThunk) = cmpvec(unthunk(x))

    # ── Builders and constants ───────────────────────────────────────────────────────────

    """
        Q(x, i=1)
        R(x, i=1)
        Ru(x, i=1)
        V(x, i=1)
        QC(x, i=1)
        L(x, i=1)

    Build a `Quaternion`, a normalized `Rotor`, an unnormalized `Rotor{eltype(x)}`, a
    `QuatVec` (ignoring `x[i]`), a `Quaternion` with complex components (real parts
    `x[i:i+3]` and imaginary parts `x[i+4:i+7]`), or a Lorentz rotor (a rapidity-form boost
    of rapidity `x[i]` along `NHAT` times the rotor of `x[i+1:i+3]`) from consecutive
    entries of `x`.
    """
    Q(x, i=1) = quaternion(x[i], x[i+1], x[i+2], x[i+3])
    R(x, i=1) = rotor(x[i], x[i+1], x[i+2], x[i+3])
    Ru(x, i=1) = Rotor{eltype(x)}(x[i], x[i+1], x[i+2], x[i+3])
    V(x, i=1) = quatvec(x[i+1], x[i+2], x[i+3])
    QC(x, i=1) = quaternion(
        complex(x[i], x[i+4]), complex(x[i+1], x[i+5]),
        complex(x[i+2], x[i+6]), complex(x[i+3], x[i+7])
    )
    L(x, i=1) = Boost(x[i], NHAT) * Lorentz(rotor(one(eltype(x)), x[i+1], x[i+2], x[i+3]))

    const P0 = quaternion(-0.5, 0.2, 0.9, -1.3)
    const R0 = rotor(0.4, -0.3, 0.8, 0.2)
    const V0 = quatvec(0.3, -0.6, 0.2)
    const NHAT = normalize([0.3, -0.5, 0.81])

    # ── Points ───────────────────────────────────────────────────────────────────────────

    const POINTS = Pair{Symbol,Vector{Float64}}[
        :generic => [1.2, -0.7, 0.5, 0.3],
        :identity => [1.0, 0.0, 0.0, 0.0],
        :nearidentity => [1.0, 1e-9, -2e-9, 3e-9],
        :purevector => [0.0, 0.6, -0.48, 0.64],
        :nearminus1 => [-1.0, 1e-7, 2e-7, -1e-7],
        :minus1 => [-1.0, 0.0, 0.0, 0.0],
        :offsphere => [0.5, 0.3, -0.2, 0.4],
        :offspherenegw => [-0.8, 0.6, 0.5, -0.7],
    ]

    const POINTKINDS = Dict{Symbol,Vector{Symbol}}(
        :all => first.(POINTS),
        :nocut => filter(!=(:minus1), first.(POINTS)),
        :nonreal => filter(n -> n ∉ (:identity, :minus1), first.(POINTS)),
        :unit => [:identity, :nearidentity, :purevector, :nearminus1, :minus1],
        :offsphere => [:generic, :offsphere, :offspherenegw],
    )

    """
        points(kind=:all; exclude=())
        points(T, kind=:all; exclude=())

    Return the points of the given kind as `name => x` pairs, in the order of `POINTS`,
    without those whose names are in `exclude`, and with entries converted to `T` if it is
    given.  The kinds are `:all`, `:nocut` (not on the negative real axis), `:nonreal`
    (nonzero vector part), `:unit` (unit norm), and `:offsphere` (norm different from 1).
    """
    function points(kind::Symbol=:all; exclude=())
        haskey(POINTKINDS, kind) || throw(ArgumentError("unknown kind of point: $kind"))
        for name in exclude
            any(p -> p.first === name, POINTS) || throw(ArgumentError("unknown point: $name"))
        end
        names = POINTKINDS[kind]
        [p.first => copy(p.second) for p in POINTS if p.first ∈ names && p.first ∉ exclude]
    end
    points(::Type{T}, kind::Symbol=:all; exclude=()) where {T} =
        [p.first => T.(p.second) for p in points(kind; exclude=exclude)]

    """
        point(name)

    Return the vector of the point called `name` in `POINTS`.
    """
    function point(name::Symbol)
        i = findfirst(p -> p.first === name, POINTS)
        i === nothing && throw(ArgumentError("unknown point: $name"))
        copy(POINTS[i].second)
    end

    # ── Random tangents for ChainRulesTestUtils ──────────────────────────────────────────

    # ChainRulesTestUtils' generic method for `Number`s calls `randn(rng, T)`, which gives a
    # normalized `Rotor` for a `Rotor` and has no method for complex quaternions.  This
    # method replaces it for every quaternion type.  It must be defined only once per
    # session; see the comment at the top of this file.
    function ChainRulesTestUtils.rand_tangent(rng::AbstractRNG, q::AbstractQuaternion)
        isdifferentiable(q) || return NoTangent()
        c = components(q)
        t(k) = ChainRulesTestUtils.rand_tangent(rng, float(c[k]))
        q isa QuatVec ? quatvec(t(2), t(3), t(4)) : quaternion(t(1), t(2), t(3), t(4))
    end
end
