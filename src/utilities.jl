# Un-normalized sinc function
# sincu(x::Number) = _sincu(float(x))
#
# The generic method is used by types such as `ForwardDiff.Dual`, whose derivatives of
# `sin(x)/x` lose accuracy to cancellation for small `x`: the error of the `k`th derivative
# is roughly `eps/x^k`.  So the series is used whenever the *value* of `x` is below the
# threshold given by `sincu_tol`, rather than only at zero.  The series runs through `x^30`,
# so that, for `Float64` and `Float32`, it can be used up to the threshold 3, beyond which
# the closed form is accurate.  Measured against a high-precision reference, with the error
# of the `k`th derivative taken relative to the larger of its magnitude and `1/(k+1)`, the
# value and its first four derivatives are then accurate to a few ulps for every value of
# `x` below 4, including both sides of the threshold, and the fifth and sixth derivatives to
# about 10 and 70 ulps near the threshold, where their errors are largest.  Types whose
# value is not a floating-point number, such as symbolic types, take the series only when
# the value is exactly zero.
_sincu(x::Number) = _sincu(x, value(x))
_sincu(x, v) = iszerovalue(x) ? 1 - x^2/6 * (1 - x^2/20) : isinf(x) ? zero(x) : sin(x)/x
function _sincu(x, v::Union{R,Complex{R}}) where {R<:AbstractFloat}
    if abs(v) < sincu_tol(R)
        sincu_series(x, R)
    else
        isinf(x) ? zero(x) : csin(x)/x
    end
end
# The Taylor series of `sin(x)/x` through `x^30`, whose coefficients are `(-1)^n/(2n+1)!`
# for `n = 0, …, 15`.  For `Float64` and `Float32`, the coefficients are stored.  For other
# floating-point types, the series is evaluated in nested form, dividing by the exact
# integers `(2n)(2n+1)`, so that no coefficients need to be computed in that type.
sincu_series(x, ::Type{Float64}) = evalpoly(x * x, sincu_coefficients(Float64))
sincu_series(x, ::Type{Float32}) = evalpoly(x * x, sincu_coefficients(Float32))
function sincu_series(x, ::Type{R}) where {R<:AbstractFloat}
    x² = x * x
    s = one(x²)
    for n ∈ 15:-1:1
        s = 1 - x² * s / ((2n) * (2n + 1))
    end
    s
end
sincu_coefficients(::Type{Float64}) = (
    1.0, -0.16666666666666666, 0.008333333333333333, -0.0001984126984126984,
    2.7557319223985893e-6, -2.505210838544172e-8, 1.6059043836821613e-10,
    -7.647163731819816e-13, 2.8114572543455206e-15, -8.22063524662433e-18,
    1.9572941063391263e-20, -3.868170170630684e-23, 6.446950284384474e-26,
    -9.183689863795546e-29, 1.1309962886447716e-31, -1.216125041553518e-34
)
sincu_coefficients(::Type{Float32}) = map(Float32, sincu_coefficients(Float64))
# For a general floating-point type, this threshold makes the first omitted term of the
# series, `x^32/33!`, about `2^-32 * eps(R)`.  It is capped at 3, where the closed form is
# accurate, as described above.  The constant is `log(factorial(33))`.  The threshold is
# computed in `Float64`, which is precise enough and avoids the cost of evaluating `exp` and
# `log` in, for example, `BigFloat` on every call; `eps(R)` is bounded below by `floatmin`,
# so that the threshold stays positive even for precisions beyond the range of `Float64`.
sincu_tol(::Type{R}) where {R<:AbstractFloat} =
    min(3.0, exp((log(max(Float64(eps(R)), floatmin(Float64))) + 85.05446701758152) / 32) / 2)
sincu_tol(::Type{Float64}) = 3.0
sincu_tol(::Type{Float32}) = 3.0f0
# _sincu(x::Number) = ifelse(
#     iszerovalue(x),
#     one(x),
#     ifelse(
#         isinf(x),
#         zero(x),
#         sin(x) / x
#     )
# )
@inline _sincu(x::Union{Float64,ComplexF64}) =
    abs(x) < 0.0031 ? evalpoly(x^2, (1.0, -0.16666666666666666, 0.008333333333333333)) :
    isinf(x) ? zero(x) : csin(x)/x
@inline _sincu(x::Union{Float32,ComplexF32}) =
    abs(x) < 0.1571f0 ? evalpoly(x^2, (1.0f0, -0.16666667f0, 0.008333334f0)) :
    isinf(x) ? zero(x) : csin(x)/x
_sincu(x::Float16) = Float16(_sincu(Float32(x)))
# _sincu(x::ComplexF16) = ComplexF16(_sincu(ComplexF32(x)))


# """
#     invsinc(x)

# Inverse of the un-normalized sinc function, defined as `x/sin(x)` for `x != 0` and `1` for
# `x == 0`.

# Note that this is only implemented accurately for ``|x| \\lesssim 1``.  When ``|x|``
# approaches ``n\\pi``, the function value blows up.  This could be avoided as something
# like ``invsinc(mod(x, π)) * x / mod(x, π)``, but in the interest of efficiency this is not
# implemented.

# This function is defined to be continuous at `x == 0`, and is implemented using a Taylor
# series for small `x`.
# """
@inline invsinc(x::T) where {T<:Union{Real,Complex{<:Real}}} =
    abs(x) < invsinc_tol(T) ? evalpoly(x^2,
        T.((1, 1//6, 7//360, 31//15120, 127//604800, 73//3421440, 1414477//653837184000))
    ) :
    isinf(x) ? x : x/sin(x)
invsinc(x::Float16) = Float16(invsinc(Float32(x)))
#invsinc(x::ComplexF16) = ComplexF16(invsinc(ComplexF32(x)))
invsinc_tol(::Type{T}) where {T} = sqrt(sqrt(sqrt(eps(T))))
invsinc_tol(::Type{Complex{T}}) where {T<:Real} = invsinc_tol(T)
@inline invsinc(x::Complex{T}) where {T<:AbstractFloat} =
    abs(x) < invsinc_tol(T) ? evalpoly(x^2,
        (one(x), one(x)/6, 7one(x)/360, 31one(x)/15120, 127one(x)/604800, 73one(x)/3421440, 1414477one(x)/653837184000)
    ) : x/sin(x)

# # Derivative of un-normalized sinc function
# function _coscu(x::Number)
#     # naive coscu formula is susceptible to catastrophic
#     # cancellation error near x=0, so we use the Taylor series
#     # for small enough |x|.
#     if abs(x) < 1.57
#         # generic Taylor series: ∑ (-1)^n (x)^{2n-1}/a(n) where
#         # a(n) = (1+2n)*(2n-1)! (= OEIS A174549)
#         s = (term = -x)/3
#         x² = term^2
#         ε = eps(abs(term)) # error threshold to stop sum
#         n = 1
#         while true
#             n += 1
#             term *= x²/((1-2n)*(2n-2))
#             s += (δs = term/(1+2n))
#             abs(δs) ≤ ε && break
#         end
#         return s
#     else
#         return isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
#     end
# end
# # hard-code Float64/Float32 Taylor series, with coefficients
# #  Float64.([(-1)^n/((2n+1)*factorial(2n-1)) for n = 1:6])
# _coscu(x::Union{Float64,ComplexF64}) =
#     abs(x) < 0.44 ? x*evalpoly(x^2, (-0.3333333333333333, 0.03333333333333333, -0.0011904761904761906, 2.2045855379188714e-5, -2.505210838544172e-7, 1.9270852604185937e-9)) :
#     isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
# _coscu(x::Union{Float32,ComplexF32}) =
#     abs(x) < 0.817f0 ? x*evalpoly(x^2, (-0.333333335f0, 0.033333335f0, -0.0011904762f0, 2.2045855f-5)) :
#     isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
# _coscu(x::Float16) = Float16(_coscu(Float32(x)))
# _coscu(x::ComplexF16) = ComplexF16(_coscu(ComplexF32(x)))

# Derivative of un-normalized sinc function divided by `x`
function _cossu(x::Number)
    if abs(x) < 0.5
        # generic Taylor series: ∑ (-1)^n (x)^{2n-2}/a(n) where
        # a(n) = (1+2n)*(2n-1)! (= OEIS A174549)
        s = (term = -one(x)) / 3
        x² = x^2
        ε = eps(abs(x)) # error threshold to stop sum
        n = 1
        while true
            n += 1
            term *= x²/((1-2n)*(2n-2))
            s += (δs = term/(1+2n))
            abs(δs) ≤ ε && break
        end
        return s
    else
        return isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^3))
    end
end
# hard-code Float64/Float32 Taylor series, with coefficients
#  Float64.([(-1)^n/((2n+1)*factorial(big(2n-1))) for n = 1:10]) The cutoffs are chosen to
# balance the truncation error of the series against the cancellation in the closed form, so
# that the relative error stays below about 2 ulps.
_cossu(x::Union{Float64,ComplexF64}) =
    abs(x) < 1.2 ? evalpoly(x^2, (
        -0.3333333333333333, 0.03333333333333333, -0.0011904761904761906,
        2.2045855379188714e-5, -2.505210838544172e-7, 1.9270852604185937e-9,
        -1.0706029224547743e-11, 4.498331606952833e-14, -1.4797143443923793e-16,
        3.9145882126782523e-19
    )) :
    isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^3))
_cossu(x::Union{Float32,ComplexF32}) =
    abs(x) < 1.2f0 ? evalpoly(x^2, (
        -0.33333334f0, 0.033333335f0, -0.0011904762f0, 2.2045855f-5, -2.5052108f-7,
        1.9270852f-9
    )) :
    isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^3))
_cossu(x::Float16) = Float16(_cossu(Float32(x)))
# _cossu(x::ComplexF16) = ComplexF16(_cossu(ComplexF32(x)))
