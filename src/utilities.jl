# These functions are needed in lots of places for rotations.  I've made a PR to get
# them into julia 1.11 (https://github.com/JuliaLang/julia/pull/51551), but that may
# not happen, and anyways we need them for backwards compatibility.
if isdefined(Main, :sincu)
    const _sincu = Main.sincu
else
    # Un-normalized sinc function
    # sincu(x::Number) = _sincu(float(x))
    _sincu(x::Number) = iszero(x) ? one(x) : isinf(x) ? zero(x) : sin(x)/x
    @inline _sincu(x::Union{Float64,ComplexF64}) =
        abs(x) < 0.0031 ? evalpoly(x^2, (1.0, -0.16666666666666666, 0.008333333333333333)) :
        isinf(x) ? zero(x) : sin(x)/x
    @inline _sincu(x::Union{Float32,ComplexF32}) =
        abs(x) < 0.1571f0 ? evalpoly(x^2, (1.0f0, -0.16666667f0, 0.008333334f0)) :
        isinf(x) ? zero(x) : sin(x)/x
    _sincu(x::Float16) = Float16(_sincu(Float32(x)))
    # _sincu(x::ComplexF16) = ComplexF16(_sincu(ComplexF32(x)))
end

# """
#     invsinc(x)

# Inverse of the un-normalized sinc function, defined as `x/sin(x)` for `x != 0` and `1` for
# `x == 0`.

# Note that this is only implemented accurately for ``|x| \\lesssim 1``.  When ``|x|``
# approaches ``n\\pi``, the function value blows up.  This could be avoided as something like
# ``invsinc(mod(x, π)) * x / mod(x, π)``, but in the interest of efficiency this is not
# implemented.

# This function is defined to be continuous at `x == 0`, and is implemented using a Taylor
# series for small `x`.
# """
@inline invsinc(x::T) where {T<:Union{Real,Complex{Real}}} =
    abs(x) < invsinc_tol(T) ? evalpoly(x^2,
        T.((1, 1//6, 7//360, 31//15120, 127//604800, 73//3421440, 1414477//653837184000))
    ) :
    isinf(x) ? x : x/sin(x)
invsinc(x::Float16) = Float16(invsinc(Float32(x)))
#invsinc(x::ComplexF16) = ComplexF16(invsinc(ComplexF32(x)))
invsinc_tol(::Type{T}) where {T} = sqrt(sqrt(sqrt(eps(T))))

# if isdefined(Main, :coscu)
#     const _coscu = Main.coscu
# else
#     # Derivative of un-normalized sinc function
#     function _coscu(x::Number)
#         # naive coscu formula is susceptible to catastrophic
#         # cancellation error near x=0, so we use the Taylor series
#         # for small enough |x|.
#         if abs(x) < 1.57
#             # generic Taylor series: ∑ (-1)^n (x)^{2n-1}/a(n) where
#             # a(n) = (1+2n)*(2n-1)! (= OEIS A174549)
#             s = (term = -x)/3
#             x² = term^2
#             ε = eps(abs(term)) # error threshold to stop sum
#             n = 1
#             while true
#                 n += 1
#                 term *= x²/((1-2n)*(2n-2))
#                 s += (δs = term/(1+2n))
#                 abs(δs) ≤ ε && break
#             end
#             return s
#         else
#             return isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
#         end
#     end
#     # hard-code Float64/Float32 Taylor series, with coefficients
#     #  Float64.([(-1)^n/((2n+1)*factorial(2n-1)) for n = 1:6])
#     _coscu(x::Union{Float64,ComplexF64}) =
#         abs(x) < 0.44 ? x*evalpoly(x^2, (-0.3333333333333333, 0.03333333333333333, -0.0011904761904761906, 2.2045855379188714e-5, -2.505210838544172e-7, 1.9270852604185937e-9)) :
#         isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
#     _coscu(x::Union{Float32,ComplexF32}) =
#         abs(x) < 0.817f0 ? x*evalpoly(x^2, (-0.333333335f0, 0.033333335f0, -0.0011904762f0, 2.2045855f-5)) :
#         isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
#     _coscu(x::Float16) = Float16(_coscu(Float32(x)))
#     _coscu(x::ComplexF16) = ComplexF16(_coscu(ComplexF32(x)))
# end

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
#  Float64.([(-1)^n/((2n+1)*factorial(2n-1)) for n = 1:6])
_cossu(x::Union{Float64,ComplexF64}) =
    abs(x) < 0.44 ? evalpoly(x^2, (-0.3333333333333333, 0.03333333333333333, -0.0011904761904761906, 2.2045855379188714e-5, -2.505210838544172e-7, 1.9270852604185937e-9)) :
    isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^3))
_cossu(x::Union{Float32,ComplexF32}) =
    abs(x) < 0.817f0 ? evalpoly(x^2, (-0.333333335f0, 0.033333335f0, -0.0011904762f0, 2.2045855f-5)) :
    isinf(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^3))
_cossu(x::Float16) = Float16(_cossu(Float32(x)))
# _cossu(x::ComplexF16) = ComplexF16(_cossu(ComplexF32(x)))
