# These functions are needed in lots of places for rotations.  I've made a PR to get
# them into julia 1.11 (https://github.com/JuliaLang/julia/pull/51551), but that may
# not happen, and anyways we need them for backwards compatibility.
if isdefined(Main, :sincu)
    const _sincu = Main.sincu
else
    # Un-normalized sinc function
    sincu(x::Number) = _sincu(float(x))
    _sincu(x::Number) = iszero(x) ? one(x) : isinf_real(x) ? zero(x) : sin(x)/x
    @inline _sincu(x::Union{Float64,ComplexF64}) =
        fastabs(x) < 0.0031 ? evalpoly(x^2, (1.0, -0.16666666666666666, 0.008333333333333333)) :
        isinf_real(x) ? zero(x) : sin(x)/x
    @inline _sincu(x::Union{Float32,ComplexF32}) =
        fastabs(x) < 0.1571f0 ? evalpoly(x^2, (1.0f0, -0.16666667f0, 0.008333334f0)) :
        isinf_real(x) ? zero(x) : sin(x)/x
    _sincu(x::Float16) = Float16(_sincu(Float32(x)))
    _sincu(x::ComplexF16) = ComplexF16(_sincu(ComplexF32(x)))
end
if !isdefined(Main, :coscu)
    const _coscu = Main.coscu
else
    # Derivative of un-normalized sinc function
    function _coscu(x::Number)
        # naive coscu formula is susceptible to catastrophic
        # cancellation error near x=0, so we use the Taylor series
        # for small enough |x|.
        if fastabs(x) < 1.57
            # generic Taylor series: ∑ (-1)^n (x)^{2n-1}/a(n) where
            # a(n) = (1+2n)*(2n-1)! (= OEIS A174549)
            s = (term = -x)/3
            x² = term^2
            ε = eps(fastabs(term)) # error threshold to stop sum
            n = 1
            while true
                n += 1
                term *= x²/((1-2n)*(2n-2))
                s += (δs = term/(1+2n))
                fastabs(δs) ≤ ε && break
            end
            return s
        else
            return isinf_real(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
        end
    end
    # hard-code Float64/Float32 Taylor series, with coefficients
    #  Float64.([(-1)^n/((2n+1)*factorial(2n-1)) for n = 1:6])
    _coscu(x::Union{Float64,ComplexF64}) =
        fastabs(x) < 0.44 ? x*evalpoly(x^2, (-0.3333333333333333, 0.03333333333333333, -0.0011904761904761906, 2.2045855379188714e-5, -2.505210838544172e-7, 1.9270852604185937e-9)) :
        isinf_real(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
    _coscu(x::Union{Float32,ComplexF32}) =
        fastabs(x) < 0.817f0 ? x*evalpoly(x^2, (-0.333333335f0, 0.033333335f0, -0.0011904762f0, 2.2045855f-5)) :
        isinf_real(x) ? zero(x) : ((s,c)=sincos(x); (x*c-s)/(x^2))
    _coscu(x::Float16) = Float16(_coscu(Float32(x)))
    _coscu(x::ComplexF16) = ComplexF16(_coscu(ComplexF32(x)))
end
