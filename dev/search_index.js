var documenterSearchIndex = {"docs":
[{"location":"#Quaternionic.jl","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"","category":"section"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"Quaternions for Julia","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"The goal of this package is to provide a simple but flexible and complete implementation of quaternions, without restricting the interpretation of quaternions to being rotations, but also providing extensive support for rotations.","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"There are numerous ways to construct a Quaternion — the simplest being to just give the components:","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> using Quaternionic\n\njulia> q = Quaternion(1.0, 2.0, 3.0, 4.0)\n1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤\njulia> p = Quaternion(4, 3, 2, 1)\n4 + 3𝐢 + 2𝐣 + 1𝐤","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"Each quaternion type is parametrized by the types of its components (which are promoted to be all the same type).  Any subtype of Real is allowed, and is detected automatically.  For example, q has type Quaternion{Float64}, while p has type Quaternion{Int64}.[1] The base type may be given explicitly if desired, to override the detected type:","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> r = Quaternion{Float64}(4, 3, 2, 1)\n4.0 + 3.0𝐢 + 2.0𝐣 + 1.0𝐤","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"The various Float and Int types work well, as do BigFloat, and the Num type from Symbolics.jl. In particular, we can use symbolic expressions as components:","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> using Quaternionic, Symbolics\n\njulia> @variables a b c d e;\n\njulia> Quaternion(a-b, b*c, c/d, d+e)\na - b + b*c𝐢 + {c*(d^-1)}𝐣 + {d + e}𝐤","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"It is also possible to construct random quaternions using randn with a Quaternion type. In analogy with the complex types, the aliases QuaternionF64, QuaternionF32, and QuaternionF16 are provided, as well as the constants imx, imy, and imz, and (for copy-paste convenience) the aliases 𝐢, 𝐣, and 𝐤 (as Unicode bold character):","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> QuaternionF64\nQuaternionF64 (alias for Quaternion{Float64})\njulia> 0.1 + 2.3imx + 4.5imz\n0.1 + 2.3𝐢 + 0.0𝐣 + 4.5𝐤\njulia> 0.1 + 2.3𝐢 + 0.0𝐣 + 4.5𝐤\n0.1 + 2.3𝐢 + 0.0𝐣 + 4.5𝐤","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"As with the complex im, the result of multiplying imx, etc., with any real number will be a quaternion with the type of the other number.","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"[1]: Note that, mathematically speaking, quaternions can only be defined over a field, which necessarily cannot be an integer type (because the multiplicative inverse of an integer is not generally an integer).  Nonetheless, it is possible to define a Quaternion{<:Integer}, which should behave as expected.  However, many functions (such as exp, log, etc.) will then return a Quaternion of some different type.","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"Components of a quaternion can be accessed as fields:","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> q.w, q.x, q.y, q.z\n(1.0, 2.0, 3.0, 4.0)","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"You can also extract the \"vector\" component (the last three elements) as","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> q.vec\n3-element Vector{Float64}:\n 2.0\n 3.0\n 4.0","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"For convenience, the scalar and vector components can also be accessed in analogy with complex numbers as","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> q.re\n1.0\njulia> q.im\n3-element Vector{Float64}:\n 2.0\n 3.0\n 4.0","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"The basic algebraic operations work as you would expect:","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"julia> p + q\n5.0 + 5.0𝐢 + 5.0𝐣 + 5.0𝐤\njulia> p - q\n3.0 + 1.0𝐢 - 1.0𝐣 - 3.0𝐤\njulia> p * q\n-12.0 + 16.0𝐢 + 4.0𝐣 + 22.0𝐤\njulia> q * p  # Note the non-commutativity\n-12.0 + 6.0𝐢 + 24.0𝐣 + 12.0𝐤\njulia> q / p\n0.6666666666666666 + 0.3333333333333333𝐢 + 0.0𝐣 + 0.6666666666666666𝐤","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"Several essential mathematical functions are also available, including","category":"page"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"abs\nabs2\nconj\nexp\nlog\nsqrt\nangle","category":"page"},{"location":"#Functions","page":"Quaternionic.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"Modules = [Quaternionic]","category":"page"},{"location":"#Quaternionic.imx","page":"Quaternionic.jl","title":"Quaternionic.imx","text":"imx\n\nThe quaternionic unit associated with rotation about the x axis.\n\nExamples\n\njulia> imx * imx\n-1 + 0𝐢 + 0𝐣 + 0𝐤\njulia> 1.2imx\n0.0 + 1.2𝐢 + 0.0𝐣 + 0.0𝐤\n\n\n\n\n\n","category":"constant"},{"location":"#Quaternionic.imy","page":"Quaternionic.jl","title":"Quaternionic.imy","text":"imy\n\nThe quaternionic unit associated with rotation about the y axis.\n\nExamples\n\njulia> imy * imy\n-1 + 0𝐢 + 0𝐣 + 0𝐤\njulia> 1.2imy\n0.0 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤\n\n\n\n\n\n","category":"constant"},{"location":"#Quaternionic.imz","page":"Quaternionic.jl","title":"Quaternionic.imz","text":"imz\n\nThe quaternionic unit associated with rotation about the z axis.\n\nExamples\n\njulia> imz * imz\n-1 + 0𝐢 + 0𝐣 + 0𝐤\njulia> 1.2imz\n0.0 + 0.0𝐢 + 0.0𝐣 + 1.2𝐤\n\n\n\n\n\n","category":"constant"},{"location":"#Quaternionic.Quaternion","page":"Quaternionic.jl","title":"Quaternionic.Quaternion","text":"Quaternion{T<:Real} <: Number\n\nQuaternionic number type with elements of type T.\n\nQuaternionF16, QuaternionF32 and QuaternionF64 are aliases for Quaternion{Float16}, Quaternion{Float32} and Quaternion{Float64} respectively.\n\nSee also: Quaternion\n\n\n\n\n\n","category":"type"},{"location":"#Quaternionic.Quaternion-NTuple{4, Real}","page":"Quaternionic.jl","title":"Quaternionic.Quaternion","text":"Quaternion(w, x, y, z)\nQuaternion(x, y, z)\nQuaternion(w)\nQuaternion(:z)\nQuaternion{T}(w, x, y, z)\n\nCreates a new quaternion with the given components.  The first argument w is the scalar component, and x, y, and z are the corresponding \"vector\" components.  The type of the returned quaternion will be inferred from the input arguments.  If numeric arguments are missing, they will be set to zero. It is also possible to pass one of the symbols :w, :x, :y, or :z to obtain a unit vector (by default, with eltype Float64) along the corresponding direction.  It is also possible to specify the element type T, by passing the type parameter as usual.\n\nExamples\n\njulia> Quaternion(1, 2, 3, 4)\n1 + 2𝐢 + 3𝐣 + 4𝐤\njulia> Quaternion{Float64}(1, 2, 3, 4)\n1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤\njulia> Quaternion(1.0, 2.0, 3.0, 4.0)\n1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤\njulia> Quaternion(2, 3, 4)\n0 + 2𝐢 + 3𝐣 + 4𝐤\njulia> Quaternion(1)\n1 + 0𝐢 + 0𝐣 + 0𝐤\njulia> Quaternion(:z)\n0.0 + 0.0𝐢 + 0.0𝐣 + 1.0𝐤\n\n\n\n\n\n\n","category":"method"},{"location":"#Base.abs-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Base.abs","text":"abs(q)\n\nSquare-root of the sum the squares of the components of the quaternion\n\nExamples\n\njulia> abs(Quaternion(1,2,4,10))\n11.0\n\n\n\n\n\n","category":"method"},{"location":"#Base.abs2-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Base.abs2","text":"abs2(q)\n\nSum the squares of the components of the quaternion\n\nExamples\n\njulia> abs2(Quaternion(1,2,4,10))\n121\n\n\n\n\n\n","category":"method"},{"location":"#Base.angle-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Base.angle","text":"angle(q)\n\nPhase angle in radians of the rotation represented by this quaternion.\n\nNote that this may be different from your interpretation of the angle of a complex number in an important way.  Because quaternions act on vectors by conjugation — as in q*v*conj(q) — there are two copies of q involved in that expression; in some sense, a quaternion acts \"twice\".  Therefore, this angle may be twice what you expect from an analogy with complex numbers — dpending on how you interpret the correspondence between complex numbers and quaternions.  Also, while rotations in the complex plane have a natural choice of axis (the positive z direction), that is not the case for quaternions, which means that the sign of this angle is arbitrary, and we always choose it to be positive.\n\nExamples\n\njulia> θ=1.2;\n\njulia> R=exp(θ * Quaternion(:z) / 2);\n\njulia> angle(R)\n1.2\n\n\n\n\n\n\n","category":"method"},{"location":"#Base.conj-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Base.conj","text":"conj(q)\n\nReturn the quaternion conjugate, which flips the sign of each \"vector\" component.\n\nExamples\n\njulia> conj(Quaternion(1,2,3,4))\n1 - 2𝐢 - 3𝐣 - 4𝐤\n\n\n\n\n\n","category":"method"},{"location":"#Base.exp-Union{Tuple{Quaternion{T}}, Tuple{T}} where T","page":"Quaternionic.jl","title":"Base.exp","text":"exp(q)\n\nExponential of a quaternion\n\nExamples\n\njulia> exp(π/4*Quaternion(:x))  # Rotation by π/2 about the x axis\n0.7071067811865476 + 0.7071067811865475𝐢 + 0.0𝐣 + 0.0𝐤\n\n\n\n\n\n","category":"method"},{"location":"#Base.log-Union{Tuple{Quaternion{T}}, Tuple{T}} where T","page":"Quaternionic.jl","title":"Base.log","text":"log(q)\n\nLogarithm of a quaternion.\n\nAs with the usual complex logarithm, the quaternion logarithm has multiple branches, though the quaternion branches are three-dimensional: for any unit \"vector\" quaternion q̂, you could add any integer multiple of 2πq̂ to the result of this function and still get the same result after exponentiating (within numerical accuracy).  This function is the principal logarithm.\n\nThis function has discontinuous (and fairly arbitrary) behavior along the negative real axis: if the \"vector\" components of the quaternion are precisely zero and the scalar component is negative, the returned quaternion will have scalar component log(-q.w), but will also have a z component of π.  The choice of the z direction is arbitrary; the \"vector\" component of the returned quaternion could be π times any unit vector.\n\nNote that this function is not specialized to unit-quaternion inputs, so the scalar component of the returned value will be nonzero unless the input has precisely unit magnitude.\n\nExamples\n\njulia> log(exp(1.2Quaternion(:y)))\n0.0 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤\n\njulia> log(Quaternion(exp(7)))\n7.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤\n\njulia> log(Quaternion(-exp(7)))\n7.0 + 0.0𝐢 + 0.0𝐣 + 3.141592653589793𝐤\n\n\n\n\n\n","category":"method"},{"location":"#Base.randn-Union{Tuple{T}, Tuple{Random.AbstractRNG, Type{Quaternion{T}}}} where T<:AbstractFloat","page":"Quaternionic.jl","title":"Base.randn","text":"randn([rng=GLOBAL_RNG], [T=Quaternion{Float64}], [dims...])\n\nGenerate a normally distributed random quaternion of type T with mean 0 and standard deviation of norm 1.  Optionally generate an array of such quaternions.  This module currently provides an implementation for the types QuaternionF16, QuaternionF32, and QuaternionF64 (the default).  The values are drawn from the spherically symmetric quaternionic normal distribution of variance 1 (corresponding to each component having independent normal distribution with mean zero and variance 1/4).\n\nSee also: randn_rotor\n\nExamples\n\njulia> randn(QuaternionF64)\n0.4336736009756228 - 0.45087190792840853𝐢 - 0.24723937675211696𝐣 - 0.4514571469326208𝐤\njulia> randn(QuaternionF16, 2, 2)\n2×2 Matrix{QuaternionF16}:\n   0.4321 + 1.105𝐢 + 0.2664𝐣 - 0.1359𝐤   0.064 + 0.9263𝐢 - 0.4138𝐣 + 0.05505𝐤\n 0.2512 - 0.2585𝐢 - 0.2803𝐣 - 0.00964𝐤  -0.1256 + 0.1848𝐢 + 0.03607𝐣 - 0.752𝐤\n\n\n\n\n\n","category":"method"},{"location":"#Base.sqrt-Union{Tuple{Quaternion{T}}, Tuple{T}} where T","page":"Quaternionic.jl","title":"Base.sqrt","text":"sqrt(q)\n\nSquare-root of a quaternion.\n\nThe general formula whenever the denominator is nonzero is\n\nsqrtq = fracq + q sqrt2q + 2qw\n\nThis can be proven by expanding q as q.w + q.vec and multiplying the expression above by itself.\n\nWhen the denominator is zero, this function has discontinuous (and fairly arbitrary) behavior, just as with the quaternion log function.  In this case, either all components are zero — in which case the result is simply the zero quaternion — or the \"vector\" components of the quaternion are precisely zero and the scalar component is negative.  If the latter is true, the denominator above will be a pure-imaginary number.  Because the quaternions come with infinitely many elements that square to -1, it is not clear which imaginary should be used, so we arbitrarily choose to set the result proportional to the z quaternion.  The choice of the z direction is arbitrary; the \"vector\" component of the returned quaternion could be in any direction.\n\nExamples\n\njulia> q = Quaternion(1.2, 3.4, 5.6, 7.8);\n\njulia> sqrtq = √q;\n\njulia> sqrtq^2 ≈ q\ntrue\n\njulia> √Quaternion(4)\n2.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤\n\njulia> √Quaternion(-4)\n0.0 + 0.0𝐢 + 0.0𝐣 + 2.0𝐤\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.abs2vec-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Quaternionic.abs2vec","text":"abs2vec(q)\n\nSum the squares of the \"vector\" components of the quaternion\n\nExamples\n\njulia> abs2vec(Quaternion(1,2,3,6))\n49\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.absvec-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Quaternionic.absvec","text":"absvec(q)\n\nSquare-root of the sum of the squares of the \"vector\" components of the quaternion\n\nExamples\n\njulia> absvec(Quaternion(1,2,3,6))\n7.0\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.as_float_array-Union{Tuple{AbstractArray{Quaternion{T}, N} where N}, Tuple{T}} where T<:AbstractFloat","page":"Quaternionic.jl","title":"Quaternionic.as_float_array","text":"as_float_array(A)\n\nView a quaternion array as an array of real numbers\n\nThis function is fast because no data is copied; the returned quantity is just a \"view\" of the original.  The output view will have an extra initial dimension (of size 4), but is otherwise the same shape as the input array.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.as_quat_array-Union{Tuple{AbstractArray{T, N} where N}, Tuple{T}} where T<:Real","page":"Quaternionic.jl","title":"Quaternionic.as_quat_array","text":"as_quat_array(A)\n\nView a real array as an array of quaternions\n\nThe input array must have an initial dimension whose size is divisible by four (or better yet is 4), because successive indices in that last dimension will be considered successive components of the output quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.from_euler_angles-Tuple{Any, Any, Any}","page":"Quaternionic.jl","title":"Quaternionic.from_euler_angles","text":"from_euler_angles(α, β, γ)\n\nImprove your life drastically.\n\nAssumes the Euler angles correspond to the quaternion R via\n\nR = exp(α 𝐤/2) * exp(β 𝐣/2) * exp(γ 𝐤/2)\n\nwhere 𝐣 and 𝐤 rotate about the fixed y and z axes, respectively, so this reprents an initial rotation about the z axis (in the positive sense) through an angle γ, followed by a rotation about the y axis by β, and a final rotation about the z axis by α.  This is equivalent to performing an initial rotation about z by α, followed by a rotation about the rotated y axis by β, followed by a rotation about the twice-rotated z axis by γ.  The angles naturally must be in radians for this to make any sense.\n\nNOTE: Before opening an issue reporting something \"wrong\" with this function, be sure to read all of this page, especially the very last section about opening issues or pull requests.\n\nSee Also\n\nto_euler_angles: Convert quaternion to Euler angles\nto_euler_phases: Convert quaternion to Euler phases\nfrom_euler_phases: Create quaternion from Euler phases\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.from_euler_phases-Tuple{Any, Any, Any}","page":"Quaternionic.jl","title":"Quaternionic.from_euler_phases","text":"from_euler_phases(zₐ, zᵦ, zᵧ)\nfrom_euler_phases(z)\n\nReturn the quaternion corresponding to these Euler phases.\n\nInterpreting the input quaternion as a rotation (though its normalization scales out), we can define the complex Euler phases from the Euler angles (α, β, γ) as\n\nzₐ ≔ exp(i*α)\nzᵦ ≔ exp(i*β)\nzᵧ ≔ exp(i*γ)\n\nThese are more useful geometric quantites than the angles themselves — being involved in computing spherical harmonics and Wigner's 𝔇 matrices — and can be computed from the components of the corresponding quaternion algebraically (without the use of transcendental functions).\n\nParameters\n\nz::Vector{Complex{T}}: complex vector of length 3, representing the complex phases (zₐ, zᵦ, zᵧ) in that order.\n\nReturns\n\nR::Quaternion{T}\n\nSee Also\n\nto_euler_phases: Convert quaternion to Euler phases\nto_euler_angles: Convert quaternion to Euler angles\nfrom_euler_angles: Create quaternion from Euler angles\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.randn_rotor-Union{Tuple{T}, Tuple{Random.AbstractRNG, Type{T}, Tuple{Vararg{Int64, N}} where N}} where T<:AbstractFloat","page":"Quaternionic.jl","title":"Quaternionic.randn_rotor","text":"randn_rotor([rng=GLOBAL_RNG], [T=Quaternion{Float64}], [dims...])\n\nGenerate a normally distributed random quaternion of type T with mean 0 and norm 1. The result is spherically symmetric, and gives rise a truly random rotation.\n\nSee also: randn\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.to_euler_angles-Tuple{Quaternion}","page":"Quaternionic.jl","title":"Quaternionic.to_euler_angles","text":"to_euler_angles(R)\n\nOpen Pandora's Box.\n\nIf somebody is trying to make you use Euler angles, tell them no, and walk away, and go and tell your mum.\n\nYou don't want to use Euler angles.  They are awful.  Stay away.  It's one thing to convert from Euler angles to quaternions; at least you're moving in the right direction.  But to go the other way?!  It's just not right.\n\nAssumes the Euler angles correspond to the quaternion R via\n\nR = exp(α 𝐤/2) * exp(β 𝐣/2) * exp(γ 𝐤/2)\n\nwhere 𝐣 and 𝐤 rotate about the fixed y and z axes, respectively, so this reprents an initial rotation about the z axis (in the positive sense) through an angle γ, followed by a rotation about the y axis by β, and a final rotation about the z axis by α.  This is equivalent to performing an initial rotation about z by α, followed by a rotation about the rotated y axis by β, followed by a rotation about the twice-rotated z axis by γ.  The angles are naturally in radians.\n\nNOTE: Before opening an issue reporting something \"wrong\" with this function, be sure to read all of this page, especially the very last section about opening issues or pull requests.\n\nReturns\n\nαβγ::Vector{T}\n\nRaises\n\nAllHell if you try to actually use Euler angles, when you could have been using quaternions like a sensible person.\n\nSee Also\n\nfrom_euler_angles: Create quaternion from Euler angles\nto_euler_phases: Convert quaternion to Euler phases\nfrom_euler_phases: Create quaternion from Euler phases\n\n\n\n\n\n","category":"method"},{"location":"#Quaternionic.to_euler_phases-Union{Tuple{Quaternion{T}}, Tuple{T}} where T","page":"Quaternionic.jl","title":"Quaternionic.to_euler_phases","text":"to_euler_phases(q)\nto_euler_phases!(z, q)\n\nConvert input quaternion to complex phases of Euler angles\n\nInterpreting the input quaternion as a rotation (though its normalization scales out), we can define the complex Euler phases from the Euler angles (α, β, γ) as\n\nzₐ ≔ exp(i*α)\nzᵦ ≔ exp(i*β)\nzᵧ ≔ exp(i*γ)\n\nThese are more useful geometric quantites than the angles themselves — being involved in computing spherical harmonics and Wigner's 𝔇 matrices — and can be computed from the components of the corresponding quaternion algebraically (without the use of transcendental functions).\n\nReturns\n\nz::Vector{Complex{T}}: complex phases (zₐ, zᵦ, zᵧ) in that order.\n\nSee Also\n\nfrom_euler_phases: Create quaternion from Euler phases\nto_euler_angles: Convert quaternion to Euler angles\nfrom_euler_angles: Create quaternion from Euler angles\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"Quaternionic.jl","title":"Index","text":"","category":"section"},{"location":"","page":"Quaternionic.jl","title":"Quaternionic.jl","text":"Modules = [Quaternionic]","category":"page"}]
}
