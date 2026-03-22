# Tests for normalization of complexified quaternions as Lorentz spinors in the STA.
#
# Complexified quaternions ℍ(ℂ) are the natural home for Lorentz-group rotors in the
# Spacetime Algebra (STA) with signature -+++.  Following Doran & Lasenby, the
# GA identification is (with 𝐞_{x,y,z} as the spatial unit vectors):
#
#   𝐢 = 𝐞𝐲 𝐞𝐳  (= yz)
#   𝐣 = 𝐞𝐳 𝐞𝐱  (= -zx, note sign!)
#   𝐤 = 𝐞𝐱 𝐞𝐲  (= xy)
#   𝒾 (complex imaginary) ↔ pseudoscalar 𝐈 = 𝐞𝐭 𝐞𝐱 𝐞𝐲 𝐞𝐳  (= txyz)
#
# From these it follows that the spacetime bivectors map as:
#
#   tx ↔ -im*𝐢   (boost in x)
#   ty ↔  im*𝐣   (boost in y)
#   tz ↔ -im*𝐤   (boost in z)
#
# The reverse in GA flips sign on all bivector parts.  In the even subalgebra, mapped
# back to complexified quaternions, this is exactly the quaternion conjugate:
# conj(w + x𝐢 + y𝐣 + z𝐤) = w - x𝐢 - y𝐣 - z𝐤, *without* conjugating the complex
# coefficients w, x, y, z ∈ ℂ.
#
# A rotor R is normalized by R R̃ = 1, which in complexified-quaternion terms means
# q * conj(q) = 1, where the norm squared is w² + x² + y² + z² ∈ ℂ (complex squares,
# NOT |w|² + |x|² + |y|² + |z|²).
#
# This is why the specialised _hypot for Complex components is needed: the standard
# hypot uses abs (absolute value), giving the Euclidean ℂ⁴ norm instead of the spinor
# norm.

@testset "Lorentz/STA normalization" begin

    # Use Float64 and Float32 only; Float16 lacks precision for cosh/sinh near the
    # cancellation cosh² - sinh² = 1, and BigFloat is tested in algebra.jl.
    LorentzTypes = [Float64, Float32]

    # ────────────────────────────────────────────────────────────────────────────────
    # 1.  The spinor norm differs from the Euclidean norm on ℂ⁴
    #
    # Use a scaled boost rotor: v = λ·(cosh(φ/2), -im·sinh(φ/2), 0, 0), λ ∈ ℝ.
    #   Spinor norm² = λ²(cosh² - sinh²) = λ²  →  _hypot = λ  (real, positive)
    #   Euclidean norm = λ√(cosh² + sinh²) > λ  for φ ≠ 0
    #
    # After spinor-normalisation we recover the original unit boost rotor.
    # With the Euclidean norm instead, the imaginary component would be too small.
    # ────────────────────────────────────────────────────────────────────────────────
    @testset "_hypot is spinor norm, not Euclidean norm, for $T" for T in LorentzTypes
        ϵ = 32eps(T)

        for φ ∈ T[0.5, 1.0, 1.5, 2.0]
            ch, sh = cosh(φ/2), sinh(φ/2)

            for λ ∈ T[0.5, 1, 2, 5]
                v = SVector{4,Complex{T}}(λ*ch, -im*λ*sh, 0, 0)

                # _hypot should return λ (real), not the Euclidean λ√(cosh²+sinh²)
                h = Quaternionic._hypot(Tuple(v))
                @test h ≈ Complex{T}(λ) rtol=ϵ
                @test isapprox(imag(h), zero(T); atol=ϵ*λ)

                # Euclidean norm is strictly larger
                @test λ * √(ch^2 + sh^2) > λ + ϵ

                # Spinor-normalisation recovers the unit boost rotor
                normalized = v / h
                @test normalized[1] ≈ Complex{T}(ch)    rtol=ϵ
                @test normalized[2] ≈ Complex{T}(-im*sh) rtol=ϵ
                @test normalized[3] ≈ zero(Complex{T})  atol=ϵ
                @test normalized[4] ≈ zero(Complex{T})  atol=ϵ

                # And its spinor norm is 1
                h2 = Quaternionic._hypot(Tuple(normalized))
                @test h2 ≈ one(Complex{T}) rtol=ϵ
            end
        end
    end

    # ────────────────────────────────────────────────────────────────────────────────
    # 2.  Spatial rotation rotors: all-real components
    #
    # A rotation by θ about axis n̂ = (nx, ny, nz):
    #   R = cos(θ/2) + sin(θ/2) (nx𝐢 + ny𝐣 + nz𝐤)
    #
    # With real components, the spinor norm² = cos²(θ/2) + sin²(θ/2) = 1, so these are
    # already unit.  When passed as Complex{T}, rotor() should return the same values.
    # ────────────────────────────────────────────────────────────────────────────────
    @testset "Spatial rotation rotors, $T" for T in LorentzTypes
        ϵ = 32eps(T)

        for θ ∈ T[0, π/7, π/4, π/3, π/2, 2π/3, π, 4π/3, 7π/4, 2π]
            for (nx, ny, nz) ∈ [
                (T(1), T(0), T(0)),
                (T(0), T(1), T(0)),
                (T(0), T(0), T(1)),
                (T(1)/√T(2), T(1)/√T(2), T(0)),
                (T(1)/√T(3), T(1)/√T(3), T(1)/√T(3)),
            ]
                c, s = cos(θ/2), sin(θ/2)
                q = rotor(Complex{T}(c), Complex{T}(nx*s), Complex{T}(ny*s), Complex{T}(nz*s))
                r = q * conj(q)
                @test r[1] ≈ one(Complex{T})       rtol=ϵ
                @test r[2] ≈ zero(Complex{T})      atol=ϵ
                @test r[3] ≈ zero(Complex{T})      atol=ϵ
                @test r[4] ≈ zero(Complex{T})      atol=ϵ
            end
        end
    end

    # ────────────────────────────────────────────────────────────────────────────────
    # 3.  Boost rotors: imaginary quaternion components
    #
    # A Lorentz boost in the j-direction with rapidity φ:
    #   R = exp(φ/2 · Bⱼ)
    # where Bⱼ is the corresponding spacetime bivector.  Using the STA mapping:
    #   R_x = cosh(φ/2) - im·sinh(φ/2)·𝐢    [tx ↔ -im𝐢]
    #   R_y = cosh(φ/2) + im·sinh(φ/2)·𝐣    [ty ↔  im𝐣]
    #   R_z = cosh(φ/2) - im·sinh(φ/2)·𝐤    [tz ↔ -im𝐤]
    #
    # The spinor norm² for each is cosh²(φ/2) + (-im)²sinh²(φ/2) = cosh² - sinh² = 1.
    # ────────────────────────────────────────────────────────────────────────────────
    @testset "Boost rotors, $T" for T in LorentzTypes
        ϵ = 32eps(T)

        for φ ∈ T[0, 0.5, 1.0, 1.5, 2.0]
            ch, sh = cosh(φ/2), sinh(φ/2)

            # Boost in x: (cosh, -im·sinh, 0, 0)
            q_x = rotor(Complex{T}(ch), -im*T(sh), zero(Complex{T}), zero(Complex{T}))
            r = q_x * conj(q_x)
            @test r[1] ≈ one(Complex{T})  rtol=ϵ
            @test r[2] ≈ zero(Complex{T}) atol=ϵ*ch
            @test r[3] ≈ zero(Complex{T}) atol=ϵ
            @test r[4] ≈ zero(Complex{T}) atol=ϵ

            # Boost in y: (cosh, 0, im·sinh, 0)
            q_y = rotor(Complex{T}(ch), zero(Complex{T}), im*T(sh), zero(Complex{T}))
            r = q_y * conj(q_y)
            @test r[1] ≈ one(Complex{T})  rtol=ϵ
            @test r[2] ≈ zero(Complex{T}) atol=ϵ
            @test r[3] ≈ zero(Complex{T}) atol=ϵ*ch
            @test r[4] ≈ zero(Complex{T}) atol=ϵ

            # Boost in z: (cosh, 0, 0, -im·sinh)
            q_z = rotor(Complex{T}(ch), zero(Complex{T}), zero(Complex{T}), -im*T(sh))
            r = q_z * conj(q_z)
            @test r[1] ≈ one(Complex{T})  rtol=ϵ
            @test r[2] ≈ zero(Complex{T}) atol=ϵ
            @test r[3] ≈ zero(Complex{T}) atol=ϵ
            @test r[4] ≈ zero(Complex{T}) atol=ϵ*ch
        end
    end

    # ────────────────────────────────────────────────────────────────────────────────
    # 4.  rotor() normalises by the spinor norm
    #
    # If we scale a (physically unit) rotor by an arbitrary complex factor λ, then call
    # rotor(), the result should recover the original (up to an overall sign/phase that
    # still leaves it a valid unit rotor, i.e. gives q*conj(q) = 1).
    # ────────────────────────────────────────────────────────────────────────────────
    @testset "rotor() normalises by spinor norm, $T" for T in LorentzTypes
        ϵ = 32eps(T)

        # Base unit rotors to scale
        φ, θ = T(0.9), T(π/5)
        ch, sh = cosh(φ/2), sinh(φ/2)
        boost_x_components = (Complex{T}(ch), -im*T(sh), zero(Complex{T}), zero(Complex{T}))
        rot_z_components   = (Complex{T}(cos(θ/2)), zero(Complex{T}), zero(Complex{T}), Complex{T}(sin(θ/2)))

        for λ ∈ Complex{T}[2, 1+im, 3+4im, -2im]
            # Scaling a unit rotor by λ gives spinor norm λ (not |λ|).
            # After rotor() normalises, q*conj(q) should be 1 again.
            for (w, x, y, z) ∈ [boost_x_components, rot_z_components]
                q = rotor(λ*w, λ*x, λ*y, λ*z)
                r = q * conj(q)
                @test r[1] ≈ one(Complex{T})  rtol=ϵ
                @test r[2] ≈ zero(Complex{T}) atol=ϵ
                @test r[3] ≈ zero(Complex{T}) atol=ϵ
                @test r[4] ≈ zero(Complex{T}) atol=ϵ
            end
        end
    end

    # ────────────────────────────────────────────────────────────────────────────────
    # 4b. Pure complex-phase factors exp[Iφ] — the tightest test
    #
    # I = txyz (grade 4) reverses to itself: Ĩ = +I.  Therefore exp[Iφ]·R₀ satisfies
    #   (exp[Iφ]·R₀)~ = exp[Iφ]·R̃₀   →   R R̃ = exp[2Iφ] ≠ 1  for φ ∉ πℤ
    #
    # In ℍ(ℂ) this is a pure complex scalar factor e^{imφ}, with |e^{imφ}| = 1.
    # For a rotation rotor R₀ (real components, Euclidean-unit):
    #
    #   Euclidean norm of e^{imφ}·R₀  = |e^{imφ}| · ‖R₀‖_euc = 1  →  no-op (wrong)
    #   Spinor norm of e^{imφ}·R₀     = e^{imφ}                →  divide out (correct)
    #
    # This is the sharpest possible distinction between the two norms.
    # ────────────────────────────────────────────────────────────────────────────────
    @testset "Pure phase exp[Iφ] is not a Lorentz transformation, $T" for T in LorentzTypes
        ϵ = 32eps(T)

        θ = T(π/5)
        # Pure rotation rotor: real components, Euclidean-unit, spinor-unit
        rot_z = (Complex{T}(cos(θ/2)), zero(Complex{T}), zero(Complex{T}), Complex{T}(sin(θ/2)))

        for φ ∈ T[π/6, π/4, π/3, π/2, 2π/3]
            λ = cos(φ) + im*sin(φ)          # e^{imφ}, |λ| = 1

            # Scaled element is NOT a valid Lorentz rotor: its conj-norm is exp[2iφ] ≠ 1
            (w, x, y, z) = rot_z
            q_scaled = Rotor{Complex{T}}(λ*w, λ*x, λ*y, λ*z)
            r_scaled = q_scaled * conj(q_scaled)
            @test !isapprox(r_scaled[1], one(Complex{T}); rtol=ϵ)  # spinor norm ≠ 1
            @test r_scaled[1] ≈ Complex{T}(cos(2φ) + im*sin(2φ)) rtol=ϵ

            # rotor() normalises and recovers a valid rotor
            q = rotor(λ*w, λ*x, λ*y, λ*z)
            r = q * conj(q)
            @test r[1] ≈ one(Complex{T})  rtol=ϵ
            @test r[2] ≈ zero(Complex{T}) atol=ϵ
            @test r[3] ≈ zero(Complex{T}) atol=ϵ
            @test r[4] ≈ zero(Complex{T}) atol=ϵ
        end
    end

    # ────────────────────────────────────────────────────────────────────────────────
    # 5.  The spinor norm is multiplicative: ‖R₁ R₂‖ = ‖R₁‖ · ‖R₂‖
    #
    # Equivalently, the composition of two unit rotors is a unit rotor.
    # Proof: (R₁R₂)(R₁R₂)̃ = R₁R₂R̃₂R̃₁ = R₁(1)R̃₁ = 1.
    #
    # We test rotation ∘ boost and boost ∘ rotation for several angles/rapidities.
    # ────────────────────────────────────────────────────────────────────────────────
    @testset "Spinor norm is multiplicative (rotation ∘ boost), $T" for T in LorentzTypes
        ϵ = 32eps(T)

        for θ ∈ T[π/6, π/4, π/3, 2π/3]
            for φ ∈ T[0.4, 0.9, 1.5]
                ch, sh = cosh(φ/2), sinh(φ/2)

                # Spatial rotation about z
                q_rot = rotor(
                    Complex{T}(cos(θ/2)), zero(Complex{T}),
                    zero(Complex{T}), Complex{T}(sin(θ/2))
                )
                # Boost in x
                q_boost = rotor(
                    Complex{T}(ch), -im*T(sh), zero(Complex{T}), zero(Complex{T})
                )

                # Both orderings
                for q_prod ∈ [q_rot * q_boost, q_boost * q_rot]
                    r = q_prod * conj(q_prod)
                    @test r[1] ≈ one(Complex{T})  rtol=ϵ
                    @test r[2] ≈ zero(Complex{T}) atol=ϵ
                    @test r[3] ≈ zero(Complex{T}) atol=ϵ
                    @test r[4] ≈ zero(Complex{T}) atol=ϵ
                end
            end
        end
    end

    @testset "Spinor norm is multiplicative (boost ∘ boost), $T" for T in LorentzTypes
        ϵ = 32eps(T)

        for φ₁ ∈ T[0.4, 1.0, 1.6]
            for φ₂ ∈ T[0.3, 0.8, 1.4]
                q1 = rotor(
                    Complex{T}(cosh(φ₁/2)), -im*T(sinh(φ₁/2)),
                    zero(Complex{T}), zero(Complex{T})
                )
                q2 = rotor(
                    Complex{T}(cosh(φ₂/2)), zero(Complex{T}),
                    im*T(sinh(φ₂/2)), zero(Complex{T})
                )
                q_prod = q1 * q2
                r = q_prod * conj(q_prod)
                @test r[1] ≈ one(Complex{T})  rtol=ϵ
                @test r[2] ≈ zero(Complex{T}) atol=ϵ
                @test r[3] ≈ zero(Complex{T}) atol=ϵ
                @test r[4] ≈ zero(Complex{T}) atol=ϵ
            end
        end
    end

end  # @testset "Lorentz/STA normalization"
