@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the
    # size of the precompile file and potentially make loading faster.
    Symbolics.@variables w x y z a b c d e
    s = randn(Float64)
    v = randn(QuatVecF64)
    r = randn(RotorF64)
    q = randn(QuaternionF64)
    𝓈 = w
    𝓋 = QuatVec(x, y, z)
    𝓇 = Rotor(a, b, c, d)
    𝓆 = Quaternion(w, x, y, z)

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether they belong to
        # your package or not (on Julia 1.8 and higher)
        r(v)
        Symbolics.simplify.(𝓇(𝓋))
        for a ∈ [s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆]
            conj(a)
            for b ∈ [s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆]
                a * b
                a / b
                a + b
                a - b
            end
        end

    end
end
