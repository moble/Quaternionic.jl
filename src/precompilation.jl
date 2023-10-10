@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the
    # size of the precompile file and potentially make loading faster.
    s = randn(Float64)
    v = randn(QuatVecF64)
    r = randn(RotorF64)
    q = randn(QuaternionF64)

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether they belong to
        # your package or not (on Julia 1.8 and higher)
        r(v)
        for a ∈ (s, v, r, q)
            conj(a)
            for b ∈ (s, v, r, q)
                a * b
                a / b
                a + b
                a - b
            end
        end

    end
end
