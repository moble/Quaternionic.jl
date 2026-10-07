@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the
    # size of the precompile file and potentially make loading faster.
    s = randn(Float64)
    v = randn(QuatVecF64)
    r = randn(RotorF64)
    r2 = randn(RotorF64)
    q = randn(QuaternionF64)
    Rs = [r, r2, r * r2, r2 * r]
    ts = [0.0, 1.0, 2.0, 3.0]

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether they belong to
        # your package or not
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

        # The first call of each of these would otherwise take tens of milliseconds (or, for
        # `squad`, about a second) to compile.
        exp(v)
        exp(q)
        log(q)
        log(r)
        sqrt(q)
        sqrt(r)
        q^2
        r^2
        q^0.3
        r^0.3
        slerp(r, r2, 0.3)
        distance(r, r2)
        distance2(r, r2)
        distance(q, q)
        to_rotation_matrix(r)
        from_rotation_matrix(to_rotation_matrix(r))
        squad(Rs, ts, [0.5, 1.5, 2.5])
    end
end
