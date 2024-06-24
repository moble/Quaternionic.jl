@testset verbose=true "Base" begin
    @testset "Numbers $T" for T in Types
        # Note that, because `Symbolics.Num` is a weird type, we have to be a little more
        # explicit below than we normally would be.  Also, because of signed zeros in the
        # float types, we have to take the absolute value of the difference before comparing
        # to zero.

        # Define basis elements
        u = Quaternion{T}(1)
        i = Quaternion{T}(𝐢)
        j = Quaternion{T}(𝐣)
        k = Quaternion{T}(𝐤)
        basis = [u, i, j, k]

        # Check basis elements
        for (index, element) in enumerate(basis)
            for component in 1:4
                if component == index
                    @test element[component] == one(T)
                else
                    @test element[component] == zero(T)
                end
            end
        end
        @test u[1]==one(T) && u[2]==zero(T) && u[3]==zero(T) && u[4]==zero(T)
        @test i[1]==zero(T) && i[2]==one(T) && i[3]==zero(T) && i[4]==zero(T)
        @test j[1]==zero(T) && j[2]==zero(T) && j[3]==one(T) && j[4]==zero(T)
        @test k[1]==zero(T) && k[2]==zero(T) && k[3]==zero(T) && k[4]==one(T)

        # Check equality with constants; note that these are *equal*, but not *identical*
        @test u == one(T) + zero(T)*𝐢 == one(T)
        @test i == 𝐢 == imx
        @test j == 𝐣 == imy
        @test k == 𝐤 == imz

        @test Quaternion{T}(1, 2, 3, 4) == Quaternion{T}(SVector{4, T}(1, 2, 3, 4))
        @test Quaternion{T}(0, 2, 3, 4) == Quaternion{T}(2, 3, 4)
        @test Quaternion{T}(1, 0, 0, 0) == Quaternion{T}(1)
        @test quaternion(T[1, 2, 3, 4]...) == quaternion(SVector{4, T}(1, 2, 3, 4))
        @test quaternion(T[0, 2, 3, 4]...) == quaternion(T(2), T(3), T(4))
        @test quaternion(T[1, 0, 0, 0]...) == quaternion(T(1))
        if !(T<:Integer)
            @test rotor(T(1), 2, 3, 4) == Rotor{T}(SVector{4, T}(1, 2, 3, 4)/√T(30))
            @test Rotor{T}(1, 0, 0, 0) == Rotor{T}(1)
            if !(T<:Symbolics.Num)
                @test rotor(T[1, 2, 3, 4]...) ≈ rotor(SVector{4, T}(1, 2, 3, 4)/√T(30)) rtol=0 atol=2eps(T)
                @test rotor(T[0, 2, 3, 4]...) ≈ rotor(T(2), T(3), T(4)) rtol=0 atol=2eps(T)
                @test rotor(T[1, 0, 0, 0]...) ≈ rotor(T(1)) rtol=0 atol=2eps(T)
            end
        end
        @test rotor(T[1, 0, 0, 0]...) == rotor(T(1))
        @test QuatVec{T}(1, 2, 3, 4) == QuatVec{T}(SVector{4, T}(1, 2, 3, 4))
        @test QuatVec{T}(0, 2, 3, 4) == QuatVec{T}(2, 3, 4)
        @test QuatVec{T}(1, 0, 0, 0) == QuatVec{T}(1)
        @test quatvec(T[0, 2, 3, 4]) == quatvec(SVector{4, T}(0, 2, 3, 4))
        @test quatvec(T[0, 2, 3, 4]) == quatvec(T(2), T(3), T(4))
        @test quatvec(T[0, 0, 0, 0]) == quatvec(T(0))
        for v ∈ basis
            @test !(quatvec(v) == one(T))
            @test !(one(T) == quatvec(v))
            @test isequal(quatvec(v), v)
            @test isequal(v, quatvec(v))
            @test isequal(quatvec(v), quatvec(v))
            for (i_v,v′) ∈ enumerate(basis)
                if v !== v′
                    @test !(quatvec(v) == quatvec(v′))
                    @test !isequal(quatvec(v), quatvec(v′))
                    @test !(quatvec(v) == v′)
                    @test !isequal(quatvec(v), v′)
                    @test !(v == quatvec(v′))
                    @test !isequal(v, quatvec(v′))
                elseif i_v>1
                    @test quatvec(v) == quatvec(v′)
                    @test isequal(quatvec(v), quatvec(v′))
                    @test quatvec(v) == v′
                    @test isequal(quatvec(v), v′)
                    @test v == quatvec(v′)
                    @test isequal(v, quatvec(v′))
                end
            end
        end

        for v ∈ [𝐢, 𝐣, 𝐤]
            @test Symbolics.Num(1)*v != one(T)
            @test !isequal(Symbolics.Num(1)*v, one(T))
            @test Quaternion(Symbolics.Num(one(T))) == one(T)
            @test Quaternion(Symbolics.Num(7one(T))) == 7one(T)
            @test one(T) == Quaternion(Symbolics.Num(one(T)))
            @test 7one(T) == Quaternion(Symbolics.Num(7one(T)))

            @test QuatVec(Symbolics.Num[1,2,3,4]) != one(T)
            @test QuatVec(Symbolics.Num[7,2,3,4]) != 7one(T)
            @test one(T) != QuatVec(Symbolics.Num[1,2,3,4])
            @test 7one(T) != QuatVec(Symbolics.Num[7,2,3,4])

            @test QuatVec(Symbolics.Num[1,2,3,4]) == QuatVec{T}(0,2,3,4)
            @test QuatVec{T}(0,2,3,4) == QuatVec(Symbolics.Num[1,2,3,4])
            @test QuatVec(Symbolics.Num[1,2,3,4]) == Quaternion{T}(0,2,3,4)
            @test Quaternion{T}(0,2,3,4) == QuatVec(Symbolics.Num[1,2,3,4])
            @test QuatVec(Symbolics.Num[1,2,3,4]) != Quaternion{T}(1,2,3,4)
            @test Quaternion{T}(1,2,3,4) != QuatVec(Symbolics.Num[1,2,3,4])
        end

        # Test indexing
        q = quaternion(T(1), T(2), T(3), T(4))
        @test q[[1, 2]] == [T(1), T(2)]
        @test q[[3, 4]] == [T(3), T(4)]
        @test q[[4, 2]] == [T(4), T(2)]
        @test q[[4, 2, 3]] == [T(4), T(2), T(3)]
        @test q[[4, 2, 3, 1]] == [T(4), T(2), T(3), T(1)]

        # Test copy constructor and self-equality
        @test quaternion(u) == Quaternion{T}(u) == u
        @test quaternion(i) == Quaternion{T}(i) == i
        @test quaternion(j) == Quaternion{T}(j) == j
        @test quaternion(k) == Quaternion{T}(k) == k
        @test Quaternionic.wrapper(quaternion(u)) == Quaternion
        @test Quaternionic.wrapper(rotor(u)) == Rotor
        @test Quaternionic.wrapper(quatvec(u)) == QuatVec
        @test u == one(T)
        @test one(T) == u
        @test i != one(T)
        @test one(T) != i
        @test j != one(T)
        @test one(T) != j
        @test k != one(T)
        @test one(T) != k
        @test isequal(u, u)
        @test !isequal(u, i)
        @test !isequal(u, j)
        @test !isequal(u, k)
        @test !isequal(i, u)
        @test isequal(i, i)
        @test !isequal(i, j)
        @test !isequal(i, k)
        @test !isequal(j, u)
        @test !isequal(j, i)
        @test isequal(j, j)
        @test !isequal(j, k)
        @test !isequal(k, u)
        @test !isequal(k, i)
        @test !isequal(k, j)
        @test isequal(k, k)

        # Check isapprox
        @test u ≈ one(T)
        @test one(T) ≈ u
        if T !== Symbolics.Num  # Symbolics.Num(1) ≉ Symbolics.Num(2) doesn't work
            @test i ≉ one(T)
            @test one(T) ≉ i
            @test j ≉ one(T)
            @test one(T) ≉ j
            @test k ≉ one(T)
            @test one(T) ≉ k
        end

        ### Actually, w need to return nonsensical basis elements because of auto-diff packages
        # # Check nonsensical basis elements are not allowed
        # @test_throws DomainError one(QuatVec)
        # @test_throws DomainError one(QuatVec{T})
        # @test_throws DomainError one(QuatVec{T}(1, 2, 3))
        # @test_throws DomainError zero(Rotor)
        # @test_throws DomainError zero(Rotor{T})
        # @test_throws DomainError zero(Rotor{T}(1))
        @test typeof(one(QuatVec{T})) === Quaternion{T}
        @test typeof(zero(Rotor{T})) === Quaternion{T}

        # Check "real" part
        @test real(u) == one(T)
        @test real(i) == zero(T)
        @test real(j) == zero(T)
        @test real(k) == zero(T)

        # Check "imag" part
        @test imag(u) == [zero(T), zero(T), zero(T)]
        @test imag(i) == [one(T), zero(T), zero(T)]
        @test imag(j) == [zero(T), one(T), zero(T)]
        @test imag(k) == [zero(T), zero(T), one(T)]

        # Check "isreal"
        @test isreal(u)
        @test !isreal(i)
        @test !isreal(j)
        @test !isreal(k)

        # Check "isinteger"
        if T != Symbolics.Num
            @test isinteger(u)
            @test !isinteger(1.2u)
        end
        @test !isinteger(i)
        @test !isinteger(j)
        @test !isinteger(k)

        if T<:AbstractFloat
            # Check "isnan"
            @test !isnan(u)
            @test !isnan(i)
            @test !isnan(j)
            @test !isnan(k)
            @test isnan(T(NaN) + 0imx)
            @test isnan(T(NaN)imx)
            @test isnan(T(NaN)imy)
            @test isnan(T(NaN)imz)

            # Check "isinf"
            @test !isinf(u)
            @test !isinf(i)
            @test !isinf(j)
            @test !isinf(k)
            @test isinf(T(Inf) + 0imx)
            @test isinf(T(Inf)imx)
            @test isinf(T(Inf)imy)
            @test isinf(T(Inf)imz)
        end

        # Check "iszero"
        @test iszero(0u)
        @test iszero(0i)
        @test iszero(0j)
        @test iszero(0k)
        @test !iszero(u)
        @test !iszero(i)
        @test !iszero(j)
        @test !iszero(k)

        # Check "isone"
        @test isone(u)
        @test !isone(2u)
        @test !isone(i)
        @test !isone(j)
        @test !isone(k)

        # Check "round"
        if T<:AbstractFloat
            @test eltype(round(T(1.2) + 0imx)) === T
            @test round(T(1.2) + 0imx) == 1 + 0imx
            @test round(T(1.2)imx) == 0 + 1imx
            @test round(T(1.2)imy) == 0 + 1imy
            @test round(T(1.2)imz) == 0 + 1imz
            let q = randn(Quaternion{T})
                @test components(round(q, digits=4)) == round.(components(q), digits=4)
            end
        end

        if T != Symbolics.Num
            # Check "in"
            @test u ∈ 0:2
            @test u ∉ 2:4
            @test i ∉ 0:2
            @test j ∉ 0:2
            @test k ∉ 0:2
        end

        # Check "flipsign"
        @test flipsign(u, 1) == u
        @test flipsign(u, -1) == -u
        @test flipsign(i, 1) == i
        @test flipsign(i, -1) == -i
        @test flipsign(j, 1) == j
        @test flipsign(j, -1) == -j
        @test flipsign(k, 1) == k
        @test flipsign(k, -1) == -k

    end

    @testset "bswap" begin
        @test bswap(quaternion(1)) == quaternion(2^56)
        @test bswap(1imx) == (2^56)imx
        @test bswap(1imy) == (2^56)imy
        @test bswap(1imz) == (2^56)imz
        @test bswap(bswap(quaternion(1))) == quaternion(1)
        @test bswap(bswap(1imx)) == 1imx
        @test bswap(bswap(1imy)) == 1imy
        @test bswap(bswap(1imz)) == 1imz
    end

    @testset "hash" begin
        for T1 in [FloatTypes...; IntTypes...]
            for T2 in [FloatTypes...; IntTypes...]
                q1 = Quaternion{T1}(1, 2, 3, 4)
                q2 = Quaternion{T2}(1, 2, 3, 4)
                @test isequal(q1, q2) && hash(q1)==hash(q2)
            end
        end

        for T1 in FloatTypes
            for T2 in FloatTypes
                q1 = quaternion(T1(0.0))
                q2 = quaternion(T2(-0.0))
                @test !isequal(q1, q2) && hash(q1)!=hash(q2)
                q1 = quaternion(T1(NaN))
                q2 = quaternion(T2(NaN))
                @test isequal(q1, q2) && hash(q1)==hash(q2)
            end
        end
    end

    @testset "io" begin
        io = IOBuffer()

        Base.show(io, MIME("text/plain"), Quaternion{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤"
        Base.show(io, MIME("text/plain"), Quaternion{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "1 + 2𝐢 + 3𝐣 + 4𝐤"
        Base.show(io, MIME("text/plain"), quaternion(a, b, c, d))
        @test String(take!(io)) == "a + b𝐢 + c𝐣 + d𝐤"
        Base.show(io, MIME("text/plain"), quaternion(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "a - b + b*c𝐢 + (c / d)𝐣 + (d + e)𝐤"
        Base.show(io, MIME("text/latex"), Quaternion{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$1.0 + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + 4.0\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), Quaternion{Float64}(1, 2, 3, 4e-9))
        @test String(take!(io)) == "\$1.0 + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + \\left(4.0e-9\\right)\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), Quaternion{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$1 + 2\\,\\mathbf{i} + 3\\,\\mathbf{j} + 4\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), quaternion(a, b, c, d))
        @test String(take!(io)) == "\$a + b\\,\\mathbf{i} + c\\,\\mathbf{j} + d\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), quaternion(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "\$a - b + b c\\,\\mathbf{i} + \\frac{c}{d}\\,\\mathbf{j} + \\left(d + e\\right)\\,\\mathbf{k}\$"

        Base.show(io, MIME("text/plain"), QuatVec{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == " + 2.0𝐢 + 3.0𝐣 + 4.0𝐤"
        Base.show(io, MIME("text/plain"), QuatVec{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == " + 2𝐢 + 3𝐣 + 4𝐤"
        Base.show(io, MIME("text/plain"), quatvec(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == " + b*c𝐢 + (c / d)𝐣 + (d + e)𝐤"
        Base.show(io, MIME("text/latex"), QuatVec{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$ + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + 4.0\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), QuatVec{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$ + 2\\,\\mathbf{i} + 3\\,\\mathbf{j} + 4\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), quatvec(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "\$ + b c\\,\\mathbf{i} + \\frac{c}{d}\\,\\mathbf{j} + \\left(d + e\\right)\\,\\mathbf{k}\$"

        Base.show(io, MIME("text/plain"), rotor(1, 5, 5, 7))
        @test String(take!(io)) == "rotor(0.1 + 0.5𝐢 + 0.5𝐣 + 0.7𝐤)"
        Base.show(io, MIME("text/latex"), rotor(1, 3, 3, 9))
        @test String(take!(io)) == "\$0.1 + 0.3\\,\\mathbf{i} + 0.3\\,\\mathbf{j} + 0.9\\,\\mathbf{k}\$"

        for T in PrimitiveTypes
            io = IOBuffer()
            q = Quaternion{T}(1, 2, 3, 4)
            write(io, q)
            seekstart(io)
            p = read(io, typeof(q))
            @test q == p
        end
    end

    @testset "Differential" begin
        Symbolics.@variables t Q(t)[1:4] R(t)[1:4] V(t)[1:3]
        ∂ₜ = Symbolics.Differential(t)
        Q = quaternion(Q...)
        R = rotor(R...)
        V = quatvec(V...)
        @test ∂ₜ(Q) == quaternion(∂ₜ(Q[1]), ∂ₜ(Q[2]), ∂ₜ(Q[3]), ∂ₜ(Q[4]))
        @test ∂ₜ(R) == quaternion(∂ₜ(R[1]), ∂ₜ(R[2]), ∂ₜ(R[3]), ∂ₜ(R[4]))
        @test ∂ₜ(V) == quatvec(∂ₜ(V[2]), ∂ₜ(V[3]), ∂ₜ(V[4]))
    end

end
