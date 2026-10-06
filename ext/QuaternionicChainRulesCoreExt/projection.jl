# Projection of tangents and cotangents onto quaternion arguments.
#
# Conventions:
#
#   * The tangent and cotangent of a `Quaternion{T}` or `Rotor{T}` primal is a
#     `Quaternion{float(T)}`.  A `Rotor`-typed cotangent would be renormalized by `Rotor`
#     arithmetic, so it is never produced.
#   * The tangent and cotangent of a `QuatVec{T}` primal is a `QuatVec{float(T)}`, whose
#     scalar part is zero by construction.
#   * A primal with `Bool` components, such as `𝐢`, `𝐣`, and `𝐤`, is not differentiable, and
#     its cotangent is `NoTangent()`, as for `ProjectTo(::Bool)`.
#   * A primal with integer components has a `Float64` (or other floating-point) cotangent,
#     as for `ProjectTo(::Integer)`.
#   * A quaternion cotangent that reaches a real or complex primal is reduced to its scalar
#     part, projected onto that primal.
#   * For complex components, a pullback applies the conjugate transpose of the Jacobian on
#     ℂ⁴, as ChainRules prescribes; see `qadjoint`, `inner`, and `conjcomponents` in
#     `helpers.jl`.
#
# The methods of `ProjectTo{Quaternion}` and `ProjectTo{QuatVec}` that take a `Tangent`
# exist partly to resolve ambiguities with ChainRulesCore's own methods for `Tangent{<:T}`
# and `Tangent{<:Number}`.

function ProjectTo(q::AbstractQuaternion)
    element = ProjectTo(zero(basetype(q)))
    element isa ProjectTo{NoTangent} && return element
    if q isa QuatVec
        return ProjectTo{QuatVec}(; element=element)
    else
        return ProjectTo{Quaternion}(; element=element)
    end
end

function (p::ProjectTo{Quaternion})(dx::AbstractQuaternion)
    return quaternion(p.element(dx[1]), p.element(dx[2]), p.element(dx[3]), p.element(dx[4]))
end
function (p::ProjectTo{QuatVec})(dx::AbstractQuaternion)
    return quatvec(p.element(dx[2]), p.element(dx[3]), p.element(dx[4]))
end

for P ∈ (:Quaternion, :QuatVec)
    @eval begin
        (p::ProjectTo{$P})(dx::Number) = p(cot(dx))
        (p::ProjectTo{$P})(dx::Tangent{<:Complex}) =
            p(Complex(zerofill(dx.re), zerofill(dx.im)))
        function (p::ProjectTo{$P})(dx::Union{Tangent,Tuple,NamedTuple,AbstractVector})
            Δ = cot(dx)
            return Δ isa AbstractZero ? Δ : p(Δ)
        end
        function (p::ProjectTo{$P})(dx::Tangent{<:$P})
            Δ = cot(dx)
            return Δ isa AbstractZero ? Δ : p(Δ)
        end
        function (p::ProjectTo{$P})(dx::Tangent{<:Number})
            Δ = cot(dx)
            return Δ isa AbstractZero ? Δ : p(Δ)
        end
    end
end

# A quaternion cotangent that reaches a real or complex primal keeps only its scalar part.
# This is not type piracy, because the argument type is ours.  The `Complex` method is needed
# by the complex `log`, whose source mixes complex scalars with quaternions.
(p::ProjectTo{<:Union{Real,Complex}})(dx::AbstractQuaternion) = p(dx[1])
function (p::ProjectTo{<:Union{Real,Complex}})(dx::Tangent{<:AbstractQuaternion})
    Δ = cot(dx)
    return Δ isa AbstractZero ? Δ : p(Δ[1])
end
