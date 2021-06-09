"""
    random_rotors()
    random_rotors(17)
    random_rotors((17, 3, 8))
    random_rotors((17, 3, 8), normalize=false)
    random_rotors(T=Quaternion{Float16})


"""
function random_rotors(
    dims::Tuple{Vararg{Int64, N}} where N=();
    normalize::Bool=true,
    T::Type{Quaternion{S}}=Quaternion{Float64}
) where {S}
    q = randn(T, dims...)
    if normalize
        return @. q / abs(q)
    end
    q
end


function Base.randn(rng::AbstractRNG, ::Type{Quaternion{T}}, dims::Dims                     ) where {T}
    A = Array{T}(undef, (4, dims...))
    randn!(rng, A)
    collect(as_quat_array(A))
end
# Note that this method explicitly does not define randn(rng, T), in order to prevent an infinite recursion.
function Base.randn(rng::AbstractRNG, ::Type{Quaternion{T}}, dim1::Integer, dims::Integer...) where {T}
    A = Array{T}(undef, 4, dim1, dims...)
    randn!(rng, A)
    collect(as_quat_array(A))
end
function Base.randn(                  ::Type{Quaternion{T}}, dims::Dims                     ) where {T}
    collect(as_quat_array(randn(default_rng(), T, (4, dims...))))
end
function Base.randn(                  ::Type{Quaternion{T}}, dims::Integer...               ) where {T}
    collect(as_quat_array(randn(default_rng(), T, (4, dims...)...)))
end
