
function Base.randn(r::AbstractRNG, T::Type{Quaternion{S}}, dims::Tuple{Vararg{Int64, N}} where N) where {S}
    random_array = randn(r, S, (4, dims...))
    as_quat_array(random_array)
end
function Base.randn(r::AbstractRNG, T::Type{Quaternion{S}}) where {S}
    Quaternion(randn(r, S, (4,)))
end
Base.randn(T::Quaternion{S}, dims::Tuple{Vararg{Int64, N}} where N) where {S} = Base.randn(default_rng(), T, dims...)
Base.randn(T::Quaternion{S}) where {S} = Base.randn(default_rng(), T)

function random_rotors(T::Type{Quaternion{S}}, normalize::Bool=true, dims::Tuple{Vararg{Int64, N}} where N=()) where {S}
    q = randn(T, dims...)
    if normalize
        return q ./ abs.(q)
    end
    q
end
random_rotors(normalize::Bool, dims::Tuple{Vararg{Int64, N}} where N) = random_rotors(Quaternion{Float64}, normalize, dims)
random_rotors(dims::Tuple{Vararg{Int64, N}} where N) = random_rotors(Quaternion{Float64}, true, dims)
function random_rotors(T::Type{Quaternion{S}}, normalize::Bool=true) where {S}
    q = randn(T)
    if normalize
        return q / abs(q)
    end
    q
end
random_rotors(normalize::Bool=true) = random_rotors(Quaternion{Float64}, normalize)
