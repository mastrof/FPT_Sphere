export DiffusiveField, source, background, excess, source_size, source_origin

abstract type AbstractSphericalField{D,T} <: AbstractChemoattractant{D,T} end

source(model::ABM) = source(chemoattractant(model))
source(f::AbstractSphericalField) = f.sphere
background(model::ABM) = background(chemoattractant(model))
background(f::AbstractSphericalField) = f.C₀
excess(model::ABM) = excess(chemoattractant(model))
excess(f::AbstractSphericalField) = f.C₁
source_size(model::ABM) = source_size(chemoattractant(model))
source_size(f::AbstractSphericalField) = source(f).radius
source_origin(model::ABM) = source_origin(chemoattractant(model))
source_origin(f::AbstractSphericalField{D,T}) where {D,T} = SVector{D,T}(source(f).center.coords)

@kwdef struct DiffusiveField{D,T} <: AbstractSphericalField{D,T}
    sphere::Meshes.Ball{D,T}
    C₀::T
    C₁::T
    concentration_field::Function = field_diffusive
    concentration_gradient::Function = gradient_diffusive
    concentration_ramp::Function = (pos, model) -> zero(T)
    diffusivity::T = T(608.0) # μm²/s
end
function DiffusiveField(C₀::Real, C₁::Real, R::Real, origin::SVector{D,T}) where {D,T<:Real}
    @assert C₀ >= zero(T)
    @assert C₁ >= zero(T)
    @assert R > zero(T)
    sphere = Meshes.Ball{D,T}(Meshes.Point{D,T}(origin), T(R))
    return DiffusiveField{D,T}(; sphere, C₀=T(C₀), C₁=T(C₁))
end

function field_diffusive(pos, model)
    C₀ = background(model)
    C₁ = excess(model)
    R = source_size(model)
    P = source_origin(model)
    r = euclidean_distance(pos, P, model)
    return C₀ + C₁*R/r
end

function gradient_diffusive(pos::SVector{D}, model) where D
    C₀ = background(model)
    C₁ = excess(model)
    R = source_size(model)
    P = source_origin(model)
    rvec = distancevector(P, pos, model)
    r² = dot(rvec, rvec)
    r³ = r²*sqrt(r²)
    return SVector{D}(-C₁*R*r/r³ for r in rvec)
end


@kwdef struct ExpField{D,T} <: AbstractSphericalField{D,T}
    sphere::Meshes.Ball{D,T}
    C₀::T
    C₁::T
    λ::T
    concentration_field::Function = field_exp
    concentration_gradient::Function = gradient_exp
    concentration_ramp::Function = (pos, model) -> zero(T)
    diffusivity::T = T(608.0) # μm²/s
end
function ExpField(C₀::Real, C₁::Real, R::Real, origin::SVector{D,T}, λ::Real) where {D,T<:Real}
    @assert C₀ >= zero(T)
    @assert C₁ >= zero(T)
    @assert R > zero(T)
    @assert λ > zero(T)
    sphere = Meshes.Ball{D,T}(Meshes.Point{D,T}(origin), T(R))
    return ExpField{D,T}(; sphere, C₀=T(C₀), C₁=T(C₁), λ=T(λ))
end
decay_length(model::ABM) = decay_length(chemoattractant(model))
decay_length(f::ExpField) = f.λ

function field_exp(pos, model)
    C₀ = background(model)
    C₁ = excess(model)
    R = source_size(model)
    P = source_origin(model)
    λ = decay_length(model)
    r = euclidean_distance(pos, P, model)
    return C₀ + C₁*R*exp(-(r-R)/λ)/r
end

function gradient_exp(pos::SVector{D}, model) where D
    C₀ = background(model)
    C₁ = excess(model)
    R = source_size(model)
    P = source_origin(model)
    λ = decay_length(model)
    rvec = distancevector(P, pos, model)
    r² = dot(rvec, rvec)
    r = sqrt(r²)
    r³ = r * r²
    return SVector{D}(-(λ+r)*exp(-(r-R)/λ)/(λ*r³)*r for r in rvec)
end
