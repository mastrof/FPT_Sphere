export AbstractConcentrationField, DiffusiveField, ExpField,
    field, background, peak, source_size, origin, gradient, time_derivative

abstract type AbstractConcentrationField{D} end

field(model::ABM) = model.ConcentrationField
background(f::AbstractConcentrationField) = f.C₀
peak(f::AbstractConcentrationField) = f.C₁
source_size(f::AbstractConcentrationField) = f.R
origin(f::AbstractConcentrationField) = f.origin
field(f::AbstractConcentrationField) = f.concentration_field
gradient(f::AbstractConcentrationField) = f.concentration_gradient
time_derivative(f::AbstractConcentrationField) = f.concentration_time_derivative


#== diffusive concentration field ==#
@kwdef struct DiffusiveField{D} <: AbstractConcentrationField{D}
    C₀::Float64
    C₁::Float64
    R::Float64
    origin::SVector{D,Float64} = zero(SVector{D,Float64})
    concentration_field::Function = field_diffusive
    concentration_gradient::Function = gradient_diffusive
    concentration_time_derivative::Function = (_,_) -> 0.0
end
function DiffusiveField(C₀::Real, C₁::Real, R::Real, origin::SVector{D,<:Real}) where D
    @assert C₀ >= 0
    @assert C₁ >= 0
    @assert R > 0
    DiffusiveField{D}(; C₀, C₁, R, origin)
end

function field_diffusive(pos, model)
    F = field(model)
    C₀ = background(F)
    C₁ = peak(F)
    R = source_size(F)
    origin = origin(F)
    r = euclidean_distance(pos, origin, model)
    C₀ + C₁*R/r
end

function gradient_diffusive(pos::SVector{D}, model) where D
    F = field(model)
    C₀ = background(F)
    C₁ = peak(F)
    R = source_size(F)
    origin = origin(F)
    rvec = distancevector(origin, pos, model)
    r² = dot(rvec, rvec)
    r³ = r²*sqrt(r²)
    SVector{D}(-C₁*R*r/r³ for r in rvec)
end


#== exponential concentration field ==#
@kwdef struct ExpField{D} <: AbstractConcentrationField{D}
    C₀::Float64
    C₁::Float64
    R::Float64
    origin::SVector{D,Float64} = zero(SVector{D,Float64})
    λ::Float64
    concentration_field::Function = field_exp
    concentration_gradient::Function = gradient_exp
    concentration_time_derivative::Function = (_,_) -> 0.0
end
function ExpField(C₀::Real, C₁::Real, R::Real, origin::SVector{D,<:Real}, λ::Real) where D
    @assert C₀ >= 0
    @assert C₁ >= 0
    @assert R > 0
    @assert λ > 0
    ExpField{D}(; C₀, C₁, R, origin, λ)
end
decay_length(F::ExpField) = F.λ

function field_exp(pos, model)
    F = field(model)
    C₀ = background(F)
    C₁ = peak(F)
    R = source_size(F)
    origin = origin(F)
    λ = decay_length(F)
    r = euclidean_distance(pos, origin, model)
    C₀ + C₁*R*exp(-(r-R)/λ)/r
end

function gradient_exp(pos::SVector{D}, model) where D
    F = field(model)
    C₀ = background(F)
    C₁ = peak(F)
    R = source_size(F)
    origin = origin(F)
    λ = decay_length(F)
    rvec = distancevector(origin, pos, model)
    r² = dot(rvec, rvec)
    r = sqrt(r²)
    r³ = r * r²
    SVector{D}(-(λ+r)*exp(-(r-R)/λ)/(λ*r³)*r for r in rvec)
end
