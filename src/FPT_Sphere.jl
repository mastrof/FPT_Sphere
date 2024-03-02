module FPT_Sphere

using DrWatson
using LinearAlgebra
using Agents
using MicrobeAgents
using StaticArrays
import Meshes # Point, Vec, Ball, Ray, intersects

include("utils.jl")
include("concentration_fields.jl")
include("encounter.jl")

end # module
