module HaloMF
using Dierckx
export tinker08MF, tinker10MF
export psMF, stMF, jenkinsMF
include("tinkerMF.jl")
include("classicMF.jl")
end # module
