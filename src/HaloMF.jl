module HaloMF
using Dierckx
export tinker08MF, tinker10MF
export bocquetMFhy, bocquetMFdm
export psMF, stMF, jenkinsMF
export stBias, stBias1, stBias2, stBias3
export tinker10Bias
include("tinkerMF.jl")
include("bocquetMF.jl")
include("classicMF.jl")
include("stBias.jl")
include("tinkerBias.jl")
end # module
